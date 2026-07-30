# Precompute script for vignettes/BATON-case-study.Rmd
#
# Runs the heavy calibration stages of the case-study vignette and stores
# compact results in case_study_results.rds (same directory). The vignette
# loads that file in its setup chunk and displays the stored output; the
# stage code below is verbatim what the vignette shows in its eval = FALSE
# chunks. To reproduce from scratch:
#
#   cd vignettes/precompute
#   Rscript case-study-precompute.R
#
# On this toy-scale simulator the whole script runs in well under a minute;
# with a realistic patient-level simulator the same stages are hours to
# days, which is why the vignette does not run them at knit time.

library(BATON)
options(BATON.q_warning_shown = TRUE)

elapsed <- list()
tic <- function() assign(".t0", Sys.time(), envir = globalenv())
toc <- function(name) {
  elapsed[[name]] <<- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
  message(sprintf("[%s] %.1f s", name, elapsed[[name]]))
}

# ---------------------------------------------------------------------------
# The evaluator (verbatim copy of the vignette's simulator chunk)
# ---------------------------------------------------------------------------

p_control <- 0.25   # control-arm response rate (null and alternative)
p_treat   <- 0.45   # treatment-arm response rate under the alternative

gs_simulator <- function(theta, fidelity = c("low", "med", "high"),
                         seed = NULL, n_rep = NULL, ...) {
  fidelity <- match.arg(fidelity)
  if (is.null(n_rep)) {
    n_rep <- switch(fidelity, low = 2000, med = 4000, high = 10000)
  }
  if (!is.null(seed)) set.seed(seed)

  # Design parameters (theta arrives fractional; realize an integer design)
  n_arm  <- floor(round(theta$n_max) / 2)            # per-arm maximum
  n1     <- max(2L, round(theta$t_interim * n_arm))  # per-arm at interim
  n2     <- n_arm - n1                               # per-arm stage-2 increment
  t_frac <- n1 / n_arm                               # realized information fraction
  c_interim <- theta$c_eff / sqrt(t_frac)  # O'Brien-Fleming-type interim bound
  c_final   <- theta$c_eff                 # final-analysis critical value
  c_fut     <- theta$c_fut                 # binding interim futility bound

  # Pooled two-sample z statistic for a difference in proportions
  z_stat <- function(x_t, x_c, n) {
    pbar <- (x_t + x_c) / (2 * n)
    se <- sqrt(pmax(pbar * (1 - pbar) * 2 / n, 1e-12))
    ifelse(pbar <= 0 | pbar >= 1, 0, (x_t / n - x_c / n) / se)
  }

  run_hypothesis <- function(p_trt) {
    xt1 <- rbinom(n_rep, n1, p_trt)
    xc1 <- rbinom(n_rep, n1, p_control)
    z1  <- z_stat(xt1, xc1, n1)
    stop_eff <- z1 >= c_interim
    stop_fut <- !stop_eff & (z1 <= c_fut)
    continue <- !stop_eff & !stop_fut
    xt2 <- rbinom(n_rep, n2, p_trt)
    xc2 <- rbinom(n_rep, n2, p_control)
    z2  <- z_stat(xt1 + xt2, xc1 + xc2, n_arm)
    reject  <- stop_eff | (continue & z2 > c_final)
    n_total <- ifelse(continue, 2 * n_arm, 2 * n1)
    list(reject = reject, n_total = n_total)
  }

  alt  <- run_hypothesis(p_treat)    # alternative: treatment works
  null <- run_hypothesis(p_control)  # null: treatment equals control

  metrics <- c(power = mean(alt$reject),
               type1 = mean(null$reject),
               EN    = mean(alt$n_total))
  attr(metrics, "variance") <- c(power = var(alt$reject)  / n_rep,
                                 type1 = var(null$reject) / n_rep,
                                 EN    = var(alt$n_total) / n_rep)
  attr(metrics, "n_rep") <- n_rep
  metrics
}

# ---------------------------------------------------------------------------
# Problem definition for the balanced run (verbatim in the vignette)
# ---------------------------------------------------------------------------

bounds_cal <- list(
  n_max     = c(150, 280),   # widened after the quick run
  t_interim = c(0.25, 0.80),
  c_eff     = c(1.96, 2.40), # floor raised to the single-look critical value
  c_fut     = c(-0.50, 1.00)
)
constraints <- list(power = c("ge", 0.85), type1 = c("le", 0.025))

# ---------------------------------------------------------------------------
# Stage 1: balanced calibration run (with the long-run controls)
# ---------------------------------------------------------------------------

tic()
ckpt_dir <- tempdir()

fit_balanced <- bo_calibrate(
  sim_fun = gs_simulator,
  bounds = bounds_cal,
  objective = "EN",
  constraints = constraints,
  n_init = 40, q = 2, budget = 140, seed = 2026,
  fidelity_levels = c(low = 2000, med = 4000, high = 10000),
  multi_seed_verify = TRUE, multi_seed_n = 5,
  max_walltime_s = 4 * 3600,
  checkpoint_fun = function(snap) {
    saveRDS(snap, file.path(ckpt_dir, "balanced-checkpoint.rds"))
  },
  checkpoint_every = 5,
  progress = FALSE
)
toc("balanced")

# ---------------------------------------------------------------------------
# Stage 2: sensitivity analysis on the balanced fit's surrogates
# ---------------------------------------------------------------------------

tic()
set.seed(303)
sens <- list(
  sobol_EN    = sa_sobol(fit_balanced$surrogates, bounds_cal,
                         outcome = "EN", n_mc = 4000),
  sobol_power = sa_sobol(fit_balanced$surrogates, bounds_cal,
                         outcome = "power", n_mc = 4000),
  grad_EN     = sa_gradients(fit_balanced$surrogates, fit_balanced$best_theta,
                             bounds_cal, outcome = "EN"),
  grad_power  = sa_gradients(fit_balanced$surrogates, fit_balanced$best_theta,
                             bounds_cal, outcome = "power"),
  cov_EN      = cov_effects(fit_balanced$surrogates, bounds_cal,
                            outcome = "EN", n_mc = 500)
)
toc("sensitivity")

# ---------------------------------------------------------------------------
# Stage 3: benchmark BO against random and grid search
# ---------------------------------------------------------------------------

tic()
bench <- benchmark_methods(
  sim_fun = gs_simulator,
  bounds = bounds_cal,
  objective = "EN",
  constraints = constraints,
  strategies = c("bo", "random", "grid"),
  bo_args = list(n_init = 40, q = 2, budget = 140, seeds = 1:3,
                 fidelity_levels = c(low = 2000, med = 4000, high = 10000),
                 slim = TRUE),
  random_args = list(n_samples = 100, seeds = 1:3),
  grid_args = list(resolution = 3),
  progress = FALSE
)
bench_summary <- summarise_benchmark(bench)
toc("benchmark")

# ---------------------------------------------------------------------------
# Stage 4: constraint reliability across calibration seeds
# ---------------------------------------------------------------------------

tic()
rel <- estimate_constraint_reliability(
  sim_fun = gs_simulator,
  bounds = bounds_cal,
  objective = "EN",
  constraints = constraints,
  strategies = c("bo", "random"),
  calibration_seeds = 1:5,
  validation_reps = 200000,
  bo_args = list(n_init = 20, q = 2, budget = 60,
                 fidelity_levels = c(low = 2000, med = 4000, high = 10000),
                 slim = TRUE),
  random_args = list(n_samples = 60),
  progress = FALSE
)
toc("reliability")

# ---------------------------------------------------------------------------
# Stage 5: calibrate three design philosophies with warm-start flow
# ---------------------------------------------------------------------------

tic()
philosophies <- list(
  balanced     = list(objective = "EN",
                      constraints = list(power = c("ge", 0.85),
                                         type1 = c("le", 0.025))),
  conservative = list(objective = "EN",
                      constraints = list(power = c("ge", 0.90),
                                         type1 = c("le", 0.025))),
  aggressive   = list(objective = "EN",
                      constraints = list(power = c("ge", 0.80),
                                         type1 = c("le", 0.025)))
)

phil <- bo_calibrate_philosophies(
  scenario_id = "resp2arm",
  sim_fun = gs_simulator,
  bounds = bounds_cal,
  philosophies = philosophies,
  multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_strict = FALSE,
  n_init = 30, q = 2, budget = 100, seed = 11,
  fidelity_levels = c(low = 2000, med = 4000, high = 10000),
  slim = TRUE, progress = FALSE
)
toc("philosophies")

# Compact per-philosophy comparison table (design + verified OCs)
phil_table <- do.call(rbind, lapply(names(phil$fits), function(nm) {
  f <- phil$fits[[nm]]
  ms <- f$multi_seed_summary
  data.frame(
    philosophy = nm,
    n_max = round(f$best_theta$n_max),
    t_interim = round(f$best_theta$t_interim, 2),
    c_eff = round(f$best_theta$c_eff, 3),
    c_fut = round(f$best_theta$c_fut, 2),
    EN = round(ms$EN_mean, 1),
    power = round(ms$power_mean, 3),
    type1 = round(ms$type1_mean, 4),
    verdict = f$verdict,
    stringsAsFactors = FALSE
  )
}))

# ---------------------------------------------------------------------------
# Assemble compact results and save
# ---------------------------------------------------------------------------

# Keep the rds small: strip model objects (surrogates, diagnostics, per-seed
# verification runs) from the balanced fit and the per-run histories from the
# benchmark. Everything the vignette displays survives; the sensitivity
# tables above were computed from the full surrogates before stripping.
fit_slim <- fit_balanced
fit_slim$surrogates <- NULL
fit_slim$diagnostics <- NULL
fit_slim$multi_seed_runs <- NULL

bench_slim <- bench
bench_slim$results$history <- NULL

cs <- list(
  meta = list(
    created = Sys.time(),
    baton_version = as.character(utils::packageVersion("BATON")),
    elapsed_s = unlist(elapsed)
  ),
  balanced = fit_slim,
  sens = sens,
  bench = bench_slim,
  bench_summary = bench_summary,
  rel = rel,
  phil_manifest = phil$manifest,
  phil_table = phil_table
)

saveRDS(cs, "case_study_results.rds", compress = "xz")
message(sprintf("Saved case_study_results.rds (%.0f KB)",
                file.size("case_study_results.rds") / 1024))
message(sprintf("Total precompute time: %.1f s", sum(unlist(elapsed))))
