# Benchmark cached vs per-candidate GP prediction in acquisition loop
#
# Run from package root:
#   Rscript inst/benchmarks/acquisition_caching_benchmark.R

suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required to load BATON for benchmarking.")
  }
  pkgload::load_all(quiet = TRUE)
})

set.seed(123)

# ---- Setup synthetic problem ----
bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
objective <- "obj"
constraints <- list(power = c("ge", 0.5))
constraint_tbl <- parse_constraints(constraints)
fidelity_levels <- c(low = 100, med = 200, high = 400)
primary_fidelity <- names(fidelity_levels)[1]
n_init <- 60
q <- 8
candidate_pool <- 2000

sim_fun <- function(theta, fidelity, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- theta$x1
  x2 <- theta$x2
  noise_sd <- switch(fidelity, low = 0.05, med = 0.02, high = 0.01, 0.05)
  obj <- (x1 - 0.2)^2 + (x2 - 0.8)^2 + stats::rnorm(1, sd = noise_sd)
  power <- x1 + x2 + stats::rnorm(1, sd = noise_sd)
  metrics <- c(obj = obj, power = power)
  variance <- c(obj = noise_sd^2, power = noise_sd^2)
  attr(metrics, "variance") <- variance
  attr(metrics, "n_rep") <- fidelity_levels[[fidelity]]
  list(metrics = metrics, variance = variance, n_rep = fidelity_levels[[fidelity]])
}

# Build an initial history to fit surrogates
history <- initialise_history()
initial_design <- lhs_design(n_init, bounds, seed = 99)
eval_counter <- 0L
for (theta in initial_design) {
  eval_counter <- eval_counter + 1L
  theta <- coerce_theta_types(theta, NULL)
  history <- record_evaluation(
    history = history,
    sim_fun = sim_fun,
    theta = theta,
    bounds = bounds,
    objective = objective,
    constraint_tbl = constraint_tbl,
    fidelity = primary_fidelity,
    fidelity_levels = fidelity_levels,
    eval_id = eval_counter,
    iter = 0L,
    seed = eval_counter,
    progress = FALSE
  )
}

surrogates <- fit_surrogates(
  history,
  objective,
  constraint_tbl,
  covtype = "matern5_2",
  use_hetgp = FALSE
)

best_feasible_value <- best_feasible_objective(history, objective)
constraint_metric_names <- constraint_tbl$metric

run_cached <- function() {
  unit_candidates <- lhs_candidate_pool(candidate_pool, bounds)
  candidate_predictions <- predict_surrogates(surrogates, unit_candidates)
  candidate_prob_feas <- prob_feasibility_batch(candidate_predictions, constraint_tbl)

  acquisition_scores <- evaluate_acquisition(
    acquisition = "eci",
    unit_candidates = unit_candidates,
    surrogates = surrogates,
    constraint_tbl = constraint_tbl,
    objective = objective,
    best_feasible = best_feasible_value,
    predictions = candidate_predictions,
    prob_feas = candidate_prob_feas
  )

  lipschitz <- estimate_lipschitz(surrogates, objective) * 1.5
  selected_idx <- select_batch_local_penalization(
    candidates = unit_candidates,
    acq_scores = acquisition_scores,
    q = q,
    lipschitz = lipschitz
  )

  # Fidelity routing using cached predictions (matches main loop)
  for (i in seq_len(length(selected_idx))) {
    prob_feas <- candidate_prob_feas[selected_idx[i]]
    obj_pred <- candidate_predictions[[objective]]
    pred_obj_mean <- obj_pred$mean[selected_idx[i]]
    pred_obj_sd <- obj_pred$sd[selected_idx[i]]
    cv_objective <- pred_obj_sd / max(abs(pred_obj_mean), 1e-6)

    cv_constraints <- vapply(constraint_metric_names, function(metric_name) {
      if (metric_name %in% names(candidate_predictions)) {
        pred <- candidate_predictions[[metric_name]]
        pred_sd <- pred$sd[selected_idx[i]]
        pred_mean <- pred$mean[selected_idx[i]]
        pred_sd / max(abs(pred_mean), 1e-6)
      } else {
        0
      }
    }, numeric(1))

    # Simple routing (no side effects; just mimics computation)
    if (prob_feas >= 0.4 && any(cv_constraints > cv_objective, na.rm = TRUE)) {
      cv_objective <- max(cv_constraints, na.rm = TRUE)
    }
    invisible(cv_objective)
  }
}

run_uncached <- function() {
  unit_candidates <- lhs_candidate_pool(candidate_pool, bounds)

  acquisition_scores <- evaluate_acquisition(
    acquisition = "eci",
    unit_candidates = unit_candidates,
    surrogates = surrogates,
    constraint_tbl = constraint_tbl,
    objective = objective,
    best_feasible = best_feasible_value
  )

  lipschitz <- estimate_lipschitz(surrogates, objective) * 1.5
  selected_idx <- select_batch_local_penalization(
    candidates = unit_candidates,
    acq_scores = acquisition_scores,
    q = q,
    lipschitz = lipschitz
  )

  # Per-candidate GP calls (previous behavior)
  for (i in seq_len(length(selected_idx))) {
    chosen_unit <- unit_candidates[[selected_idx[i]]]
    prob_feas <- estimate_candidate_feasibility(
      surrogates = surrogates,
      unit_x = list(chosen_unit),
      constraint_tbl = constraint_tbl
    )

    all_preds <- predict_surrogates(surrogates, list(chosen_unit))
    obj_pred <- all_preds[[objective]]
    cv_objective <- obj_pred$sd[[1]] / max(abs(obj_pred$mean[[1]]), 1e-6)

    cv_constraints <- sapply(constraint_metric_names, function(metric_name) {
      if (metric_name %in% names(all_preds)) {
        pred <- all_preds[[metric_name]]
        pred$sd[[1]] / max(abs(pred$mean[[1]]), 1e-6)
      } else {
        0
      }
    })

    if (prob_feas >= 0.4 && any(cv_constraints > cv_objective, na.rm = TRUE)) {
      cv_objective <- max(cv_constraints, na.rm = TRUE)
    }
    invisible(cv_objective)
  }
}

time_method <- function(fun, n = 5L) {
  times <- numeric(n)
  for (i in seq_len(n)) {
    times[[i]] <- system.time(fun())[["elapsed"]]
  }
  c(median = stats::median(times), mean = mean(times))
}

cached_time <- time_method(run_cached, n = 5L)
uncached_time <- time_method(run_uncached, n = 5L)

cat("Cached predictions (median / mean seconds):", sprintf("%.3f / %.3f\n", cached_time[["median"]], cached_time[["mean"]]))
cat("Per-candidate GP calls (median / mean seconds):", sprintf("%.3f / %.3f\n", uncached_time[["median"]], uncached_time[["mean"]]))
cat("Speedup (uncached / cached):", sprintf("%.2fx\n", uncached_time[["median"]] / cached_time[["median"]]))
