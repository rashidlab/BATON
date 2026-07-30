# Golden-master regression guard for performance refactors.
#
# A fully deterministic small calibration. Behavior-preserving performance work
# (parallelism, matrix scoring, incremental aggregation, batched append) MUST keep
# these values byte-identical. Tasks that intentionally change the search
# trajectory (warm-start, penalization, expected-violation fixes) update the
# reference in the SAME commit, with the change called out in the commit message.

golden_sim_fun <- function(theta, fidelity = c("low", "med", "high"),
                           seed = 123, n_rep = NULL, ...) {
  fidelity <- match.arg(fidelity)
  x <- vapply(theta, as.numeric, numeric(1))
  noise <- switch(fidelity, low = 0.02, med = 0.01, 0)
  en <- (x[["x1"]] - 0.5)^2 + (x[["x2"]] - 0.25)^2 + noise
  power <- 0.9 - 0.6 * (x[["x1"]] - 0.5)^2
  type1 <- 0.08 + 0.3 * (x[["x2"]] - 0.25)^2
  res <- c(power = power, type1 = type1, EN = en, ET = en * 1.5)
  attr(res, "variance") <- c(power = 0.001, type1 = 0.001, EN = 0.002, ET = 0.003)
  res
}

run_golden <- function() {
  set.seed(42)
  bo_calibrate(
    sim_fun = golden_sim_fun,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8), type1 = c("le", 0.12)),
    n_init = 6, q = 2, budget = 16, seed = 42, progress = FALSE
  )
}

test_that("golden calibration is deterministic across repeated runs", {
  a <- run_golden()
  b <- run_golden()
  expect_equal(a$history$objective, b$history$objective)
  expect_equal(a$history$feasible, b$history$feasible)
  expect_equal(unlist(a$best_theta), unlist(b$best_theta))
})

test_that("calibration results are independent of BATON.cores (reproducibility)", {
  skip_on_os("windows")  # mclapply is Unix-only; serial path used on Windows
  skip_if_not_installed("parallel")
  old <- getOption("BATON.cores"); on.exit(options(BATON.cores = old), add = TRUE)
  options(BATON.cores = 1L); serial <- run_golden()
  options(BATON.cores = 4L); parallel <- run_golden()
  # A cloud service must return identical designs regardless of how many cores
  # the host has. Deterministic per-fit seeding (run_seeded) guarantees this.
  expect_equal(serial$history$objective, parallel$history$objective)
  expect_equal(serial$history$feasible, parallel$history$feasible)
  expect_equal(unlist(serial$best_theta), unlist(parallel$best_theta))
})

test_that("golden calibration matches the reference snapshot", {
  # Reference re-baselined for v0.8 defect 1 (adaptive fidelity escalation).
  # The default adaptive method previously ran everything at low fidelity
  # (its value score had no fidelity dependence, so the cheapest level always
  # won); the fixed selector escalates high-value candidates, and this run's
  # fidelity mix is 13 low / 3 high. The high-fidelity picks carry zero
  # injected noise in golden_sim_fun (noise = 0.02/0.01/0 by tier), so the
  # incumbent is now a noise-free measurement: best_obj 0.00500 vs the
  # all-low-tier 0.02015 of the previous baseline (which included the +0.02
  # low-tier offset). The initial design (1-9) is unchanged; divergence
  # begins at index 10, the first escalated evaluation.
  fit <- run_golden()
  expect_equal(length(fit$history$objective), 16L)
  expect_equal(round(fit$best_objective, 8), 0.00500346, tolerance = 1e-8)
  expect_equal(round(unlist(fit$best_theta), 5),
               c(x1 = 0.53432, x2 = 0.18815), tolerance = 1e-5)
  # Full objective-vector fingerprint: the strongest guard against silent drift.
  expect_equal(
    round(fit$history$objective, 6),
    round(c(0.021780, 0.517127, 0.139296, 0.144183, 0.509814, 0.056903,
            0.256992, 0.252274, 0.021528, 0.032876, 0.023776, 0.005003,
            0.080482, 0.052200, 0.021191, 0.020513), 6),
    tolerance = 1e-5
  )
  # Per-evaluation feasibility classification is a core output: guard it so a
  # refactor cannot silently reclassify designs without tripping this test.
  expect_equal(
    fit$history$feasible,
    c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  )
  # Fidelity routing is the subject of the v0.8 defect-1 fix: pin the mix so
  # a regression back to all-low (or an escalation blowup) trips this test.
  expect_equal(
    fit$history$fidelity,
    c("low", "low", "low", "low", "low", "low", "low", "low",
      "low", "high", "low", "high", "low", "high", "low", "low")
  )
})
