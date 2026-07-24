# Tests for multi-seed verification (review item 1.5 / Task A7). The routine
# crashed on the documented named-vector simulator contract (res$metrics on an
# atomic vector), destroying a run AFTER the full budget was spent. It defaults
# ON in bo_calibrate_philosophies, with zero prior coverage.

test_that("multi-seed verification handles a vector-returning simulator", {
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    v <- c(power = 0.9, type1 = 0.02, EN = 40)
    attr(v, "variance") <- c(power = 1e-4, type1 = 1e-5, EN = 1)
    v
  }
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.12)))
  res <- expect_no_error(
    BATON:::run_multi_seed_verification(
      sim_fun = sim, theta = list(x = 0.5), seeds = c(1, 2, 3),
      n_rep = 100, objective = "EN", constraint_tbl = ctbl, progress = FALSE)
  )
  # All three seeds should be classified feasible with the correct metrics.
  expect_equal(nrow(res$per_seed), 3L)
  expect_true(all(res$per_seed$feasible))
  expect_equal(res$per_seed$EN, rep(40, 3))
  expect_equal(res$summary$strict_feas, 1L)
})

test_that("strict gate failure yields MULTI_SEED_WARN when multi_seed_strict = FALSE", {
  skip_if_not_installed("DiceKriging")
  # Verification-phase calls are identifiable: they come AFTER budget
  # exhaustion (budget = 6, so call 7 is the first Stage 4 per-seed
  # re-evaluation). The optimization sees feasible metrics throughout, but
  # every verification seed violates the power constraint, so the strict
  # gate fails; with multi_seed_strict = FALSE that must downgrade the
  # verdict to MULTI_SEED_WARN (not FAIL), still returning the full fit.
  n_calls <- 0
  infeasible_on_verify <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    power <- if (n_calls > 6) 0.5 else 0.95
    v <- c(power = power, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    infeasible_on_verify, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 6, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 3, multi_seed_reps = 200,
    multi_seed_strict = FALSE))
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$verdict, "MULTI_SEED_WARN")
  expect_equal(fit$multi_seed_summary$strict_feas, 0L)
  expect_false(any(fit$multi_seed_runs$feasible))
})

test_that("multi-seed verification also handles the list(metrics=) contract", {
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    list(metrics = c(power = 0.9, type1 = 0.02, EN = 40),
         variance = c(power = 1e-4, type1 = 1e-5, EN = 1))
  }
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.12)))
  res <- BATON:::run_multi_seed_verification(
    sim_fun = sim, theta = list(x = 0.5), seeds = c(1, 2),
    n_rep = 100, objective = "EN", constraint_tbl = ctbl, progress = FALSE)
  expect_true(all(res$per_seed$feasible))
  expect_equal(res$per_seed$EN, rep(40, 2))
})
