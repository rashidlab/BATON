# Tests for early-stopping robustness (review item 1.8 / Task A9).

test_that("early_stop accepts and validates acq_rel", {
  # A well-formed early_stop list including acq_rel must not error.
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = 0.9, EN = sum((x - 0.5)^2) * 10)
    attr(res, "variance") <- c(power = 0.001, EN = 0.5); res
  }
  expect_no_error(
    bo_calibrate(sim_fun = sim, bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
                 objective = "EN", constraints = list(power = c("ge", 0.8)),
                 n_init = 5, q = 1, budget = 8, progress = FALSE,
                 early_stop = list(enabled = TRUE, acq_rel = 1e-4))
  )
})

test_that("a run that becomes feasible late is not stopped at the transition", {
  # power depends on x1: infeasible for small x1, feasible for large x1. The
  # optimizer must be able to keep going after it first finds a feasible point,
  # not stop on the infeasible->feasible transition.
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = x[1], EN = sum((x - 0.5)^2) * 10)   # feasible only when x1 >= 0.8
    attr(res, "variance") <- c(power = 0.001, EN = 0.5); res
  }
  fit <- bo_calibrate(sim_fun = sim, bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
                      objective = "EN", constraints = list(power = c("ge", 0.8)),
                      n_init = 6, q = 1, budget = 20, seed = 1, progress = FALSE)
  # It should have run a meaningful number of evaluations and found a feasible design.
  expect_gt(nrow(fit$history), 8)
  expect_true(any(fit$history$feasible, na.rm = TRUE))
})
