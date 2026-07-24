# Tests for BATON_fit$status recording why a run ended (Task D2).

# Deterministic, always-feasible simulator (same shape as test-early-stop.R).
status_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
  x <- as.numeric(unlist(theta))
  res <- c(power = 0.9, EN = sum((x - 0.5)^2) * 10)
  attr(res, "variance") <- c(power = 0.001, EN = 0.5)
  res
}

test_that("status is budget_exhausted when the budget is consumed", {
  skip_if_not_installed("DiceKriging")
  fit <- suppressMessages(bo_calibrate(
    sim_fun = status_sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 5, q = 1, budget = 8, seed = 1, progress = FALSE,
    early_stop = list(enabled = FALSE)
  ))
  expect_equal(fit$status, "budget_exhausted")
})

test_that("status is early_stopped when the improvement plateau trips", {
  skip_if_not_installed("DiceKriging")
  # threshold = 10 means any relative improvement < 1000% counts as "no
  # improvement", so the plateau check trips at the first eligible iteration
  # (iteration 4, after min_iters_before_stopping = 3).
  fit <- suppressMessages(bo_calibrate(
    sim_fun = status_sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 5, q = 1, budget = 20, seed = 1, progress = FALSE,
    early_stop = list(enabled = TRUE, patience = 1, threshold = 10,
                      consecutive = 1)
  ))
  expect_equal(fit$status, "early_stopped")
  # Confirms it actually stopped early rather than exhausting the budget.
  expect_lt(nrow(fit$history), 20)
})

test_that("status is acq_flatline when the acquisition threshold stop trips", {
  skip_if_not_installed("DiceKriging")
  # acq_rel = 1e6 makes the flatline threshold enormous, so the acquisition
  # check trips after its 3-iteration patience. patience = 50 keeps the
  # improvement-plateau check from ever firing first.
  fit <- suppressMessages(bo_calibrate(
    sim_fun = status_sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 5, q = 1, budget = 20, seed = 1, progress = FALSE,
    early_stop = list(enabled = TRUE, patience = 50, acq_rel = 1e6)
  ))
  expect_equal(fit$status, "acq_flatline")
  expect_lt(nrow(fit$history), 20)
})

test_that("build_baton_fit rejects an undocumented status", {
  expect_error(
    build_baton_fit(
      history = initialise_history(),
      best_theta = NULL, best_objective = NA_real_, surrogates = NULL,
      policies = NULL, diagnostics = NULL, multi_seed_summary = NULL,
      multi_seed_runs = NULL, verdict = NA_character_, bounds = NULL,
      constraints = NULL, constraint_tbl = NULL,
      status = "not_a_real_status"
    ),
    "Invalid status"
  )
})

test_that("build_baton_fit rejects partial matches of documented statuses", {
  expect_error(
    build_baton_fit(
      history = initialise_history(),
      best_theta = NULL, best_objective = NA_real_, surrogates = NULL,
      policies = NULL, diagnostics = NULL, multi_seed_summary = NULL,
      multi_seed_runs = NULL, verdict = NA_character_, bounds = NULL,
      constraints = NULL, constraint_tbl = NULL,
      status = "budget"
    ),
    "Invalid status"
  )
})

test_that("build_baton_fit accepts every documented status value", {
  for (s in c("budget_exhausted", "early_stopped", "acq_flatline",
              "cancelled", "walltime", "errored", "checkpoint", "completed")) {
    fit <- build_baton_fit(
      history = initialise_history(),
      best_theta = NULL, best_objective = NA_real_, surrogates = NULL,
      policies = NULL, diagnostics = NULL, multi_seed_summary = NULL,
      multi_seed_runs = NULL, verdict = NA_character_, bounds = NULL,
      constraints = NULL, constraint_tbl = NULL,
      status = s
    )
    expect_s3_class(fit, "BATON_fit")
    expect_equal(fit$status, s)
  }
})

test_that("an errored philosophies fit is built by the helper with back-compat error field", {
  skip_if_not_installed("DiceKriging")
  philosophies <- list(
    BadPhil = list(objective = "EN",
                   constraints = list(missing_metric = c("ge", 0.8)))
  )
  res <- suppressMessages(expect_warning(
    bo_calibrate_philosophies(
      scenario_id = "status-test",
      sim_fun = status_sim,
      bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
      philosophies = philosophies,
      dependency_order = c(BadPhil = "donor"),
      multi_seed_verify = FALSE,
      n_init = 5, q = 1, budget = 8, progress = FALSE
    ),
    "FAILED"
  ))
  fit <- res$fits$BadPhil
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$status, "errored")
  expect_equal(fit$verdict, "ERRORED")
  # Helper-built object: full field set and c("BATON_fit", "list") class.
  expect_true(all(c("history", "surrogates", "policies", "error_message")
                  %in% names(fit)))
  expect_identical(class(fit), c("BATON_fit", "list"))
  expect_true(is.character(fit$error_message) && nzchar(fit$error_message))
  # Back-compat: the old `error` field still carries the message.
  expect_identical(fit$error, fit$error_message)
})
