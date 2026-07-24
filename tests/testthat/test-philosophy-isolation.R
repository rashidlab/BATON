# Test for per-philosophy error isolation (review item 1.10 / Task A10). One
# philosophy's failure must not discard the other completed (multi-hour) fits.

test_that("a failing philosophy does not abort the whole batch", {
  skip_if_not_installed("DiceKriging")
  good_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = 0.85, EN = sum((x - 0.5)^2) * 10)
    attr(res, "variance") <- c(power = 0.001, EN = 0.5)
    res
  }
  philosophies <- list(
    # BadPhil references a constraint metric the simulator never returns, so its
    # bo_calibrate() errors at the first evaluation.
    BadPhil  = list(objective = "EN", constraints = list(missing_metric = c("ge", 0.8))),
    GoodPhil = list(objective = "EN", constraints = list(power = c("ge", 0.8)))
  )
  res <- expect_warning(
    bo_calibrate_philosophies(
      scenario_id = "test",
      sim_fun = good_sim,
      bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
      philosophies = philosophies,
      dependency_order = c(BadPhil = "donor", GoodPhil = "donor"),
      multi_seed_verify = FALSE,
      n_init = 5, q = 1, budget = 8, progress = FALSE
    ),
    "FAILED|failed"
  )
  # The batch completes; the good philosophy has a real fit, the bad one is marked errored.
  expect_s3_class(res$fits$GoodPhil, "BATON_fit")
  expect_true(is.finite(res$fits$GoodPhil$best_objective))
  expect_equal(res$fits$BadPhil$status, "errored")
  expect_equal(res$manifest$verdict[res$manifest$philosophy == "BadPhil"], "ERRORED")
})

test_that("recipient warmstart from a donor completes (build_donor_spec theta_id path)", {
  skip_if_not_installed("DiceKriging")
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = 0.82 + 0.05 * x[1], EN = sum((x - 0.5)^2) * 10)
    attr(res, "variance") <- c(power = 0.001, EN = 0.5)
    res
  }
  philosophies <- list(
    Fleming = list(objective = "EN", constraints = list(power = c("ge", 0.8))),
    Other   = list(objective = "EN", constraints = list(power = c("ge", 0.8)))
  )
  # Fleming is a default donor; Other becomes a recipient warmstarted from it.
  # multi_seed_strict = FALSE so the NA-verdict donor is NOT filtered out and
  # build_donor_spec (the theta_id lookup) is actually exercised.
  res <- bo_calibrate_philosophies(
    scenario_id = "ws", sim_fun = sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    philosophies = philosophies,
    multi_seed_verify = FALSE, multi_seed_strict = FALSE,
    n_init = 5, q = 1, budget = 8, progress = FALSE
  )
  expect_s3_class(res$fits$Fleming, "BATON_fit")
  expect_s3_class(res$fits$Other, "BATON_fit")
  expect_true(is.finite(res$fits$Other$best_objective))
  expect_equal(res$manifest$dependency_role[res$manifest$philosophy == "Other"],
               "recipient")
  # Confirm build_donor_spec actually ran: the recipient records Fleming as a
  # warmstart donor (empty string would mean it fell back to LHS-only).
  expect_true(nzchar(res$manifest$warmstart_donors[res$manifest$philosophy == "Other"]))
})
