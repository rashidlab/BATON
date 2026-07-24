# Test for the incumbent high-fidelity fallback (review item 1.7 / Task A8).
# When feasible points exist only at low/med fidelity, best_feasible_objective
# returned Inf under high_fidelity_only, flipping acquisition into "no feasible
# solution yet" mode despite known feasible designs.

test_that("incumbent falls back to any-fidelity feasible best when no high-fid feasible", {
  history <- data.frame(
    feasible  = c(TRUE, TRUE, FALSE),
    fidelity  = c("low", "med", "high"),
    objective = c(0.50, 0.40, 0.20),
    stringsAsFactors = FALSE
  )
  # Feasible points are low/med only; the single high-fid eval is infeasible.
  val <- BATON:::best_feasible_objective(
    history, "objective", surrogates = NULL,
    high_fidelity_only = TRUE, incumbent_method = "boi")
  expect_equal(val, 0.40)          # best feasible at any fidelity, NOT Inf
  expect_true(is.finite(val))
})

test_that("incumbent uses high-fid feasible best when one exists", {
  history <- data.frame(
    feasible  = c(TRUE, TRUE, TRUE),
    fidelity  = c("low", "high", "high"),
    objective = c(0.30, 0.45, 0.50),
    stringsAsFactors = FALSE
  )
  val <- BATON:::best_feasible_objective(
    history, "objective", surrogates = NULL,
    high_fidelity_only = TRUE, incumbent_method = "boi")
  expect_equal(val, 0.45)          # min over high-fid feasible (0.45, 0.50)
})

test_that("incumbent is Inf only when no feasible point exists at any fidelity", {
  history <- data.frame(
    feasible  = c(FALSE, FALSE),
    fidelity  = c("low", "high"),
    objective = c(0.30, 0.20),
    stringsAsFactors = FALSE
  )
  val <- BATON:::best_feasible_objective(
    history, "objective", surrogates = NULL,
    high_fidelity_only = TRUE, incumbent_method = "boi")
  expect_equal(val, Inf)
})
