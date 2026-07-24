# Test for is_feasible robustness (review item D6 / A12). On an atomic named
# vector, metrics[[missing_name]] throws "subscript out of bounds" rather than
# returning NULL, which surfaced as an opaque crash deep in the BO loop instead
# of a clear name-contract error.

test_that("is_feasible returns FALSE (not error) when a metric name is missing", {
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.1)))
  # Atomic vector missing 'power' (e.g. a mis-named simulator output).
  metrics <- c(type1 = 0.05, EN = 40)
  expect_false(BATON:::is_feasible(metrics, ctbl))       # was: subscript error
})

test_that("is_feasible works for both list and atomic metrics", {
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.1)))
  ok_vec  <- c(power = 0.9, type1 = 0.05)
  ok_list <- list(power = 0.9, type1 = 0.05)
  expect_true(BATON:::is_feasible(ok_vec, ctbl))
  expect_true(BATON:::is_feasible(ok_list, ctbl))
  expect_false(BATON:::is_feasible(c(power = 0.7, type1 = 0.05), ctbl))  # power fails
})
