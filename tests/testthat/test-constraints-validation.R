# Tests for parse_constraints input validation (review item A12).

test_that("parse_constraints rejects unnamed constraint lists", {
  expect_error(BATON:::parse_constraints(list(c("ge", 0.8))), "named list")
})

test_that("parse_constraints rejects non-numeric / non-finite thresholds", {
  expect_error(BATON:::parse_constraints(list(power = c("ge", "high"))), "finite number")
  expect_error(BATON:::parse_constraints(list(power = c("ge", Inf))), "finite number")
})

test_that("parse_constraints rejects malformed (too-short) constraints", {
  expect_error(BATON:::parse_constraints(list(power = "ge")), "direction, threshold")
})

test_that("parse_constraints accepts valid input", {
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.1)))
  expect_equal(ctbl$metric, c("power", "type1"))
  expect_equal(unname(ctbl$threshold), c(0.8, 0.1))
  expect_equal(unname(ctbl$direction), c("ge", "le"))
})
