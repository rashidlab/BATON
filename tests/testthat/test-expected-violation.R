# Regression tests for compute_expected_violation (review item 1.3 / Task A4).
# Before the fix the expression collapsed to sd*dnorm(z): symmetric about the
# boundary and ~0 for deeply infeasible points, so it could not guide the search
# toward the feasibility boundary.

# Closed form: for a "ge" constraint (metric >= threshold), the expected
# violation E[max(0, threshold - Y)] with Y ~ N(mu, sd^2) is
#   sd * (z * Phi(z) + phi(z)),  z = (threshold - mu) / sd.
ev_closed <- function(mu, sd, threshold, direction) {
  z <- if (direction == "ge") (threshold - mu) / sd else (mu - threshold) / sd
  sd * (z * stats::pnorm(z) + stats::dnorm(z))
}

test_that("expected violation matches the closed form (ge and le)", {
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.1)))
  pred <- list(power = list(mean = c(0.6, 0.85), sd = c(0.05, 0.05)),
               type1 = list(mean = c(0.2, 0.05), sd = c(0.02, 0.02)))
  got <- BATON:::compute_expected_violation(pred, ctbl, c("power", "type1"))
  expected <- ev_closed(0.6, 0.05, 0.8, "ge")  + ev_closed(0.2, 0.02, 0.1, "le")
  expected2 <- ev_closed(0.85, 0.05, 0.8, "ge") + ev_closed(0.05, 0.02, 0.1, "le")
  expect_equal(got[1], expected,  tolerance = 1e-6)
  expect_equal(got[2], expected2, tolerance = 1e-6)
})

test_that("expected violation grows as a point becomes more infeasible", {
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8)))
  # Three power means: deeply infeasible (0.2) < boundary (0.8) < feasible (0.95).
  pred <- list(power = list(mean = c(0.2, 0.8, 0.95), sd = c(0.05, 0.05, 0.05)))
  ev <- BATON:::compute_expected_violation(pred, ctbl, "power")
  # Monotone: deep infeasible has the LARGEST expected violation; feasible the least.
  expect_true(ev[1] > ev[2])           # deep > boundary (bug had these ~equal/reversed)
  expect_true(ev[2] > ev[3])           # boundary > feasible
  # A deeply infeasible point's violation approaches the true gap (0.8 - 0.2 = 0.6).
  expect_gt(ev[1], 0.5)                # bug gave ~0 here
})
