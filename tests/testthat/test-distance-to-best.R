# Task C7: the hybrid_staged proximity bonus needs the unit coordinates of
# the best feasible design. The old code matched history rows by float
# equality against the BPM incumbent (a posterior mean that equals no
# observed objective), so the lookup matched nothing and distance_to_best
# was stuck at its 1.0 default. best_feasible_unit finds the row by index,
# prefers high-fidelity rows when any feasible one exists (mirroring the A8
# incumbent semantics), and rescales theta under the CURRENT bounds so
# narrowed-bounds warm-starts do not mix coordinate systems.

test_that("best_feasible_unit returns the min-objective feasible row, rescaled", {
  bounds <- list(x = c(0, 10))
  history <- tibble::tibble(
    theta = list(list(x = 1), list(x = 5), list(x = 9)),
    unit_x = list(c(x = 0.1), c(x = 0.5), c(x = 0.9)),
    objective = c(5, 2, 3),
    feasible = c(TRUE, TRUE, FALSE),
    fidelity = c("low", "low", "low")
  )
  expect_equal(BATON:::best_feasible_unit(history, "objective", bounds),
               c(x = 0.5))
})

test_that("best_feasible_unit prefers high-fidelity feasible rows", {
  bounds <- list(x = c(0, 10))
  history <- tibble::tibble(
    theta = list(list(x = 2), list(x = 8)),
    unit_x = list(c(x = 0.2), c(x = 0.8)),
    objective = c(1, 4),  # low-fidelity row has the better objective
    feasible = c(TRUE, TRUE),
    fidelity = c("low", "high")
  )
  # the high-fidelity feasible row wins despite the worse objective
  expect_equal(BATON:::best_feasible_unit(history, "objective", bounds),
               c(x = 0.8))
})

test_that("best_feasible_unit falls back to any fidelity when no high row is feasible", {
  bounds <- list(x = c(0, 10))
  history <- tibble::tibble(
    theta = list(list(x = 2), list(x = 8)),
    unit_x = list(c(x = 0.2), c(x = 0.8)),
    objective = c(3, 1),
    feasible = c(TRUE, FALSE),
    fidelity = c("low", "high")
  )
  expect_equal(BATON:::best_feasible_unit(history, "objective", bounds),
               c(x = 0.2))
})

test_that("best_feasible_unit rescales under current bounds, not stored unit_x", {
  # A previous-stage history carries unit_x computed under WIDER old bounds
  # [0, 10]; the current run narrowed to [4, 6]. theta = 5 must map to 0.5 in
  # the CURRENT unit space, not the stale stored 0.5-of-[0,10].
  bounds <- list(x = c(4, 6))
  history <- tibble::tibble(
    theta = list(list(x = 5)),
    unit_x = list(c(x = 0.5)),  # stale: computed under [0, 10]... and
    # coincidentally 0.5 here too, so make it obviously wrong instead:
    objective = 2,
    feasible = TRUE,
    fidelity = "high"
  )
  history$unit_x <- list(c(x = 0.123))  # clearly not the current-bounds value
  expect_equal(BATON:::best_feasible_unit(history, "objective", bounds),
               c(x = 0.5))
})

test_that("best_feasible_unit returns bounds-ordered coordinates for reordered thetas", {
  # initial_history matches parameters by NAME, so a donor theta stored as
  # list(x2, x1) is valid. The returned unit vector must be in names(bounds)
  # order: candidate rows are bounds-ordered and the distance subtraction is
  # positional, so theta-ordered output would compare x1 against x2.
  bounds <- list(x1 = c(0, 10), x2 = c(0, 100))
  history <- tibble::tibble(
    theta = list(list(x2 = 50, x1 = 2)),  # reversed name order
    unit_x = list(c(x2 = 0.5, x1 = 0.2)),
    objective = 1,
    feasible = TRUE,
    fidelity = "high"
  )
  got <- BATON:::best_feasible_unit(history, "objective", bounds)
  expect_equal(names(got), c("x1", "x2"))
  expect_equal(got, c(x1 = 0.2, x2 = 0.5))
})

test_that("best_feasible_unit is NULL with no feasible rows", {
  bounds <- list(x = c(0, 1))
  history <- tibble::tibble(
    theta = list(list(x = 0.1), list(x = 0.9)),
    unit_x = list(c(x = 0.1), c(x = 0.9)),
    objective = c(5, 2),
    feasible = c(FALSE, FALSE),
    fidelity = c("low", "low")
  )
  expect_null(BATON:::best_feasible_unit(history, "objective", bounds))
})

test_that("best_feasible_unit ignores NA objectives", {
  bounds <- list(x = c(0, 1))
  history <- tibble::tibble(
    theta = list(list(x = 0.1), list(x = 0.5)),
    unit_x = list(c(x = 0.1), c(x = 0.5)),
    objective = c(NA_real_, 4),
    feasible = c(TRUE, TRUE),
    fidelity = c("low", "low")
  )
  expect_equal(BATON:::best_feasible_unit(history, "objective", bounds),
               c(x = 0.5))
})

test_that("hybrid_staged calibration runs with the index-based distance", {
  sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    if (!is.null(seed)) set.seed(seed)
    x <- theta$x
    v <- c(power = 0.85 + 0.1 * x + rnorm(1, sd = 0.005),
           EN = 50 - 10 * x + rnorm(1, sd = 0.05))
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    sim_fun = sim,
    bounds = list(x = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    n_init = 5, budget = 9, seed = 12,
    fidelity_method = "hybrid_staged",
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  ))
  expect_s3_class(fit, "BATON_fit")
  expect_equal(nrow(fit$history), 9)
})
