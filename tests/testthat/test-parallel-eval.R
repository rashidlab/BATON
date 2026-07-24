# Task C2: simulator batch + initial design evaluate concurrently under
# options(BATON.cores) with results identical to the serial path.

rng_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
  # Deliberately NO set.seed here: draws come from the ambient RNG stream,
  # which invoke_evaluation seeds per evaluation via run_seeded(seed). This is
  # the hard case for serial/parallel equality.
  x <- theta$x
  v <- c(power = 0.85 + 0.1 * x + rnorm(1, sd = 0.01),
         EN = 50 - 10 * x + rnorm(1, sd = 0.1))
  attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
  v
}

run_small <- function() {
  bo_calibrate(
    sim_fun = rng_sim,
    bounds = list(x = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    n_init = 6,
    q = 2,
    budget = 12,
    seed = 99,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  )
}

test_that("simulator evaluations are identical serial vs parallel (RNG-using sim)", {
  skip_on_os("windows")
  skip_if_not_installed("parallel")
  old <- getOption("BATON.cores"); on.exit(options(BATON.cores = old), add = TRUE)
  options(BATON.cores = 1L); serial <- suppressMessages(run_small())
  options(BATON.cores = 2L); par2 <- suppressMessages(run_small())
  expect_equal(serial$history$objective, par2$history$objective)
  expect_equal(serial$history$metrics, par2$history$metrics)
  expect_equal(serial$history$feasible, par2$history$feasible)
  expect_equal(unlist(serial$best_theta), unlist(par2$best_theta))
})

test_that("history rows arrive in eval_id order after batched appends", {
  fit <- suppressMessages(run_small())
  expect_equal(fit$history$eval_id, seq_len(nrow(fit$history)))
})

test_that("a failing evaluation stops the run identically under parallel", {
  skip_on_os("windows")
  skip_if_not_installed("parallel")
  bad_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    if (theta$x > 0.9) stop("simulator exploded at x > 0.9")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  run_bad <- function() {
    bo_calibrate(
      sim_fun = bad_sim,
      bounds = list(x = c(0, 1)),
      objective = "EN",
      constraints = list(power = c("ge", 0.8)),
      n_init = 8, budget = 8, seed = 4,
      progress = FALSE, early_stop = list(enabled = FALSE)
    )
  }
  old <- getOption("BATON.cores"); on.exit(options(BATON.cores = old), add = TRUE)
  options(BATON.cores = 1L)
  expect_error(suppressMessages(run_bad()), "simulator exploded")
  options(BATON.cores = 2L)
  # suppressWarnings: mclapply warns "scheduled core encountered an error"
  # before evaluate_points re-throws; the error itself is what we assert.
  expect_error(suppressWarnings(suppressMessages(run_bad())), "simulator exploded")
})
