# Tests for RNG hygiene (Task D3): bo_calibrate() must restore the caller's
# RNG stream on exit, on both success and error paths, for both seed modes.

# Deterministic simulator with variance attr (same shape as test-status.R).
rng_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
  x <- as.numeric(unlist(theta))
  res <- c(power = 0.9, EN = sum((x - 0.5)^2) * 10)
  attr(res, "variance") <- c(power = 0.001, EN = 0.5)
  res
}

tiny_run <- function(seed) {
  invisible(suppressMessages(bo_calibrate(
    sim_fun = rng_sim,
    bounds = list(x1 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 3, q = 1, budget = 4, seed = seed, progress = FALSE,
    early_stop = list(enabled = FALSE)
  )))
}

test_that("bo_calibrate leaves the caller's RNG stream untouched (explicit seed)", {
  skip_if_not_installed("DiceKriging")
  set.seed(123)
  expected <- runif(3)
  set.seed(123)
  tiny_run(seed = 42)
  expect_equal(runif(3), expected)
})

test_that("seed = NULL advances the caller's stream by exactly one draw", {
  skip_if_not_installed("DiceKriging")
  # Expected state: the single sample.int(1e6, 1) seed draw, nothing else.
  set.seed(123)
  invisible(sample.int(1e6, 1))
  expected_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expected_next <- runif(3)

  set.seed(123)
  tiny_run(seed = NULL)
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                   expected_state)
  expect_equal(runif(3), expected_next)
})

test_that("seed = NULL runs are independent replicates", {
  skip_if_not_installed("DiceKriging")
  # If the restore hook rewound past the seed draw, both runs would draw the
  # same internal seed and be identical. The one-draw advancement guarantees
  # different seeds.
  set.seed(123)
  run_null <- function() {
    suppressMessages(bo_calibrate(
      sim_fun = rng_sim,
      bounds = list(x1 = c(0, 1)),
      objective = "EN", constraints = list(power = c("ge", 0.8)),
      n_init = 3, q = 1, budget = 4, seed = NULL, progress = FALSE,
      early_stop = list(enabled = FALSE)
    ))
  }
  fit1 <- run_null()
  fit2 <- run_null()
  expect_false(identical(fit1$policies$seed, fit2$policies$seed))
})

test_that("caller's stream is restored even when bo_calibrate errors", {
  set.seed(123)
  expected <- runif(3)
  set.seed(123)
  expect_error(suppressMessages(bo_calibrate(
    sim_fun = function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
      stop("simulator exploded")
    },
    bounds = list(x1 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 3, q = 1, budget = 4, seed = 7, progress = FALSE,
    early_stop = list(enabled = FALSE)
  )), regexp = "simulator exploded")
  expect_equal(runif(3), expected)
})
