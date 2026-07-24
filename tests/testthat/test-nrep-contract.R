# Task B1: invoke_simulator passes n_rep to simulators that accept it
# (optional-argument contract extension; PI decision 2026-07-22).

test_that("invoke_simulator forwards n_rep to a simulator with an n_rep argument", {
  received <- NULL
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL) {
    received <<- n_rep
    c(power = 0.9, EN = 40)
  }
  BATON:::invoke_simulator(sim_fun = sim, theta = list(x = 0.5),
                           fidelity = "high", n_rep = 12345L, seed = 1)
  expect_equal(received, 12345L)
})

test_that("invoke_simulator does NOT pass n_rep to a dots-only simulator", {
  # Dots-only acceptance was rejected: forwarding wrappers (e.g. the one
  # fix_parameters() returns) pass their dots verbatim to an inner simulator,
  # so injecting n_rep into dots leaks it into legacy simulators that cannot
  # accept it ("unused argument" error). Only an explicit n_rep formal opts in.
  received <- "unset"
  sim <- function(theta, fidelity = "high", seed = NULL, ...) {
    dots <- list(...)
    received <<- dots$n_rep
    c(power = 0.9, EN = 40)
  }
  BATON:::invoke_simulator(sim_fun = sim, theta = list(x = 0.5),
                           fidelity = "low", n_rep = 200L, seed = 1)
  expect_null(received)
})

test_that("fix_parameters wrapper over a legacy simulator still runs", {
  legacy_sim <- function(theta, fidelity = "high", seed = NULL) {
    c(power = 0.9, EN = 40 + theta$x)
  }
  wrapper <- suppressMessages(fix_parameters(
    sim_fun = legacy_sim,
    fixed_params = list(y = 0.25),
    original_bounds = list(x = c(0, 1), y = c(0, 1))
  ))
  res <- BATON:::invoke_simulator(sim_fun = wrapper, theta = list(x = 0.5),
                                  fidelity = "high", n_rep = 500L, seed = 1)
  expect_equal(unname(res$metrics[["EN"]]), 40.5)
})

test_that("fix_parameters wrapper forwards n_rep to an n_rep-aware simulator", {
  received <- NULL
  modern_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL) {
    received <<- n_rep
    c(power = 0.9, EN = 40)
  }
  wrapper <- suppressMessages(fix_parameters(
    sim_fun = modern_sim,
    fixed_params = list(y = 0.25),
    original_bounds = list(x = c(0, 1), y = c(0, 1))
  ))
  BATON:::invoke_simulator(sim_fun = wrapper, theta = list(x = 0.5),
                           fidelity = "high", n_rep = 750L, seed = 1)
  expect_equal(received, 750L)
})

test_that("legacy simulator without n_rep or dots still runs unchanged", {
  sim <- function(theta, fidelity = "high", seed = NULL) {
    c(power = 0.9, EN = 40)
  }
  res <- BATON:::invoke_simulator(sim_fun = sim, theta = list(x = 0.5),
                                  fidelity = "med", n_rep = 1000L, seed = 1)
  expect_equal(unname(res$metrics[["power"]]), 0.9)
  expect_equal(res$n_rep, 1000L)
})

test_that("fix_parameters wrapper preserves the inner simulator's n_rep default", {
  received <- NULL
  modern_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = 5000) {
    received <<- n_rep
    c(power = 0.9, EN = 40)
  }
  wrapper <- suppressMessages(fix_parameters(
    sim_fun = modern_sim,
    fixed_params = list(y = 0.25),
    original_bounds = list(x = c(0, 1), y = c(0, 1))
  ))
  # Called without n_rep (e.g. direct use or multi-seed verification): the
  # inner simulator's own default must apply, not a forwarded NULL.
  wrapper(theta = list(x = 0.5), fidelity = "high", seed = 1)
  expect_equal(received, 5000)
})

test_that("simulator sees the fidelity-level n_rep through a bo_calibrate run", {
  seen <- integer()
  sim <- function(theta, fidelity = "low", seed = NULL, n_rep = NULL, ...) {
    seen <<- c(seen, as.integer(n_rep))
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
    n_init = 4,
    budget = 6,
    seed = 42,
    fidelity_levels = c(low = 111, med = 222, high = 333),
    fidelity_method = "staged",
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  ))
  expect_true(length(seen) >= 5)
  expect_false(anyNA(seen))
  # every received n_rep must be one of the declared fidelity levels
  expect_true(all(seen %in% c(111L, 222L, 333L)))
  # and the recorded history n_rep must match what the simulator received
  expect_equal(fit$history$n_rep, seen[seq_len(nrow(fit$history))])
})
