# Regression tests for the warm-start fix (review item 1.1 / Task A1).
# Before the fix, parinit was passed inside control = list(...), which
# DiceKriging::km() silently ignores, so warm-start was a no-op.

test_that("km honors parinit as a top-level argument (not swallowed by control)", {
  skip_if_not_installed("DiceKriging")
  set.seed(1)
  X <- matrix(runif(30), ncol = 3)
  y <- rowSums(X) + rnorm(10, sd = 0.01)
  m1 <- DiceKriging::km(design = as.data.frame(X), response = y,
                        covtype = "matern5_2", control = list(trace = FALSE))
  hp <- as.numeric(m1@covariance@range.val)

  # Passing parinit at top level must seed the optimizer start exactly.
  m2 <- do.call(DiceKriging::km, list(
    design = X, response = y, covtype = "matern5_2",
    nugget = 1e-6, nugget.estim = FALSE,
    control = list(trace = FALSE), parinit = hp))
  expect_equal(as.numeric(m2@parinit), hp, tolerance = 1e-8)

  # Control-nested parinit (the old bug) is ignored: parinit != hp.
  m3 <- DiceKriging::km(design = X, response = y, covtype = "matern5_2",
                        nugget = 1e-6, nugget.estim = FALSE,
                        control = list(trace = FALSE, parinit = hp))
  expect_false(isTRUE(all.equal(as.numeric(m3@parinit), hp)))
})

test_that("fit_surrogates forwards prev_surrogates hyperparameters as km parinit", {
  skip_if_not_installed("DiceKriging")
  # A real (small) history via a tiny calibration; no variance attached so the
  # km (not hetGP) path is exercised. use_hetgp = FALSE forces DiceKriging.
  sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL, ...) {
    x <- vapply(theta, as.numeric, numeric(1))
    c(power = 0.9 - 0.6 * (x[["x1"]] - 0.5)^2,
      type1 = 0.08 + 0.3 * (x[["x2"]] - 0.25)^2,
      EN = (x[["x1"]] - 0.5)^2 + (x[["x2"]] - 0.25)^2)
  }
  set.seed(7)
  fit <- bo_calibrate(sim_fun = sim, bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
                      objective = "EN",
                      constraints = list(power = c("ge", 0.8), type1 = c("le", 0.12)),
                      n_init = 8, q = 1, budget = 10, seed = 7, progress = FALSE)
  ctbl <- BATON:::parse_constraints(list(power = c("ge", 0.8), type1 = c("le", 0.12)))

  s1 <- fit_surrogates(fit$history, "EN", ctbl, use_hetgp = FALSE)
  s2 <- fit_surrogates(fit$history, "EN", ctbl, use_hetgp = FALSE,
                       prev_surrogates = s1)
  # The warm fit must seed its optimizer at the cold fit's lengthscales.
  expect_s4_class(s2$EN, "km")
  expect_equal(as.numeric(s2$EN@parinit),
               as.numeric(s1$EN@covariance@range.val), tolerance = 1e-8)
})

# BATON:::run_seeded() is the primitive that makes every wrapped fit (the main km/hetGP
# calls AND the hard-to-trigger fallback fits) reproducible regardless of
# RNG-stream position or worker forking. Testing it directly is far more reliable
# than trying to trigger each fallback path, some of which fire only on degenerate
# data where the fit outcome is insensitive to the random start anyway. If this
# primitive holds, core-count independence holds for all fits that use it (the
# end-to-end guarantee is checked by the BATON.cores test in test-regression-golden.R).
test_that("run_seeded makes a random computation reproducible", {
  a <- BATON:::run_seeded(42, function() stats::runif(5))
  b <- BATON:::run_seeded(42, function() stats::runif(5))
  expect_equal(a, b)                                   # same seed -> same draws
  c <- BATON:::run_seeded(43, function() stats::runif(5))
  expect_false(isTRUE(all.equal(a, c)))                # different seed -> different
})

test_that("run_seeded restores the caller's RNG stream (no side effects)", {
  set.seed(100)
  before <- .Random.seed
  invisible(BATON:::run_seeded(7, function() stats::runif(10)))
  expect_identical(.Random.seed, before)               # stream position unchanged
  # And the caller's next draw is exactly what it would have been with no call.
  set.seed(100); expected <- stats::runif(3)
  set.seed(100); invisible(BATON:::run_seeded(7, function() stats::runif(10)))
  actual <- stats::runif(3)
  expect_equal(actual, expected)
})

test_that("run_seeded with NULL seed does not touch the RNG", {
  set.seed(100); expected <- stats::runif(3)
  set.seed(100)
  first <- BATON:::run_seeded(NULL, function() 123)            # no seeding
  expect_equal(first, 123)
  expect_equal(stats::runif(3), expected)              # stream untouched by the call
})

# Task C3: warm-start must actually reduce (or at minimum not increase) the
# cost of a surrogate fit. This is the payoff check for A1/A2: parinit seeding
# starts the MLE at the previous optimum, so the optimizer should converge in
# fewer iterations on a mid-size design.

c3_history <- function() {
  set.seed(31)
  n <- 120; d <- 6
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- paste0("x", seq_len(d))
  y <- rowSums(sin(3 * X)) + rnorm(n, sd = 0.05)
  tibble::tibble(
    unit_x = lapply(seq_len(n), function(i) X[i, ]),
    theta_id = apply(X, 1, BATON:::theta_to_id),
    metrics = lapply(y, function(v) c(EN = v)),
    variance = replicate(n, c(EN = NA_real_), simplify = FALSE)
  )
}

# Deterministic half of the payoff check: the warm fit provably starts at the
# cold fit's lengthscales, so any timing comparison is warm-vs-cold, not
# cold-vs-cold. No timing dependence; always runs.
test_that("mid-size warm fit is seeded at the cold fit's lengthscales", {
  skip_if_not_installed("DiceKriging")
  history <- c3_history()
  ctbl <- BATON:::parse_constraints(list())
  cold <- fit_surrogates(history, "EN", ctbl, use_hetgp = FALSE, fit_seed = 1)
  warm <- fit_surrogates(history, "EN", ctbl, use_hetgp = FALSE, fit_seed = 1,
                         prev_surrogates = cold)
  expect_equal(as.numeric(warm$EN@parinit),
               as.numeric(cold$EN@covariance@range.val), tolerance = 1e-8)
})

# Timing half: wall-clock assertions are only meaningful on an unloaded
# machine, so this runs on local dev boxes only (skip_on_ci; CI runners under
# shared load can blow any margin regardless of warm-start behavior). Medians
# over 5 reps and a 1.5x margin catch "warm-start makes fitting slower"
# regressions, not exact speedup factors (measured ~1-9% faster depending on
# noise structure).
test_that("warm-started km fit is not slower than a cold fit (mid-size design)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("DiceKriging")
  history <- c3_history()
  ctbl <- BATON:::parse_constraints(list())

  time_median <- function(expr_fn, reps = 5) {
    stats::median(vapply(seq_len(reps), function(i) {
      unname(system.time(expr_fn())["elapsed"])
    }, numeric(1)))
  }

  cold <- fit_surrogates(history, "EN", ctbl, use_hetgp = FALSE, fit_seed = 1)
  t_cold <- time_median(function()
    fit_surrogates(history, "EN", ctbl, use_hetgp = FALSE, fit_seed = 1))
  t_warm <- time_median(function()
    fit_surrogates(history, "EN", ctbl, use_hetgp = FALSE, fit_seed = 1,
                   prev_surrogates = cold))

  expect_lte(t_warm, t_cold * 1.5)
})
