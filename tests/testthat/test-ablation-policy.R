# Ablation policies may name only a subset of fidelity levels; bo_calibrate's
# contract requires exactly low/med/high. normalize_fidelity_policy maps a
# partial policy onto the full contract: each missing level inherits the
# nearest provided level at or below it (else the smallest provided).

test_that("single-level policies replicate their value across all levels", {
  expect_equal(BATON:::normalize_fidelity_policy(c(low = 500), "low_only"),
               c(low = 500, med = 500, high = 500))
  expect_equal(BATON:::normalize_fidelity_policy(c(med = 2000), "med_only"),
               c(low = 2000, med = 2000, high = 2000))
  expect_equal(BATON:::normalize_fidelity_policy(c(high = 10000), "high_only"),
               c(low = 10000, med = 10000, high = 10000))
})

test_that("partial policies inherit from the level below", {
  expect_equal(BATON:::normalize_fidelity_policy(c(low = 200, high = 1000), "lh"),
               c(low = 200, med = 200, high = 1000))
  expect_equal(BATON:::normalize_fidelity_policy(c(med = 2000, high = 8000), "mh"),
               c(low = 2000, med = 2000, high = 8000))
})

test_that("full policies pass through unchanged", {
  expect_equal(
    BATON:::normalize_fidelity_policy(c(low = 200, med = 1000, high = 10000), "full"),
    c(low = 200, med = 1000, high = 10000)
  )
})

test_that("invalid policy names error clearly", {
  expect_error(BATON:::normalize_fidelity_policy(c(LOW = 200), "bad"),
               "subset of 'low', 'med', 'high'")
  expect_error(BATON:::normalize_fidelity_policy(c(200), "unnamed"),
               "subset of 'low', 'med', 'high'")
  expect_error(BATON:::normalize_fidelity_policy(c(low = 200, low = 300), "dup"),
               "duplicate")
})

test_that("ablation_multifidelity runs single- and two-level policies end to end", {
  sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    if (!is.null(seed)) set.seed(seed)
    x <- theta$x
    v <- c(power = 0.9 + 0.05 * x, EN = 50 - 10 * x + rnorm(1, sd = 0.05))
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  abl <- suppressMessages(ablation_multifidelity(
    sim_fun = sim,
    bounds = list(x = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    policies = list(low_only = c(low = 200), lh = c(low = 200, high = 1000)),
    seeds = 1,
    bo_args = list(n_init = 3, q = 1, budget = 4, progress = FALSE,
                   early_stop = list(enabled = FALSE)),
    progress = FALSE
  ))
  expect_s3_class(abl, "BATON_multifidelity")
  expect_equal(sort(unique(abl$runs$policy)), c("lh", "low_only"))
})
