# Tests for the Welford / Chan variance utilities (review item 1.11 / Task A5).
# These functions feed the heteroskedastic GP noise model, so silent corruption
# here propagates everywhere. The review found zero prior coverage.

test_that("welford_mean_var matches metrics by name, not position", {
  # sample_fn returns the two metrics in ALTERNATING order each call.
  sf <- function(i, ...) if (i %% 2 == 1) c(a = 1, b = 10) else c(b = 10, a = 1)
  res <- BATON:::welford_mean_var(sf, 100)
  # Correct name-matching -> a is always 1, b always 10 (each constant).
  expect_equal(res$mean[["a"]], 1)
  expect_equal(res$mean[["b"]], 10)
  expect_lt(res$variance[["a"]], 1e-12)   # constant metric -> ~0 variance
  expect_lt(res$variance[["b"]], 1e-12)
})

test_that("welford_mean_var errors on genuinely mismatched names", {
  sf <- function(i, ...) if (i %% 2 == 1) c(a = 1, b = 2) else c(a = 1, c = 2)
  expect_error(BATON:::welford_mean_var(sf, 10), "name")
})

test_that("pool_welford_results handles n=1 chunks without producing NA", {
  set.seed(1)
  big <- BATON:::welford_mean_var(function(i, ...) c(m = stats::rnorm(1, 5, 1)), 100)
  one <- list(mean = c(m = 4.0), variance = c(m = NA_real_), n = 1)
  pooled <- BATON:::pool_welford_results(list(big, one))
  expect_false(is.na(pooled$variance[["m"]]))
  expect_equal(pooled$n, 101)
})

test_that("pooled mean and variance match direct computation (n>1 chunks)", {
  set.seed(3)
  x <- stats::rnorm(200)
  idx <- split(seq_along(x), rep(1:4, each = 50))
  chunks <- lapply(idx, function(ii) {
    vals <- x[ii]; k <- 0
    BATON:::welford_mean_var(function(i, ...) { k <<- k + 1; c(m = vals[k]) },
                             length(ii))
  })
  pooled <- BATON:::pool_welford_results(chunks)
  expect_equal(pooled$mean[["m"]], mean(x), tolerance = 1e-10)
  # welford returns variance OF THE MEAN (M2 / (n(n-1)) = sample_var / n).
  expect_equal(pooled$variance[["m"]], stats::var(x) / length(x), tolerance = 1e-8)
})
