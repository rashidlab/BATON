# Task C6: matrix-based history aggregation. The grouped means must preserve
# mean(..., na.rm = TRUE) semantics exactly, including the NaN result for
# all-NA groups (plain rowsum would propagate NA and poison surrogate inputs).

test_that("na_aware_group_mean matches per-group mean(na.rm = TRUE)", {
  set.seed(5)
  ids <- sample(c("a", "b", "c", "d"), 40, replace = TRUE)
  x <- rnorm(40)
  x[sample(40, 8)] <- NA
  ref <- tapply(x, ids, function(v) mean(v, na.rm = TRUE))
  got <- BATON:::na_aware_group_mean(x, ids)
  expect_equal(got[names(ref)], c(ref)[names(ref)])
})

test_that("all-NA groups yield NaN (matching mean(na.rm = TRUE))", {
  ids <- c("g1", "g1", "g2", "g2", "g3")
  x <- c(1, 3, NA, NA, 7)
  got <- BATON:::na_aware_group_mean(x, ids)
  expect_equal(got[["g1"]], 2)
  expect_true(is.nan(got[["g2"]]))
  expect_equal(got[["g3"]], 7)
})

test_that("groups come back ordered by sorted unique id (split()-compatible)", {
  ids <- c("zz", "aa", "zz", "mm")
  x <- c(1, 2, 3, 4)
  got <- BATON:::na_aware_group_mean(x, ids)
  expect_equal(names(got), sort(unique(ids)))
  expect_equal(names(got), names(split(seq_along(ids), ids)))
})

test_that("history without a variance column still gets a real km fit", {
  skip_if_not_installed("DiceKriging")
  # variance is documented optional: its absence must take the homoskedastic
  # nugget path, not error inside aggregation and silently fall back to a
  # constant predictor.
  set.seed(17)
  n <- 15
  X <- matrix(runif(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  y <- (X[, 1] - 0.5)^2 + rnorm(n, sd = 0.01)
  history <- tibble::tibble(
    unit_x = lapply(seq_len(n), function(i) X[i, ]),
    theta_id = apply(X, 1, BATON:::theta_to_id),
    metrics = lapply(y, function(v) c(EN = v))
    # no variance column at all
  )
  s <- fit_surrogates(history, "EN", BATON:::parse_constraints(list()),
                      use_hetgp = FALSE, fit_seed = 2)
  expect_s4_class(s$EN, "km")
})

test_that("fit_surrogates aggregates replicated thetas with NA metric and noise entries", {
  skip_if_not_installed("DiceKriging")
  set.seed(13)
  n_unique <- 12
  X <- matrix(runif(n_unique * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  # duplicate every design point, then inject NAs: one metric NA in a
  # replicated pair (mean must use the remaining value), one noise NA
  rows <- rep(seq_len(n_unique), each = 2)
  y <- (X[rows, 1] - 0.5)^2 + rnorm(length(rows), sd = 0.01)
  y[3] <- NA  # theta 2, first replicate: aggregate must equal y[4]
  noise <- rep(1e-4, length(rows))
  noise[7] <- NA  # theta 4, first replicate
  history <- tibble::tibble(
    unit_x = lapply(rows, function(i) X[i, ]),
    theta_id = apply(X, 1, BATON:::theta_to_id)[rows],
    metrics = lapply(y, function(v) c(EN = v)),
    variance = lapply(noise, function(v) c(EN = v))
  )
  s <- fit_surrogates(history, "EN", BATON:::parse_constraints(list()),
                      use_hetgp = FALSE, fit_seed = 3)
  expect_s4_class(s$EN, "km")
  # 12 unique designs, all with at least one non-NA observation -> 12 rows
  expect_equal(nrow(s$EN@X), n_unique)
  # the model interpolates near the aggregate of the half-NA pair, i.e. y[4]
  p <- BATON:::predict_surrogates(s, matrix(X[2, ], ncol = 2,
                                            dimnames = list(NULL, c("x1", "x2"))))
  expect_lt(abs(p$EN$mean - y[4]), 0.05)
})
