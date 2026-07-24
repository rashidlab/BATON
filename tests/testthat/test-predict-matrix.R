# Task C5: predict_surrogates accepts a candidate matrix directly, and
# invalid candidates error instead of being silently dropped (a drop would
# misalign predictions with the caller's candidate indices).

c5_surrogates <- function() {
  set.seed(21)
  n <- 30
  X <- matrix(runif(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  y <- (X[, 1] - 0.5)^2 + (X[, 2] - 0.25)^2
  history <- tibble::tibble(
    unit_x = lapply(seq_len(n), function(i) X[i, ]),
    theta_id = apply(X, 1, BATON:::theta_to_id),
    metrics = lapply(y, function(v) c(EN = v)),
    variance = replicate(n, c(EN = NA_real_), simplify = FALSE)
  )
  fit_surrogates(history, "EN", BATON:::parse_constraints(list()),
                 use_hetgp = FALSE, fit_seed = 5)
}

test_that("matrix and list inputs give identical predictions", {
  skip_if_not_installed("DiceKriging")
  s <- c5_surrogates()
  set.seed(9)
  M <- matrix(runif(40), ncol = 2)
  colnames(M) <- c("x1", "x2")
  as_list <- lapply(seq_len(nrow(M)), function(i) M[i, ])

  p_mat <- BATON:::predict_surrogates(s, M)
  p_list <- BATON:::predict_surrogates(s, as_list)
  expect_equal(p_mat$EN$mean, p_list$EN$mean)
  expect_equal(p_mat$EN$sd, p_list$EN$sd)
})

test_that("unnamed matrix columns are accepted in parameter order", {
  skip_if_not_installed("DiceKriging")
  s <- c5_surrogates()
  set.seed(10)
  M <- matrix(runif(20), ncol = 2)
  M_named <- M; colnames(M_named) <- c("x1", "x2")
  expect_equal(BATON:::predict_surrogates(s, M)$EN$mean,
               BATON:::predict_surrogates(s, M_named)$EN$mean)
})

test_that("shuffled matrix columns are reordered by name", {
  skip_if_not_installed("DiceKriging")
  s <- c5_surrogates()
  set.seed(11)
  M <- matrix(runif(20), ncol = 2)
  colnames(M) <- c("x1", "x2")
  M_shuffled <- M[, c("x2", "x1")]
  expect_equal(BATON:::predict_surrogates(s, M)$EN$mean,
               BATON:::predict_surrogates(s, M_shuffled)$EN$mean)
})

test_that("invalid candidates error instead of silently dropping", {
  skip_if_not_installed("DiceKriging")
  s <- c5_surrogates()
  # list path: a candidate missing a parameter
  bad_list <- list(c(x1 = 0.5, x2 = 0.5), c(x1 = 0.2))
  expect_error(BATON:::predict_surrogates(s, bad_list), "invalid candidate")
  # matrix path: non-finite coordinates
  M <- matrix(c(0.5, 0.5, NA, 0.2), ncol = 2, byrow = TRUE)
  colnames(M) <- c("x1", "x2")
  expect_error(BATON:::predict_surrogates(s, M), "non-finite")
  # list path: correctly-named but non-finite coordinates must ALSO error,
  # not silently degrade to the constant-predictor fallback
  bad_na <- list(c(x1 = 0.5, x2 = 0.5), c(x1 = NA_real_, x2 = 0.2))
  expect_error(BATON:::predict_surrogates(s, bad_na), "non-finite")
  bad_inf <- list(c(x1 = 0.5, x2 = 0.5), c(x1 = Inf, x2 = 0.2))
  expect_error(BATON:::predict_surrogates(s, bad_inf), "non-finite")
  # matrix path: wrong column count
  M3 <- matrix(runif(9), ncol = 3)
  expect_error(BATON:::predict_surrogates(s, M3), "expected 2")
})

test_that("select_batch_local_penalization accepts matrix candidates", {
  cand_mat <- matrix(c(0.00, 0.01, 0.99), ncol = 1)
  cand_list <- list(0.00, 0.01, 0.99)
  scores <- c(1.0, 0.95, 0.80)
  sel_mat <- BATON:::select_batch_local_penalization(cand_mat, scores, q = 2,
                                                     lipschitz = 10)
  sel_list <- BATON:::select_batch_local_penalization(cand_list, scores, q = 2,
                                                      lipschitz = 10)
  expect_equal(sort(sel_mat), sort(sel_list))
  expect_equal(sort(sel_mat), c(1, 3))
})
