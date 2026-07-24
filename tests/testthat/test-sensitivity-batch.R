# Task C8: batched sensitivity/diagnostics predicts. These are
# characterization tests: the reference implementations below reproduce the
# original per-point/per-call logic, and the batched code must match them
# within numerical tolerance.

c8_surrogates <- function() {
  set.seed(41)
  n <- 25
  X <- matrix(runif(n * 2), ncol = 2)
  colnames(X) <- c("x1", "x2")
  # strongly x1-driven surface so Sobol ordering is unambiguous
  y <- 3 * (X[, 1] - 0.5)^2 + 0.1 * X[, 2] + rnorm(n, sd = 0.01)
  history <- tibble::tibble(
    unit_x = lapply(seq_len(n), function(i) X[i, ]),
    theta_id = apply(X, 1, BATON:::theta_to_id),
    metrics = lapply(y, function(v) c(EN = v)),
    variance = replicate(n, c(EN = 1e-4), simplify = FALSE)
  )
  fit_surrogates(history, "EN", BATON:::parse_constraints(list()),
                 use_hetgp = FALSE, fit_seed = 8)
}

c8_bounds <- list(x1 = c(0, 2), x2 = c(-1, 1))

# Reference: the original one-predict-per-(point, param) finite difference.
ref_gradient <- function(surrogates, theta_unit, bounds, param, eps, outcome) {
  up <- down <- theta_unit
  up[param] <- min(1, up[param] + eps)
  down[param] <- max(0, down[param] - eps)
  if (abs(up[param] - down[param]) < 1e-8) return(list(gradient = 0, sd = 0))
  preds <- BATON:::predict_surrogates(surrogates[outcome], list(up, down))[[outcome]]
  theta_up <- BATON:::scale_from_unit(as.list(up), bounds)
  theta_down <- BATON:::scale_from_unit(as.list(down), bounds)
  delta <- as.numeric(theta_up[[param]] - theta_down[[param]])
  list(gradient = (preds$mean[1] - preds$mean[2]) / delta,
       sd = sqrt(sum(preds$sd^2)) / delta)
}

test_that("sa_gradients matches the per-point reference (interior point)", {
  skip_if_not_installed("DiceKriging")
  s <- c8_surrogates()
  theta <- list(x1 = 1.2, x2 = 0.4)
  theta_unit <- BATON:::scale_to_unit(theta, c8_bounds)
  got <- sa_gradients(s, theta, c8_bounds, outcome = "EN", eps = 1e-4)
  for (p in names(c8_bounds)) {
    ref <- ref_gradient(s, theta_unit, c8_bounds, p, 1e-4, "EN")
    expect_equal(got$gradient[got$parameter == p], ref$gradient, tolerance = 1e-8)
    expect_equal(got$sd[got$parameter == p], ref$sd, tolerance = 1e-8)
  }
})

test_that("sa_gradients matches the reference at a clamped boundary point", {
  skip_if_not_installed("DiceKriging")
  s <- c8_surrogates()
  theta <- list(x1 = 0, x2 = 1)  # unit coords 0 and 1: up/down get clamped
  theta_unit <- BATON:::scale_to_unit(theta, c8_bounds)
  got <- sa_gradients(s, theta, c8_bounds, outcome = "EN", eps = 1e-4)
  for (p in names(c8_bounds)) {
    ref <- ref_gradient(s, theta_unit, c8_bounds, p, 1e-4, "EN")
    expect_equal(got$gradient[got$parameter == p], ref$gradient, tolerance = 1e-8)
  }
})

test_that("cov_effects matches per-point reference gradients over the same draw", {
  skip_if_not_installed("DiceKriging")
  s <- c8_surrogates()
  set.seed(77)
  got <- cov_effects(s, c8_bounds, outcome = "EN", n_mc = 40, eps = 1e-3)
  set.seed(77)
  pts <- lhs::randomLHS(40, 2)
  colnames(pts) <- names(c8_bounds)
  ref_mat <- matrix(NA_real_, nrow = 40, ncol = 2)
  for (i in seq_len(40)) {
    tu <- pts[i, , drop = TRUE]
    for (j in seq_len(2)) {
      ref_mat[i, j] <- ref_gradient(s, tu, c8_bounds, names(c8_bounds)[j],
                                    1e-3, "EN")$gradient
    }
  }
  ref_cov <- stats::cov(ref_mat)
  dimnames(ref_cov) <- list(names(c8_bounds), names(c8_bounds))
  expect_equal(got, ref_cov, tolerance = 1e-6)
})

test_that("sa_sobol identifies the dominant parameter with valid indices", {
  skip_if_not_installed("DiceKriging")
  s <- c8_surrogates()
  set.seed(3)
  sob <- sa_sobol(s, c8_bounds, outcome = "EN", n_mc = 400)
  expect_true(all(sob$S_first >= 0 & sob$S_first <= 1))
  expect_gt(sob$S_first[sob$parameter == "x1"],
            sob$S_first[sob$parameter == "x2"])
})
