skip_on_cran()
skip_if_not_installed("DiceKriging")
skip_if_not_installed("DiceOptim")
skip_if_not_installed("lhs")
skip_if_not_installed("hetGP")

toy_sim_fun <- function(theta, fidelity = c("low", "med", "high"), seed = 123, n_rep = NULL, ...) {
  fidelity <- match.arg(fidelity)
  x <- purrr::map_dbl(theta, as.numeric)
  noise <- dplyr::case_when(
    fidelity == "low" ~ 0.02,
    fidelity == "med" ~ 0.01,
    TRUE ~ 0
  )
  en <- (x[["x1"]] - 0.5)^2 + (x[["x2"]] - 0.25)^2 + noise
  power <- 0.9 - 0.6 * (x[["x1"]] - 0.5)^2
  type1 <- 0.08 + 0.3 * (x[["x2"]] - 0.25)^2
  res <- c(power = power, type1 = type1, EN = en, ET = en * 1.5)
  attr(res, "variance") <- c(power = 0.001, type1 = 0.001, EN = 0.002, ET = 0.003)
  default_rep <- dplyr::case_when(
    fidelity == "low" ~ 200,
    fidelity == "med" ~ 500,
    TRUE ~ 1000
  )
  attr(res, "n_rep") <- if (is.null(n_rep)) default_rep else n_rep
  res
}

toy_bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
toy_constraints <- list(power = c("ge", 0.8), type1 = c("le", 0.12))

test_that("bo_calibrate returns expected structure", {
  fit <- bo_calibrate(
    sim_fun = toy_sim_fun,
    bounds = toy_bounds,
    objective = "EN",
    constraints = toy_constraints,
    n_init = 4,
    q = 1,
    budget = 8,
    progress = FALSE
  )
  expect_s3_class(fit, "BATON_fit")
  expect_equal(nrow(fit$history), 8)
  expect_true(is.list(fit$best_theta))
  expect_true(length(fit$surrogates) >= 1)
})

test_that("benchmark, reliability, and ablation helpers run", {
  bench <- benchmark_methods(
    sim_fun = toy_sim_fun,
    bounds = toy_bounds,
    objective = "EN",
    constraints = toy_constraints,
    strategies = c("bo", "random"),
    bo_args = list(n_init = 3, q = 1, budget = 5, progress = FALSE, seed = 99),
    random_args = list(n_samples = 4, seeds = 1),
    simulators_per_eval = list(random = 500),
    progress = FALSE
  )
  expect_s3_class(bench, "BATON_benchmark")
  summary_tbl <- summarise_benchmark(bench)
  expect_true(nrow(summary_tbl) >= 1)

  reliability <- estimate_constraint_reliability(
    sim_fun = toy_sim_fun,
    bounds = toy_bounds,
    objective = "EN",
    constraints = toy_constraints,
    strategies = c("bo"),
    calibration_seeds = 1:2,
    validation_reps = 1000,
    bo_args = list(n_init = 3, q = 1, budget = 5, progress = FALSE),
    progress = FALSE
  )
  expect_s3_class(reliability, "BATON_reliability")
  expect_true(nrow(reliability$summary) == 1)

  ablation <- ablation_multifidelity(
    sim_fun = toy_sim_fun,
    bounds = toy_bounds,
    objective = "EN",
    constraints = toy_constraints,
    policies = list(low_only = c(low = 200), full = c(low = 200, high = 1000)),
    seeds = 1:2,
    bo_args = list(n_init = 3, q = 1, budget = 5, progress = FALSE),
    progress = FALSE
  )
  expect_s3_class(ablation, "BATON_multifidelity")
  expect_true(nrow(ablation$summary) == 2)
})

test_that("sensitivity diagnostics and case study summaries work", {
  sens <- sensitivity_diagnostics(
    sim_fun = toy_sim_fun,
    bounds = toy_bounds,
    objective = "EN",
    constraints = toy_constraints,
    kernel_options = c("matern52", "matern32"),
    acquisition_options = c("eci", "qehvi"),
    bo_args = list(n_init = 3, q = 1, budget = 6, progress = FALSE),
    sobol_samples = 200,
    gradient_points = 10,
    progress = FALSE
  )
  expect_true(is.data.frame(sens$sobol))
  expect_true(is.data.frame(sens$gradients))

  summary <- summarise_case_study(sens$baseline_fit)
  expect_true(is.data.frame(summary$design))
  expect_true(is.data.frame(summary$operating_characteristics))

  diag <- case_study_diagnostics(sens$baseline_fit, sobol_samples = 200, gradient_points = 10)
  expect_true(is.matrix(diag$covariance))
})
