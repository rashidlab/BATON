# test-parameter-importance.R
# Unit tests for parameter importance analysis functions

library(testthat)
library(DiceKriging)

test_that("extract_lengthscales works with DiceKriging models", {
  # Create mock GP with known lengthscales
  set.seed(2025)
  n <- 20
  d <- 3
  X <- matrix(runif(n * d), ncol = d)
  y <- rowSums(X^2) + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  surrogates <- list(objective = km_model)

  # Extract lengthscales
  lengthscales <- extract_lengthscales(surrogates, metric = "objective")

  # Check structure
  expect_type(lengthscales, "double")
  expect_length(lengthscales, d)
  expect_true(all(lengthscales > 0))  # Lengthscales must be positive
})


test_that("analyze_parameter_importance returns correct structure", {
  # Create mock BO result
  set.seed(2025)
  n <- 20
  d <- 3
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- c("eff", "fut", "ev")
  y <- rowSums(X^2) + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  fit <- list(
    surrogates = list(EN = km_model),
    best_theta = list(eff = 0.95, fut = 0.15, ev = 20)
  )

  # Analyze importance
  importance <- analyze_parameter_importance(fit, metric = "EN")

  # Check structure
  expect_s3_class(importance, "data.frame")
  expect_equal(nrow(importance), d)
  expect_true(all(c("parameter", "lengthscale", "importance", "rank") %in% names(importance)))

  # Check importance sums to 1
  expect_equal(sum(importance$importance), 1.0, tolerance = 1e-6)

  # Check ranks are valid
  expect_equal(sort(importance$rank), 1:d)

  # Check sorted by importance (descending)
  expect_true(all(diff(importance$importance) <= 0))
})


test_that("compute_local_sensitivity returns correct structure", {
  # Create mock BO result with simple quadratic function
  set.seed(2025)
  n <- 20
  d <- 3
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- c("eff", "fut", "ev")

  # y = x1^2 + 2*x2^2 + 0.5*x3^2 (x2 should have highest sensitivity)
  y <- X[, 1]^2 + 2 * X[, 2]^2 + 0.5 * X[, 3]^2 + rnorm(n, sd = 0.01)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  best_theta <- list(eff = 0.5, fut = 0.5, ev = 0.5)

  fit <- list(
    surrogates = list(EN = km_model),
    best_theta = best_theta,
    bounds = list(eff = c(0, 1), fut = c(0, 1), ev = c(0, 1))
  )

  # Compute sensitivity
  sensitivity <- compute_local_sensitivity(fit, epsilon = 0.01, metric = "EN")

  # Check structure
  expect_s3_class(sensitivity, "data.frame")
  expect_equal(nrow(sensitivity), d)
  expect_true(all(c("parameter", "gradient", "abs_sensitivity", "normalized_sensitivity") %in% names(sensitivity)))

  # Check normalized sensitivity sums to 1
  expect_equal(sum(sensitivity$normalized_sensitivity), 1.0, tolerance = 1e-6)

  # Check sorted by absolute sensitivity (descending)
  expect_true(all(diff(sensitivity$abs_sensitivity) <= 0))
})


test_that("compute_local_sensitivity gradients have correct signs", {
  # Create mock BO result with linear function: y = -2*x1 + 3*x2
  set.seed(2025)
  n <- 20
  d <- 2
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- c("x1", "x2")

  y <- -2 * X[, 1] + 3 * X[, 2] + rnorm(n, sd = 0.01)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  fit <- list(
    surrogates = list(objective = km_model),
    best_theta = list(x1 = 0.5, x2 = 0.5),
    bounds = list(x1 = c(0, 1), x2 = c(0, 1))
  )

  # Compute sensitivity
  sensitivity <- compute_local_sensitivity(fit, epsilon = 0.01, metric = "objective")

  # Check gradient signs (should be approximately -2 and +3)
  x1_row <- sensitivity[sensitivity$parameter == "x1", ]
  x2_row <- sensitivity[sensitivity$parameter == "x2", ]

  expect_lt(x1_row$gradient, 0)  # Should be negative
  expect_gt(x2_row$gradient, 0)  # Should be positive

  # Check magnitudes (x2 should have larger absolute gradient)
  expect_gt(abs(x2_row$gradient), abs(x1_row$gradient))
})


test_that("generate_recommendations classifies parameters correctly", {
  # Create mock BO result with known importance structure
  set.seed(2025)
  n <- 30
  d <- 5
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- c("very_important", "important", "marginal", "unimportant", "irrelevant")

  # Create function where parameters have different importance
  # very_important: high variance (coefficient = 10)
  # important: medium variance (coefficient = 5)
  # marginal: low variance (coefficient = 2)
  # unimportant: very low variance (coefficient = 0.5)
  # irrelevant: essentially zero variance (coefficient = 0.01)
  y <- 10 * X[, 1] + 5 * X[, 2] + 2 * X[, 3] + 0.5 * X[, 4] + 0.01 * X[, 5] + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  best_theta <- as.list(setNames(rep(0.5, d), colnames(X)))

  fit <- list(
    surrogates = list(EN = km_model),
    best_theta = best_theta,
    bounds = lapply(setNames(rep(list(c(0, 1)), d), colnames(X)), identity)
  )

  # Generate recommendations
  recs <- generate_recommendations(
    fit,
    threshold_irrelevant = 0.01,
    threshold_unimportant = 0.10
  )

  # Check structure
  expect_s3_class(recs, "param_recommendations")
  expect_true(all(c("analysis", "recommendations", "thresholds") %in% names(recs)))

  # Check analysis data frame
  expect_s3_class(recs$analysis, "data.frame")
  expect_equal(nrow(recs$analysis), d)
  expect_true(all(c("parameter", "consensus_importance", "action") %in% names(recs$analysis)))

  # Check recommendations list
  expect_equal(length(recs$recommendations), d)
  expect_true(all(sapply(recs$recommendations, function(r) "action" %in% names(r))))

  # Check thresholds recorded
  expect_equal(recs$thresholds$irrelevant, 0.01)
  expect_equal(recs$thresholds$unimportant, 0.10)

  # Count actions
  actions <- sapply(recs$recommendations, function(r) r$action)
  n_optimize <- sum(actions == "OPTIMIZE")
  n_fix <- sum(actions == "FIX")
  n_remove <- sum(actions == "REMOVE")

  # Should have at least 2 parameters to optimize (very_important, important)
  expect_gte(n_optimize, 2)

  # Should have at least 1 parameter to fix or remove (unimportant, irrelevant)
  expect_gte(n_fix + n_remove, 1)
})


test_that("generate_recommendations fixes parameters at best values", {
  # Create simple mock with explicit column names matching best_theta
  set.seed(2025)
  n <- 20
  d <- 3
  X <- matrix(runif(n * d), ncol = d)
  colnames(X) <- c("x1", "x2", "x3")

  # important has high impact, others have very low impact
  y <- 100 * X[, 1] + 0.01 * X[, 2] + 0.001 * X[, 3] + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  # Use parameter names that match what extract_lengthscales returns
  best_theta <- list(x1 = 0.95, x2 = 0.42, x3 = 0.13)

  fit <- list(
    surrogates = list(EN = km_model),
    best_theta = best_theta,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1), x3 = c(0, 1))
  )

  # Generate recommendations with lenient thresholds
  recs <- generate_recommendations(
    fit,
    threshold_irrelevant = 0.001,
    threshold_unimportant = 0.30  # More lenient to ensure some params are fixed
  )

  # Check that FIX recommendations include the best value
  fix_count <- 0
  for (param in names(recs$recommendations)) {
    rec <- recs$recommendations[[param]]
    if (rec$action == "FIX") {
      fix_count <- fix_count + 1
      expect_equal(rec$value, best_theta[[param]])
      expect_true(grepl("Fix at best value", rec$reason))
    }
  }

  # Should have at least one parameter fixed
  expect_gte(fix_count, 1)
})


test_that("print.param_recommendations runs without error", {
  # Create minimal mock recommendation object
  recs <- structure(
    list(
      analysis = data.frame(
        parameter = c("x1", "x2", "x3"),
        consensus_importance = c(0.60, 0.35, 0.05),
        action = c("OPTIMIZE", "OPTIMIZE", "FIX"),
        stringsAsFactors = FALSE
      ),
      recommendations = list(
        x1 = list(action = "OPTIMIZE", reason = "Important parameter (60.0%). Continue optimizing."),
        x2 = list(action = "OPTIMIZE", reason = "Important parameter (35.0%). Continue optimizing."),
        x3 = list(action = "FIX", value = 0.123, reason = "Low importance (5.0%). Fix at best value: 0.123")
      ),
      thresholds = list(irrelevant = 0.01, unimportant = 0.05)
    ),
    class = "param_recommendations"
  )

  # Print should work without error and return the object invisibly
  result <- capture.output(returned <- print(recs))
  expect_identical(returned, recs)  # Should return input invisibly

  # Verify output contains expected text
  expect_true(length(result) > 0)  # Should produce some output

  # Verify structure is correct
  expect_equal(class(recs), "param_recommendations")
  expect_equal(length(recs$recommendations), 3)
})


test_that("extract_lengthscales handles missing metric gracefully", {
  # Create mock surrogates
  set.seed(2025)
  n <- 20
  d <- 2
  X <- matrix(runif(n * d), ncol = d)
  y <- rowSums(X) + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  surrogates <- list(EN = km_model, power = km_model)

  # Should work with default (first metric)
  lengthscales_default <- extract_lengthscales(surrogates)
  expect_length(lengthscales_default, d)

  # Should work with explicit metric
  lengthscales_power <- extract_lengthscales(surrogates, metric = "power")
  expect_length(lengthscales_power, d)

  # Should error on invalid metric
  expect_error(
    extract_lengthscales(surrogates, metric = "nonexistent"),
    "not found in surrogates"
  )
})


test_that("analyze_parameter_importance errors on missing surrogates", {
  fit_no_surrogates <- list(
    best_theta = list(x1 = 0.5, x2 = 0.5)
  )

  expect_error(
    analyze_parameter_importance(fit_no_surrogates),
    "does not contain surrogates"
  )
})


test_that("compute_local_sensitivity errors on missing best_theta", {
  set.seed(2025)
  n <- 20
  d <- 2
  X <- matrix(runif(n * d), ncol = d)
  y <- rowSums(X) + rnorm(n, sd = 0.1)

  km_model <- DiceKriging::km(
    formula = ~1,
    design = data.frame(X),
    response = y,
    covtype = "gauss",
    control = list(trace = FALSE)
  )

  fit_no_best <- list(surrogates = list(EN = km_model))

  expect_error(
    compute_local_sensitivity(fit_no_best),
    "does not contain best_theta"
  )
})
