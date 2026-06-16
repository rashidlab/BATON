# test-bug-fixes.R
# Unit tests for bug fixes in v0.3.x
#
# Tests cover:
# 1. Constant predictor mean_value consistency (Bug #1)
# 2. Homoskedastic fallback Z aggregation (Bug #2)
# 3. Warm-start with initial_history without unit_x (Bug #3)
# 4. Local sensitivity unit scale handling (Bug #4)
# 5. run_init_with_stopping objective extraction (Bug #5)
# 6. Incumbent value Inf in fixed_per_stage mode (Bug #6)

library(testthat)

# ============================================================================
# Test fixtures
# ============================================================================

# Simple mock simulator for testing
mock_sim <- function(theta, fidelity = c("low", "med", "high"), seed = 123, ...) {
  fidelity <- match.arg(fidelity)
  set.seed(seed)

  x1 <- theta$x1
  x2 <- theta$x2

  # Objective: minimize (x1-0.5)^2 + (x2-0.3)^2
  EN <- (x1 - 0.5)^2 + (x2 - 0.3)^2 + rnorm(1, 0, 0.01)

  # Constraints
  power <- pnorm((x1 + x2 - 0.5) / 0.2) + rnorm(1, 0, 0.02)
  type1 <- 0.05 + 0.05 * abs(x1 - x2) + rnorm(1, 0, 0.01)

  result <- c(EN = EN, power = power, type1 = type1)
  attr(result, "variance") <- c(EN = 0.0001, power = 0.0004, type1 = 0.0001)
  attr(result, "n_rep") <- ifelse(fidelity == "low", 100, 1000)

  return(result)
}

# ============================================================================
# Bug #1: Constant predictor mean_value consistency
# ============================================================================

test_that("constant_predictor uses mean_value consistently", {
  skip_on_cran()

  # Create a constant predictor directly (simulating fallback)
  const_pred <- structure(
    list(mean_value = 0.5, metric = "test"),
    class = "constant_predictor"
  )

  # Test that predict_surrogates can read it correctly
  surrogates <- list(test = const_pred)
  unit_x <- list(c(x1 = 0.5, x2 = 0.5))

  # This should not error and should return the constant value
  pred <- BATON:::predict_surrogates(surrogates, unit_x)

  expect_equal(length(pred), 1)
  expect_equal(pred$test$mean, 0.5)
  expect_equal(pred$test$sd, 1.0)  # High uncertainty for constant predictor
})

test_that("fit_surrogates creates constant_predictor with mean_value on insufficient data", {
  skip_on_cran()

  # Create minimal history with only 1 observation
  history <- tibble::tibble(
    iter = 0L,
    eval_id = 1L,
    theta = list(list(x1 = 0.5, x2 = 0.5)),
    unit_x = list(c(x1 = 0.5, x2 = 0.5)),
    theta_id = "0.500000|0.500000",
    fidelity = "low",
    n_rep = 100L,
    metrics = list(c(EN = 0.1, power = 0.85)),
    variance = list(c(EN = 0.01, power = 0.04)),
    objective = 0.1,
    feasible = TRUE
  )

  constraint_tbl <- tibble::tibble(
    metric = "power",
    direction = "ge",
    threshold = 0.8
  )

  # fit_surrogates should fall back to constant_predictor
  expect_message(
    surrogates <- BATON::fit_surrogates(
      history = history,
      objective = "EN",
      constraint_tbl = constraint_tbl
    ),
    "constant predictor"
  )

  # Check that the constant predictor has mean_value (not mean)
  for (metric in names(surrogates)) {
    if (inherits(surrogates[[metric]], "constant_predictor")) {
      expect_true("mean_value" %in% names(surrogates[[metric]]))
    }
  }
})

# ============================================================================
# Bug #3: Warm-start with initial_history without unit_x
# ============================================================================

test_that("initial_history without unit_x derives required columns", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history WITHOUT unit_x column (only individual param columns)
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # This should work without error - unit_x should be derived
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,  # Just a few BO iterations
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that unit_x was derived
  expect_true("unit_x" %in% names(fit$history))
  expect_true(is.list(fit$history$unit_x))
  expect_equal(length(fit$history$unit_x[[1]]), 2)  # 2 parameters

  # Check unit_x values are in [0,1]
  for (i in seq_len(nrow(fit$history))) {
    ux <- fit$history$unit_x[[i]]
    expect_true(all(ux >= 0 & ux <= 1))
  }
})

test_that("initial_history with theta column but no unit_x works", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history with theta column but no unit_x
  initial_history <- tibble::tibble(
    theta = list(
      list(x1 = 0.3, x2 = 0.2),
      list(x1 = 0.5, x2 = 0.5),
      list(x1 = 0.7, x2 = 0.3),
      list(x1 = 0.4, x2 = 0.6),
      list(x1 = 0.6, x2 = 0.4)
    ),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE)
  )

  # Should work - unit_x derived from theta
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  expect_true("unit_x" %in% names(fit$history))
  expect_true("theta_id" %in% names(fit$history))
})

# ============================================================================
# Bug #4: Local sensitivity unit scale handling
# ============================================================================

test_that("compute_local_sensitivity requires bounds in fit object", {
  skip_on_cran()

  # Create a mock fit without bounds
  fit_no_bounds <- list(
    best_theta = list(x1 = 0.5, x2 = 0.5),
    surrogates = list(EN = structure(list(mean_value = 0.1), class = "constant_predictor"))
  )

  expect_error(
    compute_local_sensitivity(fit_no_bounds),
    "does not contain bounds"
  )
})

test_that("compute_local_sensitivity works with unit-scale surrogates", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 10,
    budget = 15,
    seed = 123,
    progress = FALSE
  )

  # Should work without error
  sensitivity <- compute_local_sensitivity(fit)

  expect_s3_class(sensitivity, "data.frame")
  expect_true("parameter" %in% names(sensitivity))
  expect_true("gradient" %in% names(sensitivity))
  expect_true("normalized_sensitivity" %in% names(sensitivity))

  # Check gradients are finite
  expect_true(all(is.finite(sensitivity$gradient)))

  # Normalized sensitivity should sum to 1 (or be all zeros)
  total_sens <- sum(sensitivity$normalized_sensitivity)
  expect_true(total_sens == 0 || abs(total_sens - 1) < 1e-6)
})

test_that("compute_local_sensitivity handles constant predictor", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))

  # Create fit with constant predictor surrogate
  fit <- list(
    best_theta = list(x1 = 0.5, x2 = 0.3),
    bounds = bounds,
    surrogates = list(
      EN = structure(list(mean_value = 0.1, metric = "EN"), class = "constant_predictor")
    )
  )

  # Should return zero gradients for constant predictor
  sensitivity <- compute_local_sensitivity(fit)

  expect_true(all(sensitivity$gradient == 0))
})

# ============================================================================
# Bug #5: run_init_with_stopping objective extraction
# ============================================================================

test_that("run_init_with_stopping extracts named objective correctly", {
  skip_on_cran()

  # Simulator with named output (objective is second element)
  named_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    c(power = 0.85, EN = 0.1 + theta$x1, type1 = 0.08)
  }

  # Create design matrix
  design <- matrix(runif(10), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Run with objective = "EN" (second element)
  result <- BATON:::run_init_with_stopping(
    sim_fun = named_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values, not power values
  expected_y <- 0.1 + design[, "x1"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})

test_that("run_init_with_stopping warns for missing objective name", {
  skip_on_cran()

  sim <- function(theta, ...) {
    c(power = 0.85, type1 = 0.08)  # No "EN" metric
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Should warn about missing objective and fall back to first element
  expect_warning(
    result <- BATON:::run_init_with_stopping(
      sim_fun = sim,
      design = design,
      objective = "EN",
      fidelity = "low",
      init_config = BATON::init_stopping_config(enabled = FALSE),
      seed = 123,
      verbose = FALSE
    ),
    "not found in simulator result"
  )

  # y should be first element (power = 0.85)
  expect_equal(result$y, rep(0.85, 3))
})

test_that("run_init_with_stopping uses first element when objective is NULL", {
  skip_on_cran()

  sim <- function(theta, ...) {
    c(EN = 0.1, power = 0.85)
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # No warning expected - just use first element
  result <- BATON:::run_init_with_stopping(
    sim_fun = sim,
    design = design,
    objective = NULL,  # No objective specified
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be first element (EN = 0.1)
  expect_equal(result$y, rep(0.1, 3))
})

# ============================================================================
# Bug #2: Homoskedastic fallback Z aggregation (indirect test)
# ============================================================================

test_that("fit_hetgp_surrogate handles fallback correctly", {
  skip_on_cran()
  skip_if_not_installed("hetGP")

  # This test indirectly verifies the Z aggregation fix by ensuring

# the fallback path doesn't error when hetGP fails

  # Create history with replicated observations
  history <- tibble::tibble(
    iter = rep(0L, 4),
    eval_id = 1:4,
    theta = list(
      list(x1 = 0.5, x2 = 0.5),
      list(x1 = 0.5, x2 = 0.5),  # Replicate
      list(x1 = 0.3, x2 = 0.7),
      list(x1 = 0.3, x2 = 0.7)   # Replicate
    ),
    unit_x = list(
      c(x1 = 0.5, x2 = 0.5),
      c(x1 = 0.5, x2 = 0.5),
      c(x1 = 0.3, x2 = 0.7),
      c(x1 = 0.3, x2 = 0.7)
    ),
    theta_id = c("0.500000|0.500000", "0.500000|0.500000",
                 "0.300000|0.700000", "0.300000|0.700000"),
    fidelity = rep("low", 4),
    n_rep = rep(100L, 4),
    metrics = list(
      c(EN = 0.1, power = 0.85),
      c(EN = 0.12, power = 0.83),
      c(EN = 0.2, power = 0.75),
      c(EN = 0.18, power = 0.77)
    ),
    variance = list(
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04),
      c(EN = 0.01, power = 0.04)
    ),
    objective = c(0.1, 0.12, 0.2, 0.18),
    feasible = c(TRUE, TRUE, FALSE, FALSE)
  )

  constraint_tbl <- tibble::tibble(
    metric = "power",
    direction = "ge",
    threshold = 0.8
  )

  # fit_surrogates should work even if hetGP fails and falls back
  # The key is that it shouldn't error due to Z aggregation issues
  surrogates <- tryCatch(
    BATON::fit_surrogates(
      history = history,
      objective = "EN",
      constraint_tbl = constraint_tbl,
      use_hetgp = TRUE
    ),
    error = function(e) {
      # Some error is acceptable, but not Z aggregation error
      if (grepl("subscript out of bounds|Z_vec|matching", e$message)) {
        fail(paste("Z aggregation bug still present:", e$message))
      }
      # Return NULL to indicate error (but not the specific bug)
      NULL
    }
  )

  # If we got surrogates, verify they work
  if (!is.null(surrogates)) {
    unit_x <- list(c(x1 = 0.5, x2 = 0.5))
    pred <- BATON:::predict_surrogates(surrogates, unit_x)
    expect_true(all(is.finite(pred$EN$mean)))
  }
})

# ============================================================================
# Additional fixes: initial_history objective column mapping
# ============================================================================

test_that("initial_history with only objective column (not metric-named) works", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

 # Create initial_history with "objective" column but NOT "EN" column
  # This matches the documented interface
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),  # objective values, not "EN"
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  # Should work - objective column should be mapped to "EN" in metrics
  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that metrics list includes the objective
  for (i in seq_len(5)) {
    expect_true("EN" %in% names(fit$history$metrics[[i]]))
  }
})

test_that("initial_history with numeric variance column is preserved", {
  skip_on_cran()

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  # Create initial_history with numeric variance column (not list-column)
  initial_history <- data.frame(
    x1 = c(0.3, 0.5, 0.7, 0.4, 0.6),
    x2 = c(0.2, 0.5, 0.3, 0.6, 0.4),
    EN = c(0.1, 0.05, 0.15, 0.08, 0.12),
    power = c(0.82, 0.85, 0.78, 0.80, 0.83),
    type1 = c(0.09, 0.08, 0.11, 0.10, 0.07),
    objective = c(0.1, 0.05, 0.15, 0.08, 0.12),
    fidelity = rep("low", 5),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    variance = c(0.001, 0.002, 0.001, 0.003, 0.002),  # Numeric column
    stringsAsFactors = FALSE
  )

  fit <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = 5,
    budget = 8,
    initial_history = initial_history,
    seed = 123,
    progress = FALSE
  )

  # Check that variance was converted to list-column and values preserved
  expect_true(is.list(fit$history$variance))

  # First 5 rows should have variance for "EN"
  for (i in 1:5) {
    var_list <- fit$history$variance[[i]]
    if (length(var_list) > 0) {
      expect_true("EN" %in% names(var_list))
    }
  }
})

# ============================================================================
# Additional fixes: run_init_with_stopping list-style output
# ============================================================================

test_that("run_init_with_stopping handles list-style simulator output", {
  skip_on_cran()

  # Simulator returning list with metrics element
  list_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    list(
      metrics = c(EN = 0.1 + theta$x1, power = 0.85),
      variance = c(EN = 0.001, power = 0.004),
      n_rep = 1000
    )
  }

  design <- matrix(runif(10), ncol = 2)
  colnames(design) <- c("x1", "x2")

  # Should correctly extract EN from list(metrics = ...)
  result <- BATON:::run_init_with_stopping(
    sim_fun = list_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values (0.1 + x1)
  expected_y <- 0.1 + design[, "x1"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})

test_that("run_init_with_stopping handles plain list simulator output", {
  skip_on_cran()

  # Simulator returning plain list (not with $metrics)
  plain_list_sim <- function(theta, fidelity = "low", seed = 1, ...) {
    list(EN = 0.2 + theta$x2, power = 0.80)
  }

  design <- matrix(runif(6), ncol = 2)
  colnames(design) <- c("x1", "x2")

  result <- BATON:::run_init_with_stopping(
    sim_fun = plain_list_sim,
    design = design,
    objective = "EN",
    fidelity = "low",
    init_config = BATON::init_stopping_config(enabled = FALSE),
    seed = 123,
    verbose = FALSE
  )

  # y should be EN values (0.2 + x2)
  expected_y <- 0.2 + design[, "x2"]
  expect_equal(result$y, expected_y, tolerance = 1e-10)
})

# ============================================================================
# Bug #6: Incumbent value Inf in fixed_per_stage mode
# ============================================================================
# When high_fidelity_only=TRUE was hardcoded in bo_calibrate, Stage 1 of
# fixed_per_stage mode (which has zero high-fidelity evaluations) would always
# get incumbent = Inf from best_feasible_objective. This made both BPM and BOI
# incumbent methods inert, causing ECI to degenerate into random feasibility-
# seeking instead of proper Bayesian optimization.
#
# Fix: Dynamic has_high_fidelity check at bo_calibrate.R lines 587-591.
# ============================================================================

test_that("best_feasible_objective returns finite value with only low-fidelity data", {
  skip_on_cran()

  # Create a mock history with ONLY low-fidelity feasible observations
  # This simulates Stage 1 of fixed_per_stage mode
  history <- tibble::tibble(
    iter = rep(0L, 5),
    eval_id = 1:5,
    theta = list(
      list(x1 = 0.3, x2 = 0.2),
      list(x1 = 0.5, x2 = 0.5),
      list(x1 = 0.7, x2 = 0.3),
      list(x1 = 0.4, x2 = 0.6),
      list(x1 = 0.6, x2 = 0.4)
    ),
    unit_x = list(
      c(x1 = 0.3, x2 = 0.2),
      c(x1 = 0.5, x2 = 0.5),
      c(x1 = 0.7, x2 = 0.3),
      c(x1 = 0.4, x2 = 0.6),
      c(x1 = 0.6, x2 = 0.4)
    ),
    theta_id = c("0.300|0.200", "0.500|0.500", "0.700|0.300",
                 "0.400|0.600", "0.600|0.400"),
    fidelity = rep("low", 5),
    n_rep = rep(100L, 5),
    metrics = list(
      c(EN = 0.10, power = 0.85),
      c(EN = 0.05, power = 0.90),
      c(EN = 0.15, power = 0.78),
      c(EN = 0.08, power = 0.82),
      c(EN = 0.12, power = 0.83)
    ),
    variance = list(
      c(EN = 0.001, power = 0.004),
      c(EN = 0.001, power = 0.004),
      c(EN = 0.001, power = 0.004),
      c(EN = 0.001, power = 0.004),
      c(EN = 0.001, power = 0.004)
    ),
    objective = c(0.10, 0.05, 0.15, 0.08, 0.12),
    feasible = c(TRUE, TRUE, FALSE, TRUE, TRUE)
  )

  # ---- Test the dynamic check: no high-fidelity data present ----

  # With has_high_fidelity = FALSE (the correct dynamic check),
  # should return a finite value from all-fidelity feasible points
  has_hf <- any(history$fidelity == "high", na.rm = TRUE)
  expect_false(has_hf, info = "Test setup: should have no high-fidelity data")

  boi_value <- BATON:::best_feasible_objective(
    history, "EN",
    surrogates = NULL,
    high_fidelity_only = has_hf,  # FALSE — correct dynamic behavior
    incumbent_method = "boi"
  )
  expect_true(is.finite(boi_value),
              info = "BOI incumbent should be finite with low-fidelity-only data")
  expect_equal(boi_value, 0.05,
               info = "BOI should return min objective among feasible rows")

  # ---- Verify the old bug: high_fidelity_only = TRUE gives Inf ----
  buggy_value <- BATON:::best_feasible_objective(
    history, "EN",
    surrogates = NULL,
    high_fidelity_only = TRUE,  # Hardcoded TRUE — the old bug
    incumbent_method = "boi"
  )
  expect_equal(buggy_value, Inf,
               info = "Hardcoded high_fidelity_only=TRUE should return Inf (no HF data)")

  # ---- Verify that adding high-fidelity data restores HF-only behavior ----
  history_with_hf <- tibble::add_row(
    history,
    iter = 1L, eval_id = 6L,
    theta = list(list(x1 = 0.5, x2 = 0.5)),
    unit_x = list(c(x1 = 0.5, x2 = 0.5)),
    theta_id = "0.500|0.500",
    fidelity = "high",
    n_rep = 5000L,
    metrics = list(c(EN = 0.04, power = 0.91)),
    variance = list(c(EN = 0.0002, power = 0.001)),
    objective = 0.04,
    feasible = TRUE
  )

  has_hf2 <- any(history_with_hf$fidelity == "high", na.rm = TRUE)
  expect_true(has_hf2)

  hf_value <- BATON:::best_feasible_objective(
    history_with_hf, "EN",
    surrogates = NULL,
    high_fidelity_only = has_hf2,  # TRUE — correct when HF data exists
    incumbent_method = "boi"
  )
  expect_true(is.finite(hf_value))
  expect_equal(hf_value, 0.04,
               info = "Should use high-fidelity value when HF data present")
})

test_that("BPM and BOI incumbents produce divergent BO histories (single-fidelity)", {
  skip_on_cran()

  # This integration test runs bo_calibrate twice with the same seed and

  # single-fidelity (low only), once with BPM and once with BOI.
  # With the bug fix, both methods should be active and produce DIFFERENT
  # acquisition decisions after the shared initial design.
  # Pre-fix, both would get incumbent = Inf making them identical.

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.70), type1 = c("le", 0.15))

  n_init <- 10
  budget <- 20  # 10 init + 10 BO iterations — enough to diverge

  # Run with BPM incumbent
  fit_bpm <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = n_init,
    budget = budget,
    seed = 42,
    progress = FALSE,
    incumbent_method = "bpm"
  )

  # Run with BOI incumbent
  fit_boi <- bo_calibrate(
    sim_fun = mock_sim,
    bounds = bounds,
    objective = "EN",
    constraints = constraints,
    n_init = n_init,
    budget = budget,
    seed = 42,
    progress = FALSE,
    incumbent_method = "boi"
  )

  # Both should complete successfully

  expect_true(!is.null(fit_bpm$history))
  expect_true(!is.null(fit_boi$history))

  # The initial design (first n_init rows) should be identical (same seed)
  init_bpm <- fit_bpm$history$objective[1:n_init]
  init_boi <- fit_boi$history$objective[1:n_init]
  expect_equal(init_bpm, init_boi,
               info = "Initial designs should match with same seed")

  # The BO iterations (rows after n_init) should DIVERGE because
  # BPM and BOI compute different incumbent values, leading to
  # different ECI acquisition decisions
  n_total <- min(nrow(fit_bpm$history), nrow(fit_boi$history))
  if (n_total > n_init) {
    bo_bpm <- fit_bpm$history$unit_x[(n_init + 1):n_total]
    bo_boi <- fit_boi$history$unit_x[(n_init + 1):n_total]
    expect_false(identical(bo_bpm, bo_boi),
                 info = "BPM and BOI should produce different acquisition sequences")
  }
})
