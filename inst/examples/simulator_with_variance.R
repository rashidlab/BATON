# Example: Clinical Trial Simulator with Proper Variance Estimation
#
# This example demonstrates how to write a simulator that returns variance
# estimates using Welford's algorithm for memory efficiency.

library(BATON)

# =============================================================================
# Example 1: Simple simulator using welford_mean_var
# =============================================================================

#' Adaptive trial simulator with variance
#'
#' This simulator uses BATON's welford_mean_var() helper to efficiently
#' compute mean and variance without storing all samples.
#'
#' @param theta Named list of design parameters (e.g., c(threshold = 2.0, alpha = 0.025))
#' @param fidelity One of "low", "med", "high" controlling simulation replications
#' @param seed Random seed for reproducibility
#' @return Named vector of metrics with variance attribute
simple_adaptive_simulator <- function(theta,
                                      fidelity = c("low", "med", "high"),
                                      seed = NULL,
                                      ...) {
  fidelity <- match.arg(fidelity)

  # Determine number of replications based on fidelity
  n_rep <- switch(fidelity,
                  low = 200,
                  med = 1000,
                  high = 10000)

  if (!is.null(seed)) set.seed(seed)

  # Define function that simulates ONE trial
  simulate_one_trial <- function(i, theta) {
    # Extract parameters
    threshold <- theta$threshold
    alpha <- theta$alpha

    # Simulate patient enrollment (Poisson process)
    enrollment_rate <- 10  # patients per month
    max_time <- 24  # months
    n_patients <- rpois(1, enrollment_rate * max_time)

    # Simulate treatment effects
    control_mean <- 0
    treatment_effect <- 0.5
    sigma <- 1

    control_outcomes <- rnorm(n_patients / 2, control_mean, sigma)
    treatment_outcomes <- rnorm(n_patients / 2, control_mean + treatment_effect, sigma)

    # Conduct test
    test_result <- t.test(treatment_outcomes, control_outcomes)
    reject_null <- test_result$p.value < alpha

    # Calculate trial metrics
    c(
      power = as.numeric(reject_null),  # Indicator of rejection
      type1 = 0,  # Would need null scenario to estimate this
      EN = n_patients,
      ET = max_time
    )
  }

  # Use Welford's algorithm to compute mean and variance efficiently
  result <- welford_mean_var(
    sample_fn = simulate_one_trial,
    n_samples = n_rep,
    theta = theta
  )

  # Return metrics with variance attribute
  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- result$n

  return(metrics)
}


# =============================================================================
# Example 2: Parallel simulator using pool_welford_results
# =============================================================================

#' Adaptive trial simulator with parallel execution
#'
#' For computationally expensive simulators, run chunks in parallel and
#' pool the results.
parallel_adaptive_simulator <- function(theta,
                                        fidelity = c("low", "med", "high"),
                                        seed = NULL,
                                        n_cores = 4,
                                        ...) {
  fidelity <- match.arg(fidelity)

  n_rep <- switch(fidelity,
                  low = 200,
                  med = 1000,
                  high = 10000)

  # Ensure reproducibility with parallel execution
  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  # Define one trial simulation
  simulate_one_trial <- function(i, theta) {
    # Same as above...
    threshold <- theta$threshold
    alpha <- theta$alpha
    n_patients <- rpois(1, 10 * 24)
    treatment_effect <- 0.5

    control <- rnorm(n_patients / 2, 0, 1)
    treatment <- rnorm(n_patients / 2, treatment_effect, 1)

    reject <- t.test(treatment, control)$p.value < alpha

    c(power = as.numeric(reject), type1 = 0, EN = n_patients, ET = 24)
  }

  # Run in parallel chunks
  chunk_size <- ceiling(n_rep / n_cores)

  chunk_results <- parallel::mclapply(1:n_cores, function(core_id) {
    # Each core runs welford_mean_var independently
    welford_mean_var(
      sample_fn = simulate_one_trial,
      n_samples = chunk_size,
      theta = theta
    )
  }, mc.cores = n_cores)

  # Pool results from all cores
  pooled <- pool_welford_results(chunk_results)

  # Return metrics
  metrics <- pooled$mean
  attr(metrics, "variance") <- pooled$variance
  attr(metrics, "n_rep") <- pooled$n

  return(metrics)
}


# =============================================================================
# Example 3: Manual Welford implementation (for maximum control)
# =============================================================================

#' Adaptive trial simulator with manual Welford's algorithm
#'
#' For users who want complete control, implement Welford's algorithm directly.
manual_welford_simulator <- function(theta,
                                     fidelity = c("low", "med", "high"),
                                     seed = NULL,
                                     ...) {
  fidelity <- match.arg(fidelity)
  n_rep <- switch(fidelity, low = 200, med = 1000, high = 10000)

  if (!is.null(seed)) set.seed(seed)

  # Initialize Welford accumulators
  mean_vec <- c(power = 0, type1 = 0, EN = 0, ET = 0)
  M2_vec <- c(power = 0, type1 = 0, EN = 0, ET = 0)

  # Run simulations with online variance update
  for (i in 1:n_rep) {
    # Simulate one trial
    n_patients <- rpois(1, 10 * 24)
    treatment_effect <- 0.5
    alpha <- theta$alpha

    control <- rnorm(n_patients / 2, 0, 1)
    treatment <- rnorm(n_patients / 2, treatment_effect, 1)
    reject <- t.test(treatment, control)$p.value < alpha

    x <- c(power = as.numeric(reject), type1 = 0, EN = n_patients, ET = 24)

    # Welford's online update
    delta <- x - mean_vec
    mean_vec <- mean_vec + delta / i
    delta2 <- x - mean_vec
    M2_vec <- M2_vec + delta * delta2
  }

  # Compute variance of the mean
  variance <- M2_vec / (n_rep * (n_rep - 1))

  # Return metrics
  attr(mean_vec, "variance") <- variance
  attr(mean_vec, "n_rep") <- n_rep

  return(mean_vec)
}


# =============================================================================
# Usage Example
# =============================================================================

if (interactive()) {
  # Define design space
  bounds <- list(
    threshold = c(1.5, 3.0),
    alpha = c(0.01, 0.05)
  )

  # Define objective and constraints
  objective <- "EN"  # Minimize expected sample size
  constraints <- list(
    power = c("ge", 0.8),  # Power ≥ 0.8
    type1 = c("le", 0.05)   # Type I error ≤ 0.05
  )

  # Run Bayesian optimization with variance-aware simulator
  fit <- bo_calibrate(
    sim_fun = simple_adaptive_simulator,  # Uses Welford's algorithm
    bounds = bounds,
    objective = objective,
    constraints = constraints,
    n_init = 10,
    q = 2,
    budget = 50,
    seed = 2025
  )

  # Best design
  print(fit$best_theta)

  # Check variance was properly used
  print(head(fit$history$variance))
}


# =============================================================================
# Performance Comparison
# =============================================================================

#' Compare simulators with and without variance
#'
#' This demonstrates the performance difference between providing variance
#' vs using the nugget fallback.
if (FALSE) {  # Set to TRUE to run comparison
  library(tictoc)

  # Simulator WITHOUT variance (uses nugget)
  simulator_no_variance <- function(theta, fidelity = "high", ...) {
    n_rep <- 1000
    results <- replicate(n_rep, {
      n <- rpois(1, 240)
      reject <- runif(1) < 0.8
      c(power = reject, EN = n)
    })

    # Return only means (no variance attribute)
    rowMeans(results)
  }

  # Simulator WITH variance (uses Welford)
  simulator_with_variance <- function(theta, fidelity = "high", ...) {
    result <- welford_mean_var(
      sample_fn = function(i, theta) {
        n <- rpois(1, 240)
        reject <- runif(1) < 0.8
        c(power = reject, EN = n)
      },
      n_samples = 1000,
      theta = theta
    )

    metrics <- result$mean
    attr(metrics, "variance") <- result$variance
    return(metrics)
  }

  bounds <- list(x = c(0, 1))
  constraints <- list(power = c("ge", 0.8))

  # Run both
  tic("Without variance (nugget)")
  fit_no_var <- bo_calibrate(
    simulator_no_variance, bounds, "EN", constraints,
    n_init = 10, budget = 50, progress = FALSE
  )
  toc()

  tic("With variance (Welford)")
  fit_with_var <- bo_calibrate(
    simulator_with_variance, bounds, "EN", constraints,
    n_init = 10, budget = 50, progress = FALSE
  )
  toc()

  # Compare convergence
  cat("\nConvergence comparison:\n")
  cat("Without variance: ", nrow(fit_no_var$history), "evaluations\n")
  cat("With variance: ", nrow(fit_with_var$history), "evaluations\n")
  cat("Improvement: ",
      round((1 - nrow(fit_with_var$history) / nrow(fit_no_var$history)) * 100),
      "%\n")
}
