# Standard Test Problems for Benchmarking BATON
#
# This file defines three test problems for validating performance improvements:
# 1. toy_2d: Simple 2D problem for quick validation
# 2. high_dim_5d: 5D problem to test scalability
# 3. tight_constraints: 2D problem with difficult constraints
#
# All problems use synthetic simulators that mimic clinical trial calibration
# characteristics (power, type I error, sample size metrics).

library(BATON)

#' Simple 2D Test Problem
#'
#' Quadratic objective with feasible region.
#' Fast evaluation for quick benchmarking.
#'
#' @param theta named list with x1 and x2
#' @param fidelity simulation fidelity level
#' @param seed RNG seed
#' @return named vector of metrics with variance attributes
toy_2d_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  n_rep <- switch(fidelity,
    low = 200,
    med = 1000,
    high = 10000
  )

  # Use Welford's algorithm for variance estimation
  result <- welford_mean_var(
    sample_fn = function(i, theta) {
      x1 <- theta$x1
      x2 <- theta$x2

      # Objective: minimize distance from (0.6, 0.4)
      obj <- (x1 - 0.6)^2 + (x2 - 0.4)^2

      # Constraints (with noise)
      power <- plogis((x1 - 0.3) / 0.1 + rnorm(1, 0, 0.15))
      type1 <- plogis((x2 - 0.5) / 0.1 + rnorm(1, 0, 0.15))

      # Sample size (deterministic component + noise)
      EN <- 100 + obj * 200 + rnorm(1, 0, 10)

      c(power = power, type1 = type1, EN = EN)
    },
    n_samples = n_rep,
    theta = theta
  )

  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- n_rep

  return(metrics)
}

#' 2D Test Problem Configuration
#'
#' @export
toy_2d_problem <- function() {
  list(
    name = "toy_2d",
    sim_fun = toy_2d_sim,
    bounds = list(
      x1 = c(0, 1),
      x2 = c(0, 1)
    ),
    objective = "EN",
    constraints = list(
      power = c("ge", 0.75),
      type1 = c("le", 0.6)
    ),
    optimum = list(x1 = 0.6, x2 = 0.4),  # Approximate
    description = "Simple 2D quadratic with soft constraints"
  )
}


#' High-Dimensional 5D Test Problem
#'
#' Tests scalability of algorithms.
#' Multiplicative objective with constraint dependencies.
#'
#' @param theta named list with x1, x2, x3, x4, x5
#' @param fidelity simulation fidelity level
#' @param seed RNG seed
#' @return named vector of metrics with variance attributes
high_dim_5d_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  n_rep <- switch(fidelity,
    low = 200,
    med = 1000,
    high = 10000
  )

  result <- welford_mean_var(
    sample_fn = function(i, theta) {
      x <- unlist(theta)  # Extract as vector

      # Objective: Rosenbrock-like function
      obj <- 0
      for (j in 1:(length(x) - 1)) {
        obj <- obj + 100 * (x[j+1] - x[j]^2)^2 + (1 - x[j])^2
      }
      obj <- obj / 1000  # Scale down

      # Constraints (depend on multiple parameters)
      power <- plogis((mean(x[1:3]) - 0.5) / 0.1 + rnorm(1, 0, 0.12))
      type1 <- plogis((mean(x[3:5]) - 0.4) / 0.1 + rnorm(1, 0, 0.12))

      # Sample size increases with dimension
      EN <- 100 + obj * 500 + sum((x - 0.5)^2) * 100 + rnorm(1, 0, 15)

      c(power = power, type1 = type1, EN = EN)
    },
    n_samples = n_rep,
    theta = theta
  )

  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- n_rep

  return(metrics)
}

#' 5D Test Problem Configuration
#'
#' @export
high_dim_5d_problem <- function() {
  list(
    name = "high_dim_5d",
    sim_fun = high_dim_5d_sim,
    bounds = list(
      x1 = c(0, 1),
      x2 = c(0, 1),
      x3 = c(0, 1),
      x4 = c(0, 1),
      x5 = c(0, 1)
    ),
    objective = "EN",
    constraints = list(
      power = c("ge", 0.7),
      type1 = c("le", 0.5)
    ),
    optimum = NULL,  # Unknown
    description = "5D Rosenbrock-like function with coupled constraints"
  )
}


#' Tight Constraints 2D Test Problem
#'
#' Tests constraint handling when feasible region is small.
#' Initial random samples likely infeasible.
#'
#' @param theta named list with threshold and alpha
#' @param fidelity simulation fidelity level
#' @param seed RNG seed
#' @return named vector of metrics with variance attributes
tight_constraints_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  n_rep <- switch(fidelity,
    low = 200,
    med = 1000,
    high = 10000
  )

  result <- welford_mean_var(
    sample_fn = function(i, theta) {
      threshold <- theta$threshold
      alpha <- theta$alpha

      # Objective: minimize sample size
      # Optimal around (2.5, 0.025)
      EN <- 100 + (threshold - 2.5)^2 * 100 + (alpha - 0.025)^2 * 10000 + rnorm(1, 0, 5)

      # Tight constraints (small feasible region)
      power <- plogis((threshold - 2.0) / 0.15 + rnorm(1, 0, 0.1))
      type1 <- plogis((alpha - 0.03) / 0.005 + rnorm(1, 0, 0.1))

      c(power = power, type1 = type1, EN = EN)
    },
    n_samples = n_rep,
    theta = theta
  )

  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- n_rep

  return(metrics)
}

#' Tight Constraints 2D Test Problem Configuration
#'
#' @export
tight_constraints_problem <- function() {
  list(
    name = "tight_constraints",
    sim_fun = tight_constraints_sim,
    bounds = list(
      threshold = c(1.5, 3.5),
      alpha = c(0.01, 0.05)
    ),
    objective = "EN",
    constraints = list(
      power = c("ge", 0.85),  # High power requirement
      type1 = c("le", 0.05)   # Low type I error
    ),
    optimum = list(threshold = 2.5, alpha = 0.025),  # Approximate
    description = "2D problem with tight constraints (small feasible region)"
  )
}


#' Get All Test Problems
#'
#' Returns list of all available test problems.
#'
#' @return named list of problem configurations
#' @export
get_test_problems <- function() {
  list(
    toy_2d = toy_2d_problem(),
    high_dim_5d = high_dim_5d_problem(),
    tight_constraints = tight_constraints_problem()
  )
}


#' Evaluate True Optimum
#'
#' For problems with known optima, evaluate the simulator at the optimum
#' to get baseline performance.
#'
#' @param problem problem configuration from get_test_problems()
#' @param n_rep number of replications for evaluation
#' @return named vector of metrics at optimum
#' @export
evaluate_optimum <- function(problem, n_rep = 10000) {
  if (is.null(problem$optimum)) {
    warning("Problem '", problem$name, "' has no known optimum")
    return(NULL)
  }

  result <- problem$sim_fun(
    theta = problem$optimum,
    fidelity = "high",
    seed = 12345
  )

  return(result)
}


#' Print Problem Summary
#'
#' @param problem problem configuration
#' @export
print_problem_summary <- function(problem) {
  cat("\n=== Problem:", problem$name, "===\n")
  cat("Description:", problem$description, "\n")
  cat("Dimension:", length(problem$bounds), "\n")
  cat("Objective:", problem$objective, "(minimize)\n")
  cat("Constraints:\n")
  for (name in names(problem$constraints)) {
    constr <- problem$constraints[[name]]
    dir_str <- if (constr[1] == "ge") ">=" else "<="
    cat(sprintf("  %s %s %.2f\n", name, dir_str, as.numeric(constr[2])))
  }

  if (!is.null(problem$optimum)) {
    cat("Known optimum:\n")
    for (name in names(problem$optimum)) {
      cat(sprintf("  %s = %.3f\n", name, problem$optimum[[name]]))
    }

    cat("\nOptimum performance:\n")
    opt_metrics <- evaluate_optimum(problem)
    for (name in names(opt_metrics)) {
      cat(sprintf("  %s = %.3f\n", name, opt_metrics[[name]]))
    }
  }
  cat("\n")
}


# Example usage:
# problems <- get_test_problems()
# print_problem_summary(problems$toy_2d)
# print_problem_summary(problems$high_dim_5d)
# print_problem_summary(problems$tight_constraints)
