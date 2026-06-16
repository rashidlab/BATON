# init_stopping.R
# GP-Based Initialization Stopping for Bayesian Optimization
#
# Monitors GP model quality during initialization to determine when
# sufficient points have been collected to begin the BO loop.
#
# Key insight: If adding more initialization points doesn't significantly
# change the GP model, we've captured the landscape sufficiently.

#' Check Initialization Sufficiency via GP Diagnostics
#'
#' Fits a GP every `k` points during initialization and checks multiple
#' convergence metrics to determine if the design space is sufficiently
#' explored. This can save budget by stopping initialization early when
#' the GP model has stabilized.
#'
#' @param X Design matrix (n x d) of evaluated points
#' @param y Response vector (n) of objective values
#' @param checkpoints List of previous checkpoint results (for comparison)
#' @param metrics Character vector of metrics to track. Options:
#'   \describe{
#'     \item{"variance"}{Max/mean prediction variance on test grid (recommended)}
#'     \item{"loo"}{Leave-one-out cross-validation RMSE}
#'     \item{"hyperparams"}{GP hyperparameter stability}
#'     \item{"loglik"}{Log marginal likelihood}
#'   }
#' @param test_grid Optional test grid for variance computation. If NULL,
#'   a Sobol sequence is generated.
#' @param n_test Number of test points for variance computation (default: 200)
#' @param threshold Relative change threshold for convergence (default: 0.10)
#' @param window Number of consecutive stable checkpoints required (default: 2
#' @param verbose Logical: print diagnostic messages (default: FALSE)
#'
#' @return List with:
#'   \describe{
#'     \item{metrics}{Named list of computed metric values}
#'     \item{converged}{Logical: whether convergence detected}
#'     \item{reason}{Character: reason for convergence, or NULL}
#'     \item{recommendation}{Character: suggested action}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # During initialization loop
#' checkpoints <- list()
#' for (i in seq(20, n_init, by = 20)) {
#'   result <- check_init_sufficiency(
#'     X = design[1:i, ],
#'     y = responses[1:i],
#'     checkpoints = checkpoints,
#'     metrics = c("variance", "loo"),
#'     threshold = 0.10
#'   )
#'   checkpoints <- c(checkpoints, list(result$metrics))
#'
#'   if (result$converged) {
#'     message("Initialization sufficient at ", i, " points: ", result$reason)
#'     break
#'   }
#' }
#' }
check_init_sufficiency <- function(
    X,
    y,
    checkpoints = list(),
    metrics = c("variance"),
    test_grid = NULL,
    n_test = 200,
    threshold = 0.10,
    window = 2,
    verbose = FALSE
) {

  # Input validation
  if (!is.matrix(X)) X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  d <- ncol(X)

  if (n < 5) {
    return(list(
      metrics = NULL,
      converged = FALSE,
      reason = NULL,
      recommendation = "Need at least 5 points"
    ))
  }

  # Fit GP model
  gp <- tryCatch({
    DiceKriging::km(
      design = X,
      response = y,
      covtype = "matern5_2",
      control = list(trace = FALSE, pop.size = 50, max.generations = 20),
      optim.method = "gen"
    )
  }, error = function(e) {
    if (verbose) message("GP fitting failed: ", e$message)
    NULL
  })

  if (is.null(gp)) {
    return(list(
      metrics = NULL,
      converged = FALSE,
      reason = NULL,
      recommendation = "GP fitting failed"
    ))
  }

  # Initialize results
  current_metrics <- list(n_points = n)

  # Metric 1: Prediction variance on test grid (recommended)
  if ("variance" %in% metrics) {
    if (is.null(test_grid)) {
      # Generate Sobol test grid
      test_grid <- generate_test_grid(X, n_test)
    }

    pred <- tryCatch({
      predict(gp, newdata = as.data.frame(test_grid), type = "UK")
    }, error = function(e) NULL)

    if (!is.null(pred)) {
      current_metrics$max_var <- max(pred$sd^2)
      current_metrics$mean_var <- mean(pred$sd^2)
      current_metrics$median_var <- median(pred$sd^2)
      # Coefficient of variation of variances (heterogeneity)
      current_metrics$cv_var <- sd(pred$sd^2) / mean(pred$sd^2)
    }
  }

  # Metric 2: Leave-one-out cross-validation
  if ("loo" %in% metrics) {
    loo <- tryCatch({
      DiceKriging::leaveOneOut.km(gp, type = "UK")
    }, error = function(e) NULL)

    if (!is.null(loo)) {
      current_metrics$loo_rmse <- sqrt(mean(loo$mean^2))
      current_metrics$loo_mae <- mean(abs(loo$mean))
      # Standardized LOO (normalized by response range)
      y_range <- diff(range(y))
      if (y_range > 0) {
        current_metrics$loo_nrmse <- current_metrics$loo_rmse / y_range
      }
    }
  }

  # Metric 3: Hyperparameter values
  if ("hyperparams" %in% metrics) {
    current_metrics$lengthscales <- gp@covariance@range.val
    current_metrics$process_var <- gp@covariance@sd2
    current_metrics$nugget <- gp@covariance@nugget
  }

  # Metric 4: Log marginal likelihood
  if ("loglik" %in% metrics) {
    current_metrics$loglik <- gp@logLik
    # Normalized by n (per-point loglik)
    current_metrics$loglik_per_point <- gp@logLik / n
  }

  # Check convergence against previous checkpoints
  converged <- FALSE
  reason <- NULL

  if (length(checkpoints) >= window) {
    conv_result <- check_metric_convergence(
      current = current_metrics,
      previous = checkpoints,
      window = window,
      threshold = threshold,
      verbose = verbose
    )
    converged <- conv_result$converged
    reason <- conv_result$reason
  }

  # Generate recommendation
  recommendation <- if (converged) {
    sprintf("Initialization sufficient at %d points. %s", n, reason)
  } else if (length(checkpoints) < window) {
    sprintf("Need %d more checkpoints for comparison", window - length(checkpoints))
  } else {
    "Continue initialization"
  }

  list(
    metrics = current_metrics,
    converged = converged,
    reason = reason,
    recommendation = recommendation
  )
}


#' Check if Metrics Have Converged
#'
#' Compares current metrics to previous checkpoints to detect stability.
#'
#' @param current Named list of current metric values
#' @param previous List of previous checkpoint metric lists
#' @param window Number of checkpoints to compare
#' @param threshold Relative change threshold
#' @param verbose Print diagnostic info
#'
#' @return List with converged (logical) and reason (character)
#' @keywords internal
check_metric_convergence <- function(current, previous, window, threshold, verbose) {

  # Get the last 'window' checkpoints
  recent <- utils::tail(previous, window)

  # Metrics to check (in order of priority)
  check_metrics <- c("max_var", "mean_var", "loo_rmse", "loglik_per_point")

  for (metric in check_metrics) {
    if (!is.null(current[[metric]]) &&
        all(sapply(recent, function(x) !is.null(x[[metric]])))) {

      current_val <- current[[metric]]
      prev_vals <- sapply(recent, `[[`, metric)

      # Compute relative changes
      rel_changes <- abs(current_val - prev_vals) / pmax(abs(prev_vals), 1e-10)

      if (all(rel_changes < threshold)) {
        reason <- sprintf("%s stable (changes: %s < %.1f%%)",
                          metric,
                          paste(sprintf("%.1f%%", rel_changes * 100), collapse = ", "),
                          threshold * 100)

        if (verbose) {
          message(sprintf("[Init stopping] %s", reason))
        }

        return(list(converged = TRUE, reason = reason))
      }
    }
  }

  list(converged = FALSE, reason = NULL)
}


#' Generate Test Grid for Variance Computation
#'
#' Creates a space-filling test grid using Sobol sequences.
#'
#' @param X Design matrix to extract bounds from
#' @param n Number of test points
#'
#' @return Matrix of test points (n x d)
#' @keywords internal
generate_test_grid <- function(X, n = 200) {
  d <- ncol(X)
  lower <- apply(X, 2, min)
  upper <- apply(X, 2, max)

  # Slightly expand bounds to check edges
  range_vec <- upper - lower
  lower <- lower - 0.05 * range_vec
  upper <- upper + 0.05 * range_vec

  # Generate Sobol sequence (or fallback to random)
  if (requireNamespace("randtoolbox", quietly = TRUE)) {
    sobol_01 <- randtoolbox::sobol(n, dim = d)
  } else {
    # Fallback: stratified random
    sobol_01 <- matrix(runif(n * d), nrow = n, ncol = d)
  }

  # Scale to bounds
  test_grid <- sweep(sobol_01, 2, upper - lower, "*")
  test_grid <- sweep(test_grid, 2, lower, "+")

  # Preserve column names from X (needed for DiceKriging::predict)
  test_grid <- as.matrix(test_grid)
  colnames(test_grid) <- colnames(X)

  test_grid
}


#' Create Initialization Stopping Configuration
#'
#' Helper function to create a configuration list for initialization stopping.
#'
#' @param enabled Logical: enable initialization stopping (default: TRUE)
#' @param min_init Minimum initialization points before checking (default: 20)
#' @param check_every Check every k points (default: 20)
#' @param metrics Character vector of metrics to use (default: "variance")
#' @param threshold Relative change threshold (default: 0.10)
#' @param window Consecutive stable checkpoints required (default: 2)
#' @param max_init Maximum initialization even if not converged (default: NULL, use n_init)
#' @param verbose Print messages (default: FALSE)
#'
#' @return List of configuration parameters
#' @export
#'
#' @examples
#' \dontrun{
#' # Default configuration
#' config <- init_stopping_config()
#'
#' # Aggressive early stopping
#' config <- init_stopping_config(
#'   min_init = 15,
#'   check_every = 10,
#'   threshold = 0.15
#' )
#'
#' # Conservative (high precision)
#' config <- init_stopping_config(
#'   min_init = 40,
#'   threshold = 0.05,
#'   window = 3
#' )
#' }
init_stopping_config <- function(
    enabled = TRUE,
    min_init = 20,
    check_every = 20,
    metrics = c("variance"),
    threshold = 0.10,
    window = 2,
    max_init = NULL,
    verbose = FALSE
) {
  list(
    enabled = enabled,
    min_init = min_init,
    check_every = check_every,
    metrics = metrics,
    threshold = threshold,
    window = window,
    max_init = max_init,
    verbose = verbose
  )
}


#' Run Initialization with Optional Early Stopping
#'
#' Wrapper for initialization phase that supports GP-based early stopping.
#'
#' @param sim_fun Simulation function
#' @param design Initial design matrix (from LHS or similar)
#' @param objective Name of the objective metric (character). If NULL, uses
#'   the first element of the simulator result.
#' @param fidelity Fidelity level for evaluations
#' @param init_config Configuration from init_stopping_config()
#' @param seed Random seed
#' @param verbose Print progress
#'
#' @return List with:
#'   \describe{
#'     \item{X}{Final design matrix}
#'     \item{y}{Objective values}
#'     \item{results}{Full results list}
#'     \item{n_used}{Number of initialization points used}
#'     \item{stopped_early}{Logical: whether early stopping triggered}
#'     \item{checkpoints}{List of checkpoint metrics}
#'   }
#'
#' @keywords internal
run_init_with_stopping <- function(
    sim_fun,
    design,
    objective = NULL,
    fidelity = "low",
    init_config = init_stopping_config(),
    seed = NULL,
    verbose = TRUE
) {

  n_total <- nrow(design)
  d <- ncol(design)

  # Storage
  results <- vector("list", n_total)
  y <- numeric(n_total)
  checkpoints <- list()
  stopped_early <- FALSE
  n_used <- 0

  for (i in seq_len(n_total)) {
    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed + i)

    # Evaluate point
    theta <- as.list(design[i, ])
    names(theta) <- colnames(design)

    result <- sim_fun(theta, fidelity = fidelity, seed = seed + i)
    results[[i]] <- result

    # Extract metrics from result - handle both vector and list-style outputs
    # Simulators can return:
    #   1. Named vector: c(EN = 0.1, power = 0.85)
    #   2. List with metrics element: list(metrics = c(EN = 0.1, ...), variance = ...)
    #   3. List with named elements: list(EN = 0.1, power = 0.85)
    if (is.list(result) && !is.null(result$metrics)) {
      # Case 2: list with metrics element
      metrics <- result$metrics
    } else {
      # Case 1 or 3: named vector or plain list
      metrics <- result
    }

    # Extract objective value - use named metric if specified, else first element
    if (!is.null(objective) && objective %in% names(metrics)) {
      y[i] <- as.numeric(metrics[[objective]])
    } else if (!is.null(objective)) {
      # Objective name specified but not found in result
      warning(sprintf(
        "Objective '%s' not found in simulator result at iteration %d. Using first element.",
        objective, i
      ))
      y[i] <- as.numeric(metrics[[1]])
    } else {
      # No objective specified - use first element (legacy behavior)
      y[i] <- as.numeric(metrics[[1]])
    }
    n_used <- i

    if (verbose && i %% 10 == 0) {
      message(sprintf("  Init %d/%d: obj = %.4f", i, n_total, y[i]))
    }

    # Check stopping criteria
    if (init_config$enabled &&
        i >= init_config$min_init &&
        i %% init_config$check_every == 0) {

      check <- check_init_sufficiency(
        X = design[1:i, , drop = FALSE],
        y = y[1:i],
        checkpoints = checkpoints,
        metrics = init_config$metrics,
        threshold = init_config$threshold,
        window = init_config$window,
        verbose = init_config$verbose
      )

      checkpoints <- c(checkpoints, list(check$metrics))

      if (check$converged) {
        if (verbose) {
          message(sprintf("[Init stopping] %s at %d/%d points",
                          check$reason, i, n_total))
        }
        stopped_early <- TRUE
        break
      }
    }

    # Respect max_init if set
    if (!is.null(init_config$max_init) && i >= init_config$max_init) {
      if (verbose) {
        message(sprintf("[Init stopping] Reached max_init = %d", init_config$max_init))
      }
      break
    }
  }

  list(
    X = design[1:n_used, , drop = FALSE],
    y = y[1:n_used],
    results = results[1:n_used],
    n_used = n_used,
    stopped_early = stopped_early,
    checkpoints = checkpoints,
    saved = n_total - n_used
  )
}


#' Visualize Initialization Convergence
#'
#' Creates diagnostic plots showing how metrics evolved during initialization.
#'
#' @param checkpoints List of checkpoint metrics from check_init_sufficiency()
#' @param threshold Convergence threshold (for reference line)
#'
#' @return ggplot object (if ggplot2 available) or NULL
#' @export
plot_init_convergence <- function(checkpoints, threshold = 0.10) {

  if (length(checkpoints) < 2) {
    message("Need at least 2 checkpoints to plot")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 required for plotting")
    return(invisible(NULL))
  }

  # Extract metrics into data frame
  plot_data <- do.call(rbind, lapply(seq_along(checkpoints), function(i) {
    cp <- checkpoints[[i]]
    data.frame(
      checkpoint = i,
      n_points = cp$n_points,
      max_var = cp$max_var %||% NA,
      mean_var = cp$mean_var %||% NA,
      loo_rmse = cp$loo_rmse %||% NA,
      stringsAsFactors = FALSE
    )
  }))

  # Compute relative changes
  plot_data$max_var_change <- c(NA, abs(diff(plot_data$max_var)) /
                                   head(plot_data$max_var, -1))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n_points)) +
    ggplot2::geom_line(ggplot2::aes(y = max_var), color = "blue", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = max_var), color = "blue", size = 3) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                        color = "red", alpha = 0.7) +
    ggplot2::labs(
      title = "Initialization Convergence: Max Prediction Variance",
      subtitle = sprintf("Threshold: %.0f%% relative change", threshold * 100),
      x = "Number of Initialization Points",
      y = "Max Prediction Variance"
    ) +
    ggplot2::theme_minimal()

  p
}


# Null-coalescing operator (if not already defined)
`%||%` <- function(x, y) if (is.null(x)) y else x
