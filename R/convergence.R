# convergence.R
# Convergence detection utilities for Bayesian Optimization
#
# Functions for detecting optimization convergence within stages,
# enabling early stopping to save budget.

#' Check Optimization Convergence
#'
#' Detects if optimization has converged based on recent history.
#' Useful for early stopping in multi-stage BO.
#'
#' Checks:
#' - Improvement stagnation: best objective not improving for `patience` steps
#' - Acquisition flatline: acquisition scores below threshold (optional)
#'
#' @param history Data frame: optimization history with objective values
#' @param objective Character: name of objective column (default: "objective")
#' @param patience Integer: number of recent iterations to check for improvement
#'   (default: 3)
#' @param improvement_threshold Numeric: relative improvement threshold - smaller
#'   improvements are considered stagnation (default: 1e-3)
#' @param acq_flatline_threshold Numeric or NULL: if acquisition scores (column
#'   "acq_score") fall below this for `patience` iterations, consider converged.
#'   NULL disables this check (default: NULL)
#'
#' @return List with:
#'   \describe{
#'     \item{converged}{Logical: whether convergence detected}
#'     \item{reason}{Character: reason for convergence, or NULL if not converged}
#'     \item{best_objective}{Numeric: best feasible objective value}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After each BO iteration
#' conv <- check_convergence(
#'   history = fit$history,
#'   objective = "EN",
#'   patience = 5,
#'   improvement_threshold = 0.001
#' )
#'
#' if (conv$converged) {
#'   cat("Converged:", conv$reason, "\n")
#'   cat("Best objective:", conv$best_objective, "\n")
#' }
#' }
check_convergence <- function(
  history,
  objective = "objective",
  patience = 3,
  improvement_threshold = 1e-3,
  acq_flatline_threshold = NULL
) {

  if (is.null(history) || nrow(history) == 0) {
    return(list(converged = FALSE, reason = NULL, best_objective = NA_real_))
  }

  # Check for feasible column
  if (!"feasible" %in% names(history)) {
    # Assume all feasible if column missing
    feasible <- history
  } else {
    feasible <- history[history$feasible == TRUE, , drop = FALSE]
  }

  if (nrow(feasible) == 0) {
    return(list(converged = FALSE, reason = NULL, best_objective = NA_real_))
  }

  # Get objective column
  obj_col <- if (objective %in% names(feasible)) objective else "objective"
  if (!obj_col %in% names(feasible)) {
    return(list(converged = FALSE, reason = NULL, best_objective = NA_real_))
  }

  obj_values <- feasible[[obj_col]]
  obj_values <- obj_values[!is.na(obj_values)]

  # Guard: all objective values are NA

  if (length(obj_values) == 0) {
    return(list(converged = FALSE, reason = NULL, best_objective = NA_real_))
  }

  best_objective <- min(obj_values)

  if (length(obj_values) < patience + 1) {
    return(list(converged = FALSE, reason = NULL, best_objective = best_objective))
  }

  # Check improvement stagnation
  n <- length(obj_values)
  cummin_obj <- cummin(obj_values)

  # Improvement in last 'patience' steps
  recent_best <- cummin_obj[n]
  earlier_best <- cummin_obj[n - patience]

  # Compute relative improvement (works for any objective including zero/negative)
  # Only skip if earlier_best is exactly zero (would cause division by zero)
  if (earlier_best != 0) {
    rel_improvement <- (earlier_best - recent_best) / abs(earlier_best)

    if (rel_improvement < improvement_threshold) {
      return(list(
        converged = TRUE,
        reason = sprintf("improvement stagnation (%.2e < %.2e)",
                         rel_improvement, improvement_threshold),
        best_objective = best_objective
      ))
    }
  } else if (recent_best == earlier_best) {
    # Both are zero - stagnation
    return(list(
      converged = TRUE,
      reason = sprintf("improvement stagnation (no change from zero)"),
      best_objective = best_objective
    ))
  }

  # Check acquisition flatline (if available and threshold set)
  if (!is.null(acq_flatline_threshold) && "acq_score" %in% names(history)) {
    recent_acq <- utils::tail(history$acq_score, patience)
    recent_acq <- recent_acq[!is.na(recent_acq)]

    if (length(recent_acq) >= patience && all(recent_acq < acq_flatline_threshold)) {
      return(list(
        converged = TRUE,
        reason = sprintf("acquisition flatline (max=%.2e < %.2e)",
                         max(recent_acq), acq_flatline_threshold),
        best_objective = best_objective
      ))
    }
  }

  list(converged = FALSE, reason = NULL, best_objective = best_objective)
}


#' Early Stopping Callback for BO
#'
#' Creates a callback function for use with bo_calibrate that checks
#' convergence after each iteration.
#'
#' @param patience Integer: number of iterations to check (default: 5)
#' @param improvement_threshold Numeric: minimum relative improvement (default: 1e-3)
#' @param verbose Logical: print message when stopping early (default: TRUE)
#'
#' @return Function suitable for use as a callback in bo_calibrate
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # For future integration with bo_calibrate callback system
#' callback <- early_stopping_callback(patience = 5)
#' fit <- bo_calibrate(..., callbacks = list(callback))
#' }
early_stopping_callback <- function(
  patience = 5,
  improvement_threshold = 1e-3,
  verbose = TRUE
) {

  function(history, iter, ...) {
    if (is.null(history) || nrow(history) < patience + 1) {
      return(list(stop = FALSE))
    }

    result <- check_convergence(
      history = history,
      patience = patience,
      improvement_threshold = improvement_threshold
    )

    if (result$converged && verbose) {
      cat(sprintf("\n[Early stopping] %s at iteration %d\n",
                  result$reason, iter))
      cat(sprintf("  Best objective: %.4f\n", result$best_objective))
    }

    list(stop = result$converged, reason = result$reason)
  }
}
