# bounds.R
# Bound manipulation utilities for Bayesian Optimization
#
# Functions for narrowing, refining, and validating parameter bounds
# during multi-stage optimization.

#' Conservative Quantile-Based Bound Narrowing
#'
#' Narrows parameter bounds based on observed feasible designs, with safeguards:
#' - Uses quantiles (not min/max) to avoid over-narrowing from outliers
#' - Enforces minimum width relative to original bounds
#' - Requires minimum number of feasible designs before narrowing
#' - Softens narrowing when feasible rate is low
#' - Boundary-aware expansion when clamping reduces width
#'
#' @param history Data frame from previous BO stage with columns including
#'   parameter values and `feasible` (logical)
#' @param bounds Current bounds (named list, each element: c(lower, upper))
#' @param original_bounds Original wide bounds for minimum width calculation
#'   (default: same as bounds)
#' @param quantile_range Numeric vector of length 2: quantiles to use for
#'   narrowing (default: c(0.05, 0.95))
#' @param min_width_frac Numeric: minimum width as fraction of original bounds
#'   (default: 0.30 = at least 30\% of original range)
#' @param min_feasible Integer: minimum feasible designs required before
#'   narrowing is applied (default: 10)
#' @param soften_threshold Numeric: if feasible rate < this threshold, soften
#'   the narrowing (default: 0.10 = 10\%)
#' @param soften_factor Numeric: blend factor when softening - 0 = no narrowing,
#'   1 = full narrowing (default: 0.5)
#' @param verbose Logical: print narrowing details?
#'
#' @return List with:
#'   \describe{
#'     \item{bounds}{Narrowed bounds (named list)}
#'     \item{narrowed}{Logical: whether any narrowing was applied}
#'     \item{diagnostics}{List with narrowing decision details:
#'       n_total, n_feasible, feasible_rate, softened, params_narrowed,
#'       params_skipped, reason_skipped}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After Stage 1 BO
#' history_s1 <- fit_s1$history
#'
#' # Narrow bounds for Stage 2
#' result <- narrow_bounds_conservative(
#'   history = history_s1,
#'   bounds = bounds_wide,
#'   original_bounds = bounds_wide,
#'   min_width_frac = 0.30,
#'   min_feasible = 10
#' )
#'
#' if (result$narrowed) {
#'   bounds_s2 <- result$bounds
#' }
#' }
narrow_bounds_conservative <- function(
  history,
  bounds,
  original_bounds = bounds,
  quantile_range = c(0.05, 0.95),
  min_width_frac = 0.30,
  min_feasible = 10,
  soften_threshold = 0.10,
  soften_factor = 0.5,
  verbose = TRUE
) {

  # Guard: NULL or empty history
  if (is.null(history) || nrow(history) == 0) {
    return(list(
      bounds = bounds,
      narrowed = FALSE,
      diagnostics = list(
        n_total = 0,
        n_feasible = 0,
        feasible_rate = NA_real_,
        softened = FALSE,
        params_narrowed = character(),
        params_skipped = character(),
        reason_skipped = "empty_history"
      )
    ))
  }

  diagnostics <- list(
    n_total = nrow(history),
    n_feasible = sum(history$feasible, na.rm = TRUE),
    feasible_rate = NA_real_,
    softened = FALSE,
    params_narrowed = character(),
    params_skipped = character()
  )

  # Check feasible designs
  feasible <- history[history$feasible == TRUE, , drop = FALSE]
  n_feasible <- nrow(feasible)
  n_total <- nrow(history)
  feasible_rate <- n_feasible / n_total

  diagnostics$feasible_rate <- feasible_rate

  # Guard: not enough feasible designs
  if (n_feasible < min_feasible) {
    if (verbose) {
      cat(sprintf("  Skipping bound narrowing: only %d feasible designs (need %d)\n",
                  n_feasible, min_feasible))
    }
    diagnostics$reason_skipped <- "insufficient_feasible"
    return(list(bounds = bounds, narrowed = FALSE, diagnostics = diagnostics))
  }

  # Determine if we should soften narrowing
  soften <- feasible_rate < soften_threshold
  if (soften) {
    diagnostics$softened <- TRUE
    if (verbose) {
      cat(sprintf("  Softening narrowing: feasible rate %.1f%% < %.1f%% threshold\n",
                  100 * feasible_rate, 100 * soften_threshold))
    }
  }

  # Narrow each parameter
  new_bounds <- bounds
  any_narrowed <- FALSE

  for (param in names(bounds)) {
    if (!param %in% names(feasible)) {
      diagnostics$params_skipped <- c(diagnostics$params_skipped, param)
      next
    }

    values <- feasible[[param]]
    if (all(is.na(values))) {
      diagnostics$params_skipped <- c(diagnostics$params_skipped, param)
      next
    }

    # Current and original widths
    current_lower <- bounds[[param]][1]
    current_upper <- bounds[[param]][2]
    current_width <- current_upper - current_lower

    orig_lower <- original_bounds[[param]][1]
    orig_upper <- original_bounds[[param]][2]
    orig_width <- orig_upper - orig_lower

    # Compute quantile-based range
    q_range <- stats::quantile(values, probs = quantile_range, na.rm = TRUE)
    q_lower <- q_range[1]
    q_upper <- q_range[2]
    q_width <- q_upper - q_lower

    # Apply softening if needed (narrow less aggressively)
    if (soften) {
      # Blend towards current bounds
      q_lower <- soften_factor * q_lower + (1 - soften_factor) * current_lower
      q_upper <- soften_factor * q_upper + (1 - soften_factor) * current_upper
      q_width <- q_upper - q_lower
    }

    # Enforce minimum width
    min_width <- min_width_frac * orig_width
    if (q_width < min_width) {
      # Expand symmetrically to meet minimum width
      center <- (q_lower + q_upper) / 2
      q_lower <- center - min_width / 2
      q_upper <- center + min_width / 2
    }

    # Clamp to original bounds (never expand beyond original)
    q_lower <- max(q_lower, orig_lower)
    q_upper <- min(q_upper, orig_upper)

    # If clamping reduced width below minimum, expand the other direction
    actual_width <- q_upper - q_lower
    if (actual_width < min_width) {
      shortfall <- min_width - actual_width
      # Try to expand on the side that wasn't clamped
      if (q_lower == orig_lower && q_upper < orig_upper) {
        # Lower was clamped, expand upper
        q_upper <- min(q_upper + shortfall, orig_upper)
      } else if (q_upper == orig_upper && q_lower > orig_lower) {
        # Upper was clamped, expand lower
        q_lower <- max(q_lower - shortfall, orig_lower)
      }
      # Note: if both were clamped, we accept the narrower width
    }

    # Only update if actually narrower
    if (q_lower > current_lower || q_upper < current_upper) {
      new_bounds[[param]] <- c(q_lower, q_upper)
      any_narrowed <- TRUE
      diagnostics$params_narrowed <- c(diagnostics$params_narrowed, param)

      if (verbose) {
        cat(sprintf("  %s: [%.3f, %.3f] -> [%.3f, %.3f] (%.0f%% of original)\n",
                    param, current_lower, current_upper, q_lower, q_upper,
                    100 * (q_upper - q_lower) / orig_width))
      }
    }
  }

  if (verbose && !any_narrowed) {
    cat("  No parameters narrowed (all within acceptable range)\n")
  }

  list(bounds = new_bounds, narrowed = any_narrowed, diagnostics = diagnostics)
}


#' Validate Parameter Bounds
#'
#' Checks that bounds are properly structured with valid lower < upper values.
#'
#' @param bounds Named list of bounds (each element: c(lower, upper))
#' @param param_names Optional character vector of expected parameter names
#'
#' @return Logical: TRUE if valid, otherwise throws error
#' @keywords internal
validate_bounds <- function(bounds, param_names = NULL) {
  if (!is.list(bounds)) {
    stop("`bounds` must be a list")
  }

  if (is.null(names(bounds)) || any(names(bounds) == "")) {
    stop("`bounds` must be a named list")
  }

  for (param in names(bounds)) {
    b <- bounds[[param]]
    if (length(b) != 2) {
      stop(sprintf("Bounds for '%s' must have length 2, got %d", param, length(b)))
    }
    if (!is.numeric(b)) {
      stop(sprintf("Bounds for '%s' must be numeric", param))
    }
    if (b[1] >= b[2]) {
      stop(sprintf("Invalid bounds for '%s': lower (%.4f) >= upper (%.4f)",
                   param, b[1], b[2]))
    }
  }

  if (!is.null(param_names)) {
    missing <- setdiff(param_names, names(bounds))
    if (length(missing) > 0) {
      stop(sprintf("Missing bounds for parameters: %s", paste(missing, collapse = ", ")))
    }
  }

  TRUE
}
