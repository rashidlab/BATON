#' Normalise user-supplied constraint specification
#' @keywords internal
parse_constraints <- function(constraints) {
  if (is.null(constraints) || length(constraints) == 0) {
    return(tibble::tibble(metric = character(), direction = character(), threshold = numeric()))
  }
  if (!is.list(constraints)) {
    stop("`constraints` must be a named list.", call. = FALSE)
  }
  # Use base R vapply instead of purrr::map_chr/map_dbl to avoid "In index: X" errors
  constraint_names <- names(constraints)
  directions <- vapply(constraints, function(x) {
    dir <- x[[1]]
    if (!dir %in% c("ge", "le")) {
      stop("Constraint directions must be 'ge' or 'le'.", call. = FALSE)
    }
    dir
  }, FUN.VALUE = character(1))
  thresholds <- vapply(constraints, function(x) as.numeric(x[[2]]), FUN.VALUE = numeric(1))

  tibble::tibble(
    metric = constraint_names,
    direction = directions,
    threshold = thresholds
  )
}

#' Deterministic feasibility check for observed metrics
#' @keywords internal
is_feasible <- function(metrics, constraint_tbl) {
  if (nrow(constraint_tbl) == 0L) {
    return(TRUE)
  }
  # Use base R vapply instead of purrr::pmap_lgl to avoid "In index: X" errors
  results <- vapply(seq_len(nrow(constraint_tbl)), function(i) {
    metric <- constraint_tbl$metric[i]
    direction <- constraint_tbl$direction[i]
    threshold <- constraint_tbl$threshold[i]

    value <- metrics[[metric]]
    if (is.null(value) || is.na(value)) {
      return(FALSE)
    }
    if (direction == "ge") {
      value >= threshold
    } else {
      value <= threshold
    }
  }, FUN.VALUE = logical(1))
  all(results)
}

#' Probability of feasibility under Gaussian predictive distribution
#' @keywords internal
prob_feasibility <- function(mean, sd, constraint_tbl) {
  if (nrow(constraint_tbl) == 0L) {
    return(1)
  }
  # Use base R vapply instead of purrr::pmap_dbl to avoid "In index: X" errors
  probs <- vapply(seq_len(nrow(constraint_tbl)), function(i) {
    metric <- constraint_tbl$metric[i]
    direction <- constraint_tbl$direction[i]
    threshold <- constraint_tbl$threshold[i]

    mu <- mean[[metric]]
    sigma <- sd[[metric]]
    if (is.null(mu) || is.na(mu) || is.null(sigma) || is.na(sigma) || sigma <= 0) {
      return(0)
    }
    if (direction == "ge") {
      stats::pnorm((mu - threshold) / sigma, lower.tail = TRUE)
    } else {
      stats::pnorm((threshold - mu) / sigma, lower.tail = TRUE)
    }
  }, FUN.VALUE = numeric(1))
  prod(probs)
}


#' Batch probability of feasibility (vectorized across candidates)
#'
#' PERFORMANCE: Computes P(feasible) for all candidates at once instead of
#' calling prob_feasibility in a loop. This reduces purrr::map overhead from
#' O(n_candidates × n_constraints) to O(n_constraints).
#'
#' @param pred list of predictions from predict_surrogates (contains mean and sd for each metric)
#' @param constraint_tbl tibble from parse_constraints
#'
#' @return numeric vector of feasibility probabilities for each candidate
#' @keywords internal
prob_feasibility_batch <- function(pred, constraint_tbl) {
  n_candidates <- length(pred[[1]]$mean)

  if (nrow(constraint_tbl) == 0L) {
    return(rep(1, n_candidates))
  }

  # Initialize cumulative probability (product across constraints)
  prob_feas <- rep(1, n_candidates)

  # Loop over constraints (typically 2-4), vectorize across all candidates
  for (j in seq_len(nrow(constraint_tbl))) {
    metric <- constraint_tbl$metric[j]
    direction <- constraint_tbl$direction[j]
    threshold <- constraint_tbl$threshold[j]

    if (!metric %in% names(pred)) {
      # Missing metric: set probability to 0 for all candidates
      prob_feas <- rep(0, n_candidates)
      break
    }

    mu_vec <- pred[[metric]]$mean
    sd_vec <- pred[[metric]]$sd

    # Handle invalid sd values (NA, NULL, or <= 0)
    invalid <- is.na(sd_vec) | sd_vec <= 0
    sd_vec[invalid] <- 1e-10  # Small positive to avoid division by zero

    if (direction == "ge") {
      # P(metric >= threshold) = P(Z >= (threshold - mu) / sd) = Phi((mu - threshold) / sd)
      constraint_prob <- stats::pnorm((mu_vec - threshold) / sd_vec)
    } else {
      # P(metric <= threshold) = P(Z <= (threshold - mu) / sd) = Phi((threshold - mu) / sd)
      constraint_prob <- stats::pnorm((threshold - mu_vec) / sd_vec)
    }

    # Set invalid entries to 0
    constraint_prob[invalid] <- 0

    # Multiply into cumulative (assuming independence)
    prob_feas <- prob_feas * constraint_prob
  }

  prob_feas
}


#' Filter History for New Constraints
#'
#' When adding constraints in later stages, filters history to only include
#' rows that have the required metrics (not NA). This prevents poisoning
#' the GP surrogate with incomplete data from earlier stages that didn't
#' evaluate the new constraint metrics.
#'
#' @param history Data frame from previous stage(s)
#' @param new_constraints Named list of new constraints being added
#'   (e.g., list(power = c("ge", 0.80)))
#' @param verbose Logical: print filtering details?
#'
#' @return List with:
#'   \describe{
#'     \item{history}{Filtered history data frame, or NULL if all rows dropped}
#'     \item{n_original}{Integer: original row count}
#'     \item{n_retained}{Integer: retained row count}
#'     \item{metrics_checked}{Character: names of metrics that were checked}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Stage 2 adds PET constraint not evaluated in Stage 1
#' new_constraints_s2 <- list(PET_fut_null = c("ge", 0.30))
#'
#' result <- filter_history_for_constraints(
#'   history = fit_s1$history,
#'   new_constraints = new_constraints_s2
#' )
#'
#' if (result$n_retained > 0) {
#'   # Use filtered history for Stage 2 warmstart
#'   initial_history_s2 <- result$history
#' } else {
#'   # Stage 2 must start fresh (no usable history)
#'   initial_history_s2 <- NULL
#' }
#' }
filter_history_for_constraints <- function(history, new_constraints, verbose = TRUE) {

  if (is.null(history) || nrow(history) == 0) {
    return(list(history = NULL, n_original = 0, n_retained = 0, metrics_checked = character()))
  }

  n_original <- nrow(history)
  required_metrics <- names(new_constraints)

  if (length(required_metrics) == 0) {
    return(list(history = history, n_original = n_original,
                n_retained = n_original, metrics_checked = character()))
  }

  # Check which metrics are available
  metrics_present <- intersect(required_metrics, names(history))
  metrics_missing <- setdiff(required_metrics, names(history))

  if (length(metrics_missing) > 0 && verbose) {
    cat(sprintf("  Warning: Constraint metrics not in history columns: %s\n",
                paste(metrics_missing, collapse = ", ")))
    cat("  Will check metrics list-column if available\n")
  }

  # Build mask: TRUE if row has all required metrics (not NA)
  keep_mask <- rep(TRUE, n_original)

  for (metric in required_metrics) {
    if (metric %in% names(history)) {
      # Check column directly
      keep_mask <- keep_mask & !is.na(history[[metric]])
    } else if ("metrics" %in% names(history) && is.list(history$metrics)) {
      # Check metrics list-column
      has_metric <- sapply(history$metrics, function(m) {
        if (is.null(m)) return(FALSE)
        metric %in% names(m) && !is.na(m[[metric]])
      })
      keep_mask <- keep_mask & has_metric
    } else {
      # Metric not available anywhere - drop all rows
      if (verbose) {
        cat(sprintf("  Warning: Metric '%s' not found in history - dropping all rows\n", metric))
      }
      keep_mask <- rep(FALSE, n_original)
      break
    }
  }

  filtered <- history[keep_mask, , drop = FALSE]
  n_retained <- nrow(filtered)

  if (verbose) {
    cat(sprintf("  Filtered history for new constraints: %d -> %d rows (%.1f%% retained)\n",
                n_original, n_retained, 100 * n_retained / n_original))
    if (n_retained == 0) {
      cat("  Warning: No rows have all required constraint metrics - Stage will use random init\n")
    }
  }

  list(
    history = if (n_retained > 0) filtered else NULL,
    n_original = n_original,
    n_retained = n_retained,
    metrics_checked = required_metrics
  )
}
