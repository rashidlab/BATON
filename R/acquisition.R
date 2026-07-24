#' Expected constrained improvement acquisition function
#'
#' @param unit_x candidate points on the unit hypercube: a numeric matrix
#'   (rows = candidates, named columns) or a list of named numeric vectors.
#' @param surrogates named list of fitted surrogate models (as returned by [fit_surrogates()]).
#' @param constraint_tbl tibble from [parse_constraints()].
#' @param objective name of the objective metric (character scalar).
#' @param best_feasible best observed value of the objective among feasible points
#'   (numeric scalar or `Inf` if none observed).
#' @param pred optional list of surrogate predictions (as from [predict_surrogates()]).
#'   When supplied, avoids recomputing GP predictions for acquisition scoring.
#' @param prob_feas optional vector of feasibility probabilities aligned with `unit_x`.
#'   If NULL, computed internally from `pred` and `constraint_tbl`.
#'
#' @return numeric vector of acquisition scores for each candidate.
#' @export
acq_eci <- function(unit_x,
                    surrogates,
                    constraint_tbl,
                    objective,
                    best_feasible,
                    pred = NULL,
                    prob_feas = NULL) {
  if (!objective %in% names(surrogates)) {
    stop("Objective surrogate not available.", call. = FALSE)
  }
  if (is.null(pred)) {
    pred <- predict_surrogates(surrogates, unit_x)
  }

  obj_pred <- pred[[objective]]
  mu_obj <- obj_pred$mean
  sd_obj <- obj_pred$sd

  metric_names <- names(pred)

  # PERFORMANCE: Use vectorized batch computation instead of per-candidate loop
  if (is.null(prob_feas)) {
    prob_feas <- prob_feasibility_batch(pred, constraint_tbl)
  }

  has_feasible <- is.finite(best_feasible)

  if (!has_feasible) {
    # No feasible solution yet
    # Strategy: Minimize expected constraint violation weighted by feasibility probability
    violations <- compute_expected_violation(pred, constraint_tbl, metric_names)

    # Normalize violations to [0, 1] range
    max_viol <- max(violations, na.rm = TRUE)
    if (is.finite(max_viol) && max_viol > 0) {
      normalized_viol <- violations / max_viol
    } else {
      normalized_viol <- violations
    }

    # Acquisition = (1 - violation) + exploration bonus, weighted by feasibility
    # Scale to be comparable to EI values
    acq <- (1 - normalized_viol + 0.3 * sd_obj) * (0.3 + 0.7 * prob_feas)
    return(acq)
  }

  # Standard ECI: EI × P(feasible)
  ei <- compute_ei(mu_obj, sd_obj, best_feasible)
  ei * prob_feas
}

#' Quasi Expected Hypervolume Improvement (constraint-aware)
#'
#' @description
#' \strong{Note:} This function is currently a placeholder that calls
#' [acq_eci()]. A full qEHVI implementation with proper hypervolume computation
#' is planned for a future release.
#'
#' For batch size greater than one, the function uses a sequential (Kriging
#' believer) approximation via constrained expected improvement.
#'
#' @inheritParams acq_eci
#' @return numeric vector of acquisition scores.
#' @keywords internal
acq_qehvi <- function(unit_x,
                      surrogates,
                      constraint_tbl,
                      objective,
                      best_feasible,
                      pred = NULL,
                      prob_feas = NULL) {
  warning("acq_qehvi is currently an alias for acq_eci. Full qEHVI implementation is planned for a future release.",
          call. = FALSE)
  acq_eci(unit_x, surrogates, constraint_tbl, objective, best_feasible,
          pred = pred, prob_feas = prob_feas)
}

#' @keywords internal
compute_ei <- function(mu, sd, best_feasible) {
  ei <- numeric(length(mu))
  has_feasible <- is.finite(best_feasible)
  if (!has_feasible) {
    # No feasible solution yet - use exploration weighted by uncertainty
    # Add small epsilon to prevent division by zero
    return(pmax(sd + 1e-10, 0))
  }
  improvement <- best_feasible - mu
  positive <- sd > 1e-10  # Add epsilon for numerical stability
  z <- improvement[positive] / (sd[positive] + 1e-10)
  phi <- stats::dnorm(z)
  Phi <- stats::pnorm(z)
  ei[positive] <- improvement[positive] * Phi + (sd[positive] + 1e-10) * phi
  ei[!positive] <- pmax(improvement[!positive], 0)
  ei
}

#' Compute expected constraint violation
#'
#' For each candidate point, computes the expected magnitude of constraint
#' violations under the GP posterior. Used to guide search toward feasible
#' region when no feasible solutions have been found yet.
#'
#' @param pred list of predictions from predict_surrogates (contains mean and sd for each metric)
#' @param constraint_tbl tibble from parse_constraints
#' @param metric_names character vector of metric names
#'
#' @return numeric vector of expected violations for each candidate
#' @keywords internal
compute_expected_violation <- function(pred, constraint_tbl, metric_names) {
  n_candidates <- length(pred[[1]]$mean)
  violations <- numeric(n_candidates)

  # PERFORMANCE: Loop over constraints (typically 2-4), vectorize across all candidates
  for (j in seq_len(nrow(constraint_tbl))) {
    metric <- constraint_tbl$metric[j]
    direction <- constraint_tbl$direction[j]
    threshold <- constraint_tbl$threshold[j]

    if (!metric %in% metric_names) next

    # Extract mu and sd vectors for ALL candidates at once
    mu_vec <- pred[[metric]]$mean
    sd_vec <- pred[[metric]]$sd + 1e-10  # Add epsilon once

    # Expected violation magnitude E[max(0, gap)] for a Gaussian predictive
    # distribution is sd * (z * Phi(z) + phi(z)), where z is the standardized
    # signed distance into the infeasible region. This grows without bound as a
    # point becomes more infeasible and decays to 0 well inside the feasible
    # region - the property needed to steer search toward the boundary. (The
    # previous code reduced to sd * phi(z), which peaks AT the boundary and
    # vanishes for deeply infeasible points.)
    if (direction == "ge") {
      # Constraint: metric >= threshold; violation when metric < threshold
      z <- (threshold - mu_vec) / sd_vec
    } else {  # "le"
      # Constraint: metric <= threshold; violation when metric > threshold
      z <- (mu_vec - threshold) / sd_vec
    }
    violations <- violations + sd_vec * (z * stats::pnorm(z) + stats::dnorm(z))
  }

  violations
}

#' Select diverse batch using local penalization
#'
#' Implements the local penalization strategy of González et al. (2016).
#' Iteratively selects points by penalizing acquisition near previously
#' selected points, ensuring spatial diversity in batch selection.
#'
#' @param candidates list of candidate points (unit scale)
#' @param acq_scores numeric vector of acquisition values for each candidate
#' @param q batch size (number of points to select)
#' @param lipschitz Lipschitz constant for penalization. Higher values enforce
#'   greater diversity. Default estimated from typical BO landscapes.
#'
#' @return integer vector of indices of selected candidates
#' @keywords internal
#'
#' @references
#' González, J., Dai, Z., Hennig, P., & Lawrence, N. (2016).
#' Batch Bayesian Optimization via Local Penalization. AISTATS.
select_batch_local_penalization <- function(candidates, acq_scores, q,
                                             lipschitz = 10) {
  n_candidates <- length(acq_scores)
  if (q >= n_candidates) {
    return(seq_len(n_candidates))
  }

  selected_indices <- integer(q)
  penalized_scores <- acq_scores
  # Task C5: the candidate pool arrives as a matrix; keep the list path for
  # back-compat with direct callers.
  if (is.matrix(candidates)) {
    candidates_matrix <- candidates
  } else {
    n_dims <- length(candidates[[1]])
    candidates_matrix <- matrix(NA_real_, nrow = n_candidates, ncol = n_dims)
    for (i in seq_len(n_candidates)) {
      candidates_matrix[i, ] <- as.numeric(candidates[[i]])
    }
  }

  for (i in seq_len(q)) {
    # Select point with highest penalized acquisition
    best_idx <- which.max(penalized_scores)
    selected_indices[i] <- best_idx

    if (i < q) {
      # Penalize acquisition near selected point
      selected_point <- candidates_matrix[best_idx, , drop = FALSE]

      # Compute distances to selected point
      distances <- compute_distances(candidates_matrix, selected_point)

      # Penalization function (Gonzalez et al. 2016): max(0, acq_best - L * r)
      # where r is distance and L is the Lipschitz constant. The penalty is
      # LARGEST at the selected point (r = 0) and decays to 0 at radius
      # acq_best / L, so candidates NEAR an already-selected design are
      # suppressed and the batch spreads out. (The previous form used
      # L*r - acq_best, which penalized FAR candidates and clustered the batch.)
      penalty <- pmax(0, penalized_scores[best_idx] - lipschitz * distances)

      # Apply penalty
      penalized_scores <- penalized_scores - penalty

      # Ensure we don't select same point again
      penalized_scores[best_idx] <- -Inf
    }
  }

  selected_indices
}

#' Compute Euclidean distances from points to reference
#'
#' @param points n × d matrix of points
#' @param reference 1 × d matrix (or vector) of reference point
#'
#' @return numeric vector of n distances
#' @keywords internal
compute_distances <- function(points, reference) {
  # Ensure reference is a matrix
  if (is.vector(reference)) {
    reference <- matrix(reference, nrow = 1)
  }

  # Compute differences
  diff <- sweep(points, 2, reference[1, ], "-")

  # Euclidean distance
  sqrt(rowSums(diff^2))
}

#' Estimate Lipschitz constant from GP lengthscales
#'
#' Uses GP lengthscales to estimate a reasonable Lipschitz constant
#' for local penalization. Conservative estimate to ensure diversity.
#'
#' @param surrogates list of GP models
#' @param objective name of objective metric
#'
#' @return numeric scalar, Lipschitz constant estimate
#' @keywords internal
estimate_lipschitz <- function(surrogates, objective) {
  model <- surrogates[[objective]]

  tryCatch({
    lengthscales <- NULL

    # Extract lengthscales based on model type
    if (inherits(model, "km")) {
      # DiceKriging model (S4 slots)
      lengthscales <- model@covariance@range.val
    } else if (inherits(model, c("hetGP", "homGP"))) {
      # hetGP/homGP models (list elements)
      lengthscales <- model$theta
    }

    if (!is.null(lengthscales) && length(lengthscales) > 0 &&
        all(is.finite(lengthscales)) && all(lengthscales > 0)) {
      # Lipschitz constant inversely related to lengthscale
      # Use minimum lengthscale (most sensitive direction)
      L <- 1 / min(lengthscales)
      # Conservative factor: multiply by 2
      return(L * 2)
    }

    # Default if extraction fails
    return(10)
  }, error = function(e) {
    # If extraction fails, use default
    return(10)
  })
}
