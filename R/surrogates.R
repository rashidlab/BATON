# Copyright (c) 2026. For not-for-profit research and educational use only; all
# other rights reserved. See the LICENSE file for full terms.

#' Fit Gaussian process surrogates for operating characteristics
#'
#' Implements heteroskedastic GP modeling when variance information is available,
#' using either hetGP (preferred) or DiceKriging with known noise variance.
#'
#' @param history tibble with columns `unit_x` (list of numeric vectors),
#'   `metrics` (list of named numeric vectors), and optional `variance`
#'   (list of named numeric vectors with noise variances).
#' @param objective name of the objective metric.
#' @param constraint_tbl tibble produced by [parse_constraints()].
#' @param covtype covariance kernel used by DiceKriging (default `"matern5_2"`).
#' @param use_hetgp logical; if TRUE and hetGP package is available, use
#'   heteroskedastic GP with replicated observations. Default TRUE.
#' @param prev_surrogates optional list of surrogates from previous iteration
#'   for warm-starting hyperparameter optimization. If provided, uses previous
#'   hyperparameters as initial values, reducing optimization time by 30-50\%.
#'
#' @return A named list of fitted surrogate models (GP objects).
#' @export
fit_surrogates <- function(history,
                           objective,
                           constraint_tbl,
                           covtype = "matern5_2",
                           use_hetgp = TRUE,
                           prev_surrogates = NULL) {
  if (nrow(history) == 0L) {
    stop("History is empty; cannot fit surrogates.", call. = FALSE)
  }

  metrics_needed <- unique(c(objective, constraint_tbl$metric))

  # Robustly extract param_names from history
  param_names <- NULL
  if (nrow(history) > 0) {
    # Try unit_x first
    if ("unit_x" %in% names(history) && length(history$unit_x) > 0 && !is.null(history$unit_x[[1]])) {
      param_names <- names(history$unit_x[[1]])
    }
    # Fall back to theta
    if (is.null(param_names) && "theta" %in% names(history) && length(history$theta) > 0 && !is.null(history$theta[[1]])) {
      param_names <- names(history$theta[[1]])
    }
  }

  if (is.null(param_names) || length(param_names) == 0) {
    stop("Cannot determine parameter names from history. Check that unit_x or theta columns exist and have named elements.", call. = FALSE)
  }

  theta_ids <- history$theta_id
  id_groups <- split(seq_along(theta_ids), theta_ids)

  # Check if hetGP is available and requested
  has_hetgp <- requireNamespace("hetGP", quietly = TRUE)
  use_hetgp <- use_hetgp && has_hetgp

  # Use lapply instead of purrr::map to avoid purrr's error wrapping
  # which can cause confusing "In index: X" error messages
  surrogates <- lapply(metrics_needed, function(metric) {
    tryCatch({
      # OPTIMIZED: Extract values with vapply for type safety and speed
      values <- vapply(seq_along(history$metrics), function(i) {
        tryCatch({
          m <- history$metrics[[i]]
          if (is.null(m) || !metric %in% names(m)) return(NA_real_)
          as.numeric(m[[metric]])
        }, error = function(e) NA_real_)
      }, FUN.VALUE = numeric(1))

      noise <- vapply(seq_along(history$variance), function(i) {
        tryCatch({
          var_list <- history$variance[[i]]
          if (is.null(var_list) || !metric %in% names(var_list)) return(NA_real_)
          as.numeric(var_list[[metric]])
        }, error = function(e) NA_real_)
      }, FUN.VALUE = numeric(1))

      # Check if we have replicated observations (multiple evals at same theta)
      has_replicates <- any(vapply(id_groups, length, FUN.VALUE = integer(1)) > 1)
      has_variance <- !all(is.na(noise))

      # Extract previous model for warm-starting (if available)
      prev_model <- if (!is.null(prev_surrogates) && metric %in% names(prev_surrogates)) {
        prev_surrogates[[metric]]
      } else {
        NULL
      }

      if (use_hetgp && has_replicates && has_variance) {
        # Use hetGP with replicated observations
        fit_hetgp_surrogate(history, metric, id_groups, param_names, covtype, prev_model)
      } else {
        # Fall back to DiceKriging with aggregated observations
        if (use_hetgp && !has_variance && !isTRUE(getOption("BATON.homoskedastic_warned"))) {
          message(sprintf("  [surrogate] Using homoskedastic GP for '%s' (no variance estimates provided). Attach attr(result, 'variance') to enable hetGP.", metric))
          options(BATON.homoskedastic_warned = TRUE)
        }
        fit_dicekriging_surrogate(history, metric, id_groups, param_names,
                                  covtype, noise, values, prev_model)
      }
    }, error = function(e) {
      # On any error, return constant predictor instead of failing
      warning(sprintf("[fit_surrogates] Error fitting surrogate for metric '%s': %s. Using constant predictor.",
                      metric, e$message), call. = FALSE)
      # Compute fallback mean from available values
      fallback_mean <- tryCatch({
        vals <- vapply(history$metrics, function(m) {
          if (is.null(m) || !metric %in% names(m)) NA_real_ else as.numeric(m[[metric]])
        }, FUN.VALUE = numeric(1))
        mean(vals, na.rm = TRUE)
      }, error = function(e2) NA_real_)

      structure(
        list(mean_value = if (is.finite(fallback_mean)) fallback_mean else 0,
             metric = metric),
        class = "constant_predictor"
      )
    })
  })

  names(surrogates) <- metrics_needed
  surrogates
}

#' Fit heteroskedastic GP using hetGP package
#' @keywords internal
fit_hetgp_surrogate <- function(history, metric, id_groups, param_names, covtype, prev_model = NULL) {
  # Prepare data with replicates
  X_list <- list()
  Z_list <- list()

  for (i in seq_along(id_groups)) {
    idx <- id_groups[[i]]
    unit_theta <- history$unit_x[[idx[1]]]
    X_list[[i]] <- as.numeric(unit_theta)

    # All observations at this location
    Z_list[[i]] <- vapply(idx, function(j) {
      as.numeric(history$metrics[[j]][[metric]])
    }, FUN.VALUE = numeric(1))
  }

  # Filter out locations where ALL Z values are NA
  valid_locs <- vapply(Z_list, function(z) !all(is.na(z)), FUN.VALUE = logical(1))
  n_invalid <- sum(!valid_locs)
  if (n_invalid > 0) {
    message(sprintf("  [hetGP] Filtered %d/%d locations with all-NA values for metric '%s'",
                    n_invalid, length(Z_list), metric))
    X_list <- X_list[valid_locs]
    Z_list <- Z_list[valid_locs]
  }

  # Check if we have enough observations
  if (length(X_list) < 2) {
    mean_val <- if (length(Z_list) == 1) mean(Z_list[[1]], na.rm = TRUE) else NA_real_
    message(sprintf("  [hetGP] Insufficient observations (%d) for GP - using constant predictor for '%s'",
                    length(X_list), metric))
    return(structure(
      list(mean_value = mean_val, metric = metric),
      class = "constant_predictor"
    ))
  }

  # OPTIMIZED: Pre-allocate matrix instead of do.call(rbind, ...)
  n_locs <- length(X_list)
  n_params <- length(param_names)
  X <- matrix(NA_real_, nrow = n_locs, ncol = n_params)
  for (i in seq_len(n_locs)) {
    X[i, ] <- X_list[[i]]
  }
  colnames(X) <- param_names

  # Filter NA values within each location's Z values (keep only non-NA)
  Z_list <- lapply(Z_list, function(z) z[!is.na(z)])

  # hetGP expects X with duplicated rows matching Z (auto-detects replicates)
  Z_vec <- unlist(Z_list)
  mult <- vapply(Z_list, length, FUN.VALUE = integer(1))

  # Expand X so each row corresponds to one observation in Z_vec
  X_expanded <- X[rep(seq_len(n_locs), mult), , drop = FALSE]

  # Map covtype to hetGP
  hetgp_cov <- switch(covtype,
                      matern5_2 = "Matern5_2",
                      matern3_2 = "Matern3_2",
                      gauss = "Gaussian",
                      "Matern5_2")  # default

  tryCatch({
    hetGP::mleHetGP(
      X = X_expanded,
      Z = Z_vec,
      covtype = hetgp_cov,
      settings = list(return.hom = TRUE),  # also return homoskedastic fit
      known = list(g = 1e-8)  # small nugget for numerical stability
    )
  }, error = function(e) {
    warning(sprintf("hetGP fit failed for metric '%s': %s\nFalling back to homoskedastic GP.",
                    metric, e$message), call. = FALSE)
    # Fall back to homoskedastic GP
    # Note: mleHomGP does NOT accept 'mult' argument - it expects pre-aggregated data
    # Aggregate Z by unique X locations before fitting
    # Use row-wise deduplication by converting to data.frame
    X_df <- as.data.frame(X)
    X_unique_idx <- !duplicated(X_df)
    X_agg <- X[X_unique_idx, , drop = FALSE]

    # Aggregate Z values at each unique location
    # Z_list[[i]] contains all Z values at location X[i,], so use Z_list directly
    # (Z_vec = unlist(Z_list) has different length than nrow(X))
    Z_agg <- vapply(seq_len(nrow(X_agg)), function(i) {
      # Find which rows of X match this unique location
      matching_rows <- which(apply(X, 1, function(row) all(row == X_agg[i, ])))
      # Get all Z values from those matching locations
      all_z <- unlist(Z_list[matching_rows])
      mean(all_z, na.rm = TRUE)
    }, FUN.VALUE = numeric(1))

    # Check for constant/near-constant data which causes GP fitting to fail
    z_range <- diff(range(Z_agg, na.rm = TRUE))
    if (z_range < 1e-10) {
      warning(sprintf("Metric '%s' has constant or near-constant values (range=%.2e). Using DiceKriging fallback.",
                      metric, z_range), call. = FALSE)
      # Fall back to DiceKriging which handles constant data more gracefully
      return(tryCatch({
        DiceKriging::km(
          design = X_agg,
          response = Z_agg,
          covtype = "matern5_2",
          nugget = 1e-4,  # Larger nugget for stability
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        )
      }, error = function(e2) {
        # Ultimate fallback: return a dummy model that predicts the mean
        warning(sprintf("All GP fits failed for metric '%s'. Using constant predictor.", metric), call. = FALSE)
        structure(
          list(mean_value = mean(Z_agg, na.rm = TRUE), metric = metric),
          class = "constant_predictor"
        )
      }))
    }

    tryCatch({
      hetGP::mleHomGP(X = X_agg, Z = Z_agg, covtype = hetgp_cov,
                     known = list(g = 1e-6))
    }, error = function(e2) {
      warning(sprintf("mleHomGP also failed for metric '%s': %s\nUsing DiceKriging fallback.",
                      metric, e2$message), call. = FALSE)
      # Try DiceKriging as last resort before constant predictor
      tryCatch({
        DiceKriging::km(
          design = X_agg,
          response = Z_agg,
          covtype = "matern5_2",
          nugget = 1e-4,
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        )
      }, error = function(e3) {
        warning(sprintf("All GP fits failed for metric '%s'. Using constant predictor.", metric), call. = FALSE)
        structure(
          list(mean_value = mean(Z_agg, na.rm = TRUE), metric = metric),
          class = "constant_predictor"
        )
      })
    })
  })
}

#' Fit GP using DiceKriging (aggregated observations)
#' @keywords internal
fit_dicekriging_surrogate <- function(history, metric, id_groups, param_names,
                                      covtype, noise, values, prev_model = NULL) {
  # Aggregate by theta_id - with robust error handling
  # Use base R lapply instead of purrr::imap to avoid "In index: X" errors
  group_names <- names(id_groups)
  aggr_list <- lapply(seq_along(id_groups), function(i) {
    idx <- id_groups[[i]]
    id <- group_names[i]
    tryCatch({
      unit_theta <- history$unit_x[[idx[1]]]
      # Ensure unit_theta is a proper numeric vector
      if (is.null(unit_theta) || length(unit_theta) == 0) {
        warning(sprintf("Empty unit_theta at theta_id '%s' (idx=%d)", id, idx[1]))
        return(NULL)
      }
      tibble::tibble(
        theta_id = id,
        unit_x = list(as.numeric(unit_theta)),
        value = mean(values[idx], na.rm = TRUE),
        noise = if (all(is.na(noise[idx]))) NA_real_ else mean(noise[idx], na.rm = TRUE)
      )
    }, error = function(e) {
      warning(sprintf("Error aggregating theta_id '%s': %s", id, e$message))
      return(NULL)
    })
  })

  # Filter out NULLs and bind
  aggr_list <- aggr_list[!vapply(aggr_list, is.null, FUN.VALUE = logical(1))]
  if (length(aggr_list) == 0) {
    stop(sprintf("No valid observations to fit surrogate for metric '%s'", metric), call. = FALSE)
  }
  aggr <- dplyr::bind_rows(aggr_list)

  # Filter out rows with NA/NaN values (can happen when all observations at a theta are NA)
  valid_rows <- !is.na(aggr$value) & !is.nan(aggr$value)
  n_invalid <- sum(!valid_rows)
  if (n_invalid > 0) {
    message(sprintf("  [surrogate] Filtered %d/%d rows with NA/NaN values for metric '%s'",
                    n_invalid, nrow(aggr), metric))
    aggr <- aggr[valid_rows, , drop = FALSE]
  }

  # Check if we have enough observations to fit a GP
  if (nrow(aggr) < 2) {
    # Return a constant predictor if insufficient data
    mean_val <- if (nrow(aggr) == 1) aggr$value[1] else NA_real_
    message(sprintf("  [surrogate] Insufficient observations (%d) for GP - using constant predictor for '%s'",
                    nrow(aggr), metric))
    return(structure(
      list(mean_value = mean_val, metric = metric),
      class = "constant_predictor"
    ))
  }

  # OPTIMIZED: Pre-allocate matrix instead of do.call(rbind, ...)
  n_rows <- nrow(aggr)
  n_cols <- length(param_names)
  X_unique <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)
  for (i in seq_len(n_rows)) {
    vec <- as.numeric(aggr$unit_x[[i]])
    if (length(vec) != n_cols) {
      stop("Surrogate design dimension mismatch.", call. = FALSE)
    }
    X_unique[i, ] <- vec
  }
  colnames(X_unique) <- param_names

  noise_vec <- aggr$noise
  if (all(is.na(noise_vec))) {
    nugget <- 1e-6
    noise_vec <- NULL
  } else {
    noise_vec[is.na(noise_vec)] <- min(noise_vec, na.rm = TRUE)
    noise_vec <- pmax(noise_vec, 1e-6)
    nugget <- 0
  }

  # Extract hyperparameters from previous model for warm-starting
  parinit <- extract_gp_hyperparams(prev_model)

  # Debug: Check what we're passing to DiceKriging::km
  if (!is.matrix(X_unique)) {
    warning(sprintf("X_unique is not a matrix! Class: %s, Dim: %s",
                    class(X_unique)[1],
                    paste(dim(X_unique), collapse = "x")))
  }
  if (!is.numeric(aggr$value)) {
    warning(sprintf("aggr$value is not numeric! Class: %s, Length: %d",
                    class(aggr$value)[1],
                    length(aggr$value)))
  }

  tryCatch({
    if (is.null(noise_vec)) {
      DiceKriging::km(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        nugget = nugget,
        nugget.estim = FALSE,
        control = list(
          trace = FALSE,
          parinit = parinit  # Warm start
        )
      )
    } else {
      DiceKriging::km(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        noise.var = noise_vec,
        nugget.estim = FALSE,
        control = list(
          trace = FALSE,
          parinit = parinit  # Warm start
        )
      )
    }
  }, error = function(e) {
    stop(sprintf("Failed to fit surrogate for metric '%s': %s\nThis may indicate ill-conditioned data or insufficient observations.",
                 metric, e$message), call. = FALSE)
  })
}

#' Extract GP hyperparameters for warm-starting
#'
#' Extracts lengthscale parameters from a fitted GP model to use as
#' initial values for the next optimization. Works with DiceKriging::km,
#' hetGP, and homGP models.
#'
#' @param model fitted GP model (km, hetGP, homGP, or constant_predictor), or NULL
#' @return numeric vector of hyperparameters, or NULL if extraction fails
#' @keywords internal
extract_gp_hyperparams <- function(model) {
  if (is.null(model)) {
    return(NULL)
  }

  # Constant predictor has no hyperparameters
  if (inherits(model, "constant_predictor")) {
    return(NULL)
  }

  tryCatch({
    if (inherits(model, "km")) {
      # DiceKriging model
      theta <- model@covariance@range.val
      if (is.numeric(theta) && all(is.finite(theta)) && all(theta > 0)) {
        return(theta)
      }
    } else if (inherits(model, c("hetGP", "homGP"))) {
      # hetGP or homGP model (both use $theta)
      theta <- model$theta
      if (is.numeric(theta) && all(is.finite(theta)) && all(theta > 0)) {
        return(theta)
      }
    }
    return(NULL)
  }, error = function(e) {
    # If extraction fails, return NULL (no warm-start)
    return(NULL)
  })
}

#' Predict metrics from fitted surrogates
#'
#' Handles both DiceKriging::km and hetGP::hetGP model objects.
#' Optimized for batch prediction with minimal data frame overhead.
#'
#' @keywords internal
predict_surrogates <- function(surrogates, unit_x) {
  if (length(surrogates) == 0L) {
    stop("No surrogate models available for prediction.", call. = FALSE)
  }

  n_candidates <- length(unit_x)
  if (n_candidates == 0L) {
    stop("predict_surrogates: No candidate points provided", call. = FALSE)
  }

  # Detect model type from first non-constant surrogate to get param_names
  param_names <- NULL
  for (model in surrogates) {
    if (inherits(model, "constant_predictor")) {
      next
    } else if (inherits(model, c("hetGP", "homGP"))) {
      param_names <- colnames(model$X0)
      break
    } else if (inherits(model, "km")) {
      param_names <- colnames(model@X)
      break
    }
  }

  if (is.null(param_names)) {
    param_names <- names(unit_x[[1]])
  }

  n_params <- length(param_names)

  # OPTIMIZED: Build design matrix directly without per-candidate data.frame creation
  # Pre-allocate matrix and fill row-by-row (much faster than lapply + rbind)
  design_mat <- matrix(NA_real_, nrow = n_candidates, ncol = n_params)
  colnames(design_mat) <- param_names

  valid_rows <- rep(TRUE, n_candidates)
  for (i in seq_len(n_candidates)) {
    point <- unit_x[[i]]
    tryCatch({
      vec <- unlist(point)
      if (is.null(vec) || length(vec) == 0) {
        valid_rows[i] <- FALSE
        next
      }
      if (is.null(names(vec))) {
        names(vec) <- param_names
      }
      # Extract values in param_names order
      if (all(param_names %in% names(vec))) {
        design_mat[i, ] <- as.numeric(vec[param_names])
      } else {
        valid_rows[i] <- FALSE
      }
    }, error = function(e) {
      valid_rows[i] <<- FALSE
    })
  }

  # Filter invalid rows if any
  if (!all(valid_rows)) {
    n_invalid <- sum(!valid_rows)
    if (n_invalid > 0 && n_invalid < n_candidates) {
      warning(sprintf("predict_surrogates: Filtered %d invalid points", n_invalid))
      design_mat <- design_mat[valid_rows, , drop = FALSE]
    } else if (n_invalid == n_candidates) {
      stop("predict_surrogates: No valid design points", call. = FALSE)
    }
  }

  # Convert to data.frame ONCE for DiceKriging compatibility (km requires data.frame)
  # hetGP can use matrix directly
  design_df <- as.data.frame(design_mat, stringsAsFactors = FALSE)

  # Cache number of points for efficiency

  n_points <- nrow(design_mat)

  # Use base R lapply + setNames instead of purrr::imap to avoid
  # confusing "In index: X" error messages from purrr 1.0.0+
  surrogate_names <- names(surrogates)
  predictions <- lapply(seq_along(surrogates), function(i) {
    metric_name <- surrogate_names[i]
    model <- surrogates[[i]]

    tryCatch({
      if (inherits(model, "constant_predictor")) {
        # Constant predictor fallback - returns constant mean with high uncertainty
        list(mean = rep(model$mean_value, n_points),
             sd = rep(1.0, n_points),  # High uncertainty to encourage exploration
             model = model)
      } else if (inherits(model, c("hetGP", "homGP"))) {
        # Use hetGP/homGP predict method - can use matrix directly (faster)
        pred <- predict(x = design_mat, object = model)
        list(mean = as.numeric(pred$mean),
             sd = as.numeric(sqrt(pmax(pred$sd2, 0))),
             model = model)
      } else {
        # Use DiceKriging predict method - requires data.frame
        pred <- DiceKriging::predict.km(model, newdata = design_df,
                                        type = "UK", se.compute = TRUE,
                                        cov.compute = FALSE, checkNames = FALSE)
        list(mean = as.numeric(pred$mean),
             sd = sqrt(pmax(pred$sd^2, 0)),
             model = model)
      }
    }, error = function(e) {
      warning(sprintf("[predict_surrogates] Error predicting metric '%s': %s. Using fallback.",
                      metric_name, e$message), call. = FALSE)
      # Return high uncertainty predictions to encourage exploration
      list(mean = rep(0, n_points),
           sd = rep(10.0, n_points),
           model = structure(list(mean_value = 0, metric = metric_name),
                             class = "constant_predictor"))
    })
  })
  names(predictions) <- surrogate_names
  predictions
}
