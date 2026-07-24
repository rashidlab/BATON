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
#' @param fit_seed optional integer base seed for reproducible hyperparameter
#'   fitting. When supplied, each metric's fit runs under a deterministic seed
#'   (`fit_seed` + metric index), making results identical regardless of core
#'   count when fits are parallelised via `options(BATON.cores)`. `NULL` (default)
#'   leaves the RNG untouched.
#'
#' @return A named list of fitted surrogate models (GP objects).
#' @export
fit_surrogates <- function(history,
                           objective,
                           constraint_tbl,
                           covtype = "matern5_2",
                           use_hetgp = TRUE,
                           prev_surrogates = NULL,
                           fit_seed = NULL) {
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

  # The m per-metric fits are independent. They run concurrently when the caller
  # opts in via options(BATON.cores = k) on Unix; the default k=1 keeps the serial
  # path. Reproducibility is guaranteed regardless of core count by giving each fit
  # a deterministic seed (fit_seed + metric index) so km/hetGP random starts do not
  # depend on RNG-stream position or worker forking (see run_seeded()).
  metric_index <- stats::setNames(seq_along(metrics_needed), metrics_needed)
  .cores <- effective_cores()
  .map <- if (.cores > 1L) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = .cores)
  } else {
    lapply
  }

  # Use an explicit mapper (not purrr::map) to avoid purrr's error wrapping
  # which can cause confusing "In index: X" error messages.
  surrogates <- .map(metrics_needed, function(metric) {
    metric_seed <- if (!is.null(fit_seed)) fit_seed + metric_index[[metric]] else NULL
    tryCatch({
      # OPTIMIZED: Extract values with vapply for type safety and speed
      values <- vapply(seq_along(history$metrics), function(i) {
        tryCatch({
          m <- history$metrics[[i]]
          if (is.null(m) || !metric %in% names(m)) return(NA_real_)
          as.numeric(m[[metric]])
        }, error = function(e) NA_real_)
      }, FUN.VALUE = numeric(1))

      # variance is optional: a history without the column must take the
      # homoskedastic nugget path (all-NA noise), not produce a length-0
      # vector that breaks grouped aggregation downstream.
      noise <- if (!"variance" %in% names(history)) {
        rep(NA_real_, nrow(history))
      } else {
        vapply(seq_along(history$variance), function(i) {
          tryCatch({
            var_list <- history$variance[[i]]
            if (is.null(var_list) || !metric %in% names(var_list)) return(NA_real_)
            as.numeric(var_list[[metric]])
          }, error = function(e) NA_real_)
        }, FUN.VALUE = numeric(1))
      }

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
        fit_hetgp_surrogate(history, metric, id_groups, param_names, covtype,
                            prev_model, fit_seed = metric_seed)
      } else {
        # Fall back to DiceKriging with aggregated observations
        if (use_hetgp && !has_variance && !isTRUE(getOption("BATON.homoskedastic_warned"))) {
          message(sprintf("  [surrogate] Using homoskedastic GP for '%s' (no variance estimates provided). Attach attr(result, 'variance') to enable hetGP.", metric))
          options(BATON.homoskedastic_warned = TRUE)
        }
        fit_dicekriging_surrogate(history, metric, id_groups, param_names,
                                  covtype, noise, values, prev_model,
                                  fit_seed = metric_seed)
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
fit_hetgp_surrogate <- function(history, metric, id_groups, param_names, covtype,
                                prev_model = NULL, fit_seed = NULL) {
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
    run_seeded(fit_seed, function() hetGP::mleHetGP(
      X = X_expanded,
      Z = Z_vec,
      covtype = hetgp_cov,
      settings = list(return.hom = TRUE),  # also return homoskedastic fit
      known = list(g = 1e-8)  # small nugget for numerical stability
    ))
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
        run_seeded(fit_seed, function() DiceKriging::km(
          design = X_agg,
          response = Z_agg,
          covtype = "matern5_2",
          nugget = 1e-4,  # Larger nugget for stability
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        ))
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
      run_seeded(fit_seed, function() hetGP::mleHomGP(X = X_agg, Z = Z_agg,
                     covtype = hetgp_cov, known = list(g = 1e-6)))
    }, error = function(e2) {
      warning(sprintf("mleHomGP also failed for metric '%s': %s\nUsing DiceKriging fallback.",
                      metric, e2$message), call. = FALSE)
      # Try DiceKriging as last resort before constant predictor
      tryCatch({
        run_seeded(fit_seed, function() DiceKriging::km(
          design = X_agg,
          response = Z_agg,
          covtype = "matern5_2",
          nugget = 1e-4,
          nugget.estim = FALSE,
          control = list(trace = FALSE)
        ))
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

#' Run a fitting expression under a fixed RNG seed, restoring the stream after
#'
#' GP hyperparameter optimizers use random starts (km cold-start, mleHetGP), so
#' without a fixed seed the result depends on RNG-stream position - which differs
#' between a serial lapply and forked mclapply workers, making calibration
#' core-count-dependent. This isolates each fit under a deterministic seed and
#' restores the caller's RNG state on exit, so results are identical regardless
#' of execution order or parallelism. When `seed` is NULL, behaves as a plain call.
#' @keywords internal
run_seeded <- function(seed, fn) {
  if (is.null(seed)) return(fn())
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }
  set.seed(seed)
  fn()
}

#' NA-aware grouped means in one matrix pass
#'
#' Task C6: replaces the per-theta one-row-tibble + `bind_rows` aggregation.
#' Matches `mean(x[group], na.rm = TRUE)` per group exactly, including the
#' NaN result for all-NA groups (a plain `rowsum` would propagate NA into
#' the surrogate inputs). Groups are returned ordered by sorted unique id,
#' the same order `split()` produces.
#'
#' @keywords internal
na_aware_group_mean <- function(x, ids) {
  ok <- !is.na(x)
  num <- rowsum(ifelse(ok, x, 0), ids, reorder = TRUE)
  den <- rowsum(as.numeric(ok), ids, reorder = TRUE)
  res <- as.numeric(num / den)  # 0/0 = NaN, matching mean(all-NA, na.rm=TRUE)
  names(res) <- rownames(num)
  res
}

#' Fit GP using DiceKriging (aggregated observations)
#' @keywords internal
fit_dicekriging_surrogate <- function(history, metric, id_groups, param_names,
                                      covtype, noise, values, prev_model = NULL,
                                      fit_seed = NULL) {
  # Task C6: aggregate by theta_id in one matrix pass (grouped means via
  # rowsum) instead of building one tibble per unique theta and bind_rows-ing
  # them (~1200 allocations per iteration on a mid-size history).
  group_names <- names(id_groups)
  n_groups <- length(group_names)
  n_cols <- length(param_names)

  agg_value <- na_aware_group_mean(values, history$theta_id)
  agg_noise <- na_aware_group_mean(noise, history$theta_id)
  # Defensive alignment: both na_aware_group_mean and split() order groups by
  # sorted unique id, but index by name so a divergence cannot misalign rows.
  agg_value <- agg_value[group_names]
  agg_noise <- agg_noise[group_names]

  # First-occurrence design point per group, with the previous per-group
  # guards (warn and drop the group rather than abort the whole fit).
  X_all <- matrix(NA_real_, nrow = n_groups, ncol = n_cols)
  ok_group <- rep(TRUE, n_groups)
  for (i in seq_len(n_groups)) {
    idx <- id_groups[[i]]
    unit_theta <- history$unit_x[[idx[1]]]
    if (is.null(unit_theta) || length(unit_theta) == 0) {
      warning(sprintf("Empty unit_theta at theta_id '%s' (idx=%d)",
                      group_names[i], idx[1]))
      ok_group[i] <- FALSE
      next
    }
    vec <- tryCatch(as.numeric(unit_theta), error = function(e) NULL)
    if (is.null(vec)) {
      warning(sprintf("Error aggregating theta_id '%s': non-numeric unit_x",
                      group_names[i]))
      ok_group[i] <- FALSE
      next
    }
    if (length(vec) != n_cols) {
      stop("Surrogate design dimension mismatch.", call. = FALSE)
    }
    X_all[i, ] <- vec
  }

  if (!any(ok_group)) {
    stop(sprintf("No valid observations to fit surrogate for metric '%s'", metric), call. = FALSE)
  }

  # Filter out groups with NA/NaN aggregate values (all observations NA)
  value_ok <- !is.na(agg_value) & !is.nan(agg_value)
  n_invalid <- sum(ok_group & !value_ok)
  if (n_invalid > 0) {
    message(sprintf("  [surrogate] Filtered %d/%d rows with NA/NaN values for metric '%s'",
                    n_invalid, sum(ok_group), metric))
  }
  keep <- ok_group & value_ok

  # Check if we have enough observations to fit a GP
  if (sum(keep) < 2) {
    # Return a constant predictor if insufficient data
    mean_val <- if (sum(keep) == 1) unname(agg_value[keep][1]) else NA_real_
    message(sprintf("  [surrogate] Insufficient observations (%d) for GP - using constant predictor for '%s'",
                    sum(keep), metric))
    return(structure(
      list(mean_value = mean_val, metric = metric),
      class = "constant_predictor"
    ))
  }

  X_unique <- X_all[keep, , drop = FALSE]
  colnames(X_unique) <- param_names
  aggr <- list(value = unname(agg_value[keep]))

  noise_vec <- unname(agg_noise[keep])
  if (all(is.na(noise_vec))) {
    nugget <- 1e-6
    noise_vec <- NULL
  } else {
    # Impute unknown (NA) noise with the LARGEST observed noise, not the
    # smallest: an unknown-variance observation should be trusted least by the
    # GP, not most. (min-imputation gave such points maximal weight.)
    noise_vec[is.na(noise_vec)] <- max(noise_vec, na.rm = TRUE)
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

  # Warm-start: parinit is a TOP-LEVEL km() argument, not a control entry
  # (km ignores control$parinit entirely). Only pass it when it is a valid
  # lengthscale vector of the right dimension; otherwise km cold-starts.
  km_control <- list(trace = FALSE)
  use_parinit <- !is.null(parinit) && is.numeric(parinit) &&
    length(parinit) == ncol(X_unique) &&
    all(is.finite(parinit)) && all(parinit > 0)

  tryCatch({
    if (is.null(noise_vec)) {
      km_args <- list(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        nugget = nugget,
        nugget.estim = FALSE,
        control = km_control
      )
      if (use_parinit) km_args$parinit <- parinit
      run_seeded(fit_seed, function() do.call(DiceKriging::km, km_args))
    } else {
      km_args <- list(
        design = X_unique,
        response = aggr$value,
        covtype = covtype,
        noise.var = noise_vec,
        nugget.estim = FALSE,
        control = km_control
      )
      if (use_parinit) km_args$parinit <- parinit
      run_seeded(fit_seed, function() do.call(DiceKriging::km, km_args))
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
#' `unit_x` is either a numeric matrix (rows = candidates, columns in the
#' surrogate's parameter order, or named via `colnames`) or a list of named
#' numeric vectors. The matrix path (Task C5) skips the per-candidate
#' rebuild entirely and is what the acquisition hot loop uses. Invalid
#' candidates are an ERROR on both paths: silently dropping rows would
#' misalign the returned predictions with the caller's candidate indices
#' (acquisition scores, batch selection).
#'
#' @keywords internal
predict_surrogates <- function(surrogates, unit_x) {
  if (length(surrogates) == 0L) {
    stop("No surrogate models available for prediction.", call. = FALSE)
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

  if (is.matrix(unit_x)) {
    # Task C5 fast path: candidates already in matrix form.
    n_candidates <- nrow(unit_x)
    if (n_candidates == 0L) {
      stop("predict_surrogates: No candidate points provided", call. = FALSE)
    }
    if (is.null(param_names)) {
      param_names <- colnames(unit_x)
    }
    design_mat <- unit_x
    if (!is.null(param_names)) {
      if (!is.null(colnames(design_mat))) {
        if (!all(param_names %in% colnames(design_mat))) {
          stop(sprintf(
            "predict_surrogates: candidate matrix is missing column(s): %s",
            paste(setdiff(param_names, colnames(design_mat)), collapse = ", ")
          ), call. = FALSE)
        }
        design_mat <- design_mat[, param_names, drop = FALSE]
      } else if (ncol(design_mat) == length(param_names)) {
        colnames(design_mat) <- param_names
      } else {
        stop(sprintf(
          "predict_surrogates: candidate matrix has %d columns; expected %d (%s)",
          ncol(design_mat), length(param_names),
          paste(param_names, collapse = ", ")
        ), call. = FALSE)
      }
    }
    if (!all(is.finite(design_mat))) {
      bad <- which(!apply(is.finite(design_mat), 1, all))
      stop(sprintf(
        "predict_surrogates: non-finite candidate coordinates in row(s): %s",
        paste(utils::head(bad, 5), collapse = ", ")
      ), call. = FALSE)
    }
  } else {
    n_candidates <- length(unit_x)
    if (n_candidates == 0L) {
      stop("predict_surrogates: No candidate points provided", call. = FALSE)
    }
    if (is.null(param_names)) {
      param_names <- names(unit_x[[1]])
    }
    n_params <- length(param_names)

    # Build design matrix directly without per-candidate data.frame creation
    design_mat <- matrix(NA_real_, nrow = n_candidates, ncol = n_params)
    colnames(design_mat) <- param_names

    invalid_rows <- integer(0)
    for (i in seq_len(n_candidates)) {
      point <- unit_x[[i]]
      ok <- tryCatch({
        vec <- unlist(point)
        if (is.null(vec) || length(vec) == 0) {
          FALSE
        } else {
          if (is.null(names(vec))) {
            names(vec) <- param_names
          }
          # Extract values in param_names order
          if (all(param_names %in% names(vec))) {
            design_mat[i, ] <- as.numeric(vec[param_names])
            TRUE
          } else {
            FALSE
          }
        }
      }, error = function(e) FALSE)
      if (!ok) invalid_rows <- c(invalid_rows, i)
    }

    # Invalid candidates are an error, never a silent drop: filtering rows
    # would misalign predictions with the caller's candidate indices.
    if (length(invalid_rows) > 0) {
      stop(sprintf(
        "predict_surrogates: %d invalid candidate point(s) (indices: %s). Each candidate must be a named numeric vector covering: %s.",
        length(invalid_rows),
        paste(utils::head(invalid_rows, 5), collapse = ", "),
        paste(param_names, collapse = ", ")
      ), call. = FALSE)
    }

    # Same finite check as the matrix path: correctly-named NA/Inf coordinates
    # would otherwise reach the predictors, whose error fallback silently
    # replaces predictions with constant high-uncertainty scores.
    if (!all(is.finite(design_mat))) {
      bad <- which(!apply(is.finite(design_mat), 1, all))
      stop(sprintf(
        "predict_surrogates: non-finite candidate coordinates in row(s): %s",
        paste(utils::head(bad, 5), collapse = ", ")
      ), call. = FALSE)
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
