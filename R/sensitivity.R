#' Global sensitivity via Saltelli Sobol indices
#'
#' @param surrogates list of surrogate models (from [fit_surrogates()]).
#' @param bounds parameter bounds used to sample the design space.
#' @param outcome operating characteristic of interest.
#' @param n_mc number of Monte Carlo samples.
#'
#' @return tibble with first-order Sobol indices per parameter.
#' @importFrom rlang .data
#' @export
sa_sobol <- function(surrogates,
                     bounds,
                     outcome = "EN",
                     n_mc = 1000) {
  if (!outcome %in% names(surrogates)) {
    stop(sprintf("Outcome '%s' not available in surrogates.", outcome), call. = FALSE)
  }
  d <- length(bounds)
  if (n_mc < 100) {
    stop("`n_mc` should be at least 100 for stable Sobol estimates.", call. = FALSE)
  }
  param_names <- names(bounds)

  A <- lhs::randomLHS(n_mc, d)
  B <- lhs::randomLHS(n_mc, d)
  colnames(A) <- param_names
  colnames(B) <- param_names

  A_list <- unit_matrix_to_list(A)
  B_list <- unit_matrix_to_list(B)

  # Create subset with the actual outcome name, not the literal string "outcome"
  surrogate_subset <- surrogates[outcome]  # Keep the original name
  pred_A <- predict_surrogates(surrogate_subset, A_list)[[outcome]]$mean
  pred_B <- predict_surrogates(surrogate_subset, B_list)[[outcome]]$mean

  if (is.null(pred_A) || is.null(pred_B) || length(pred_A) == 0 || length(pred_B) == 0) {
    warning("Surrogate predictions returned NULL or empty; Sobol indices cannot be computed.",
            call. = FALSE)
    return(tibble::tibble(parameter = param_names, S_first = NA_real_))
  }

  var_y <- stats::var(c(pred_A, pred_B))
  if (!is.finite(var_y) || var_y <= 1e-10) {
    warning("Surrogate has near-zero variance; Sobol indices may be unreliable.",
            call. = FALSE)
    var_y <- max(var_y, 1e-8)
  }

  indices <- purrr::map_dbl(seq_len(d), function(j) {
    C <- A
    C[, j] <- B[, j]
    C_list <- unit_matrix_to_list(C)
    pred_C <- predict_surrogates(surrogate_subset, C_list)[[outcome]]$mean
    if (is.null(pred_C) || length(pred_C) == 0) {
      return(NA_real_)
    }
    estimate <- mean(pred_B * (pred_C - pred_A)) / var_y
    max(0, min(1, estimate))
  })

  tibble::tibble(parameter = param_names, S_first = indices)
}

#' Local gradients of emulator mean
#'
#' @param surrogates list of surrogates from [fit_surrogates()].
#' @param theta design point (named list of parameter values).
#' @param bounds parameter bounds.
#' @param outcome metric to differentiate.
#' @param eps step size on the unit scale.
#' @param return_sd logical; if `TRUE` also return approximate standard
#'   deviations from the surrogate uncertainty.
#'
#' @return tibble of gradients by parameter.
#' @export
sa_gradients <- function(surrogates,
                         theta,
                         bounds,
                         outcome = "EN",
                         eps = 1e-4,
                         return_sd = TRUE) {
  if (!outcome %in% names(surrogates)) {
    stop(sprintf("Outcome '%s' not available in surrogates.", outcome), call. = FALSE)
  }
  theta_unit <- scale_to_unit(theta, bounds)
  param_names <- names(bounds)
  model_list <- list(outcome = surrogates[[outcome]])

  # OPTIMIZED: Use lapply + bind_rows instead of map_dfr
  grad_list <- lapply(param_names, function(param) {
    grad_info <- surrogate_gradient(model_list, theta_unit, bounds, param, eps)
    tibble::tibble(
      parameter = param,
      gradient = grad_info$gradient,
      sd = if (return_sd) grad_info$sd else NA_real_
    )
  })
  gradients <- dplyr::bind_rows(grad_list)

  gradients
}

#' Covariance of gradient effects across the design space
#'
#' Approximates the covariance matrix of local effects by sampling gradients at
#' random points and computing their empirical covariance.
#'
#' @param surrogates list of surrogates from [fit_surrogates()].
#' @param bounds parameter bounds.
#' @param outcome metric of interest.
#' @param n_mc number of Monte Carlo locations used to estimate covariance.
#' @param eps step size for finite-difference gradients.
#'
#' @return covariance matrix (parameters × parameters).
#' @export
cov_effects <- function(surrogates,
                        bounds,
                        outcome = "EN",
                        n_mc = 500,
                        eps = 1e-4) {
  if (!outcome %in% names(surrogates)) {
    stop(sprintf("Outcome '%s' not available in surrogates.", outcome), call. = FALSE)
  }
  param_names <- names(bounds)
  d <- length(bounds)

  points <- lhs::randomLHS(n_mc, d)
  colnames(points) <- param_names
  # OPTIMIZED: Pre-allocate matrix instead of do.call(rbind, ...)
  gradient_mat <- matrix(NA_real_, nrow = n_mc, ncol = d)
  for (i in seq_len(n_mc)) {
    theta_unit <- points[i, , drop = TRUE]
    theta_list <- scale_from_unit(as.list(theta_unit), bounds)
    gradient_mat[i, ] <- sa_gradients(surrogates, theta_list, bounds, outcome = outcome, eps = eps, return_sd = FALSE)$gradient
  }
  cov_mat <- stats::cov(gradient_mat)
  dimnames(cov_mat) <- list(param_names, param_names)
  cov_mat
}

#' Comprehensive sensitivity diagnostics
#'
#' @inheritParams bo_calibrate
#' @param kernel_options character vector of covariance kernels to compare.
#' @param acquisition_options character vector of acquisition strategies.
#' @param prior_specs optional named list describing alternative priors (stored
#'   for reporting purposes).
#' @param bo_args list of arguments forwarded to [bo_calibrate()].
#' @param sobol_samples number of Monte Carlo points for Sobol indices.
#' @param gradient_points number of gradient evaluation points.
#'
#' @return List containing Sobol indices, gradient draws, kernel comparison
#'   metrics, and acquisition summaries.
#' @export
sensitivity_diagnostics <- function(sim_fun,
                                    bounds,
                                    objective,
                                    constraints,
                                    kernel_options = c("matern5_2"),
                                    acquisition_options = c("eci"),
                                    prior_specs = list(),
                                    bo_args = list(),
                                    sobol_samples = 2000,
                                    gradient_points = 200,
                                    integer_params = NULL,
                                    progress = TRUE,
                                    ...) {
  primary_kernel <- kernel_options[[1]]
  primary_acquisition <- acquisition_options[[1]]

  args <- modify_list(list(
    sim_fun = sim_fun,
    bounds = bounds,
    objective = objective,
    constraints = constraints,
    covtype = convert_kernel(primary_kernel),
    acquisition = primary_acquisition,
    integer_params = integer_params,
    progress = progress
  ), bo_args)

  fit <- do.call(bo_calibrate, c(args, list(...)))

  sobol_tbl <- sa_sobol(fit$surrogates, bounds, outcome = objective, n_mc = sobol_samples)

  gradient_samples <- gradient_sampling(fit$surrogates, bounds, outcome = objective,
                                        gradient_points = gradient_points, eps = 1e-3)

  kernel_cmp <- kernel_comparison(fit$history, bounds, objective, constraints, kernel_options)

  acquisition_cmp <- acquisition_comparison(fit$surrogates, bounds, objective, constraints, acquisition_options)

  list(
    baseline_fit = fit,
    sobol = sobol_tbl,
    gradients = gradient_samples,
    kernel_comparison = kernel_cmp,
    acquisition_comparison = acquisition_cmp,
    prior_specs = prior_specs
  )
}

# --- plotting helpers -------------------------------------------------------

#' Plot Sobol indices
#' @param sobol_tbl tibble returned by [sa_sobol()].
#' @return ggplot object.
#' @export
plot_sobol_indices <- function(sobol_tbl) {
  ggplot2::ggplot(sobol_tbl, ggplot2::aes(x = .data$parameter, y = .data$S_first, fill = .data$parameter)) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Parameter", y = "First-order Sobol index") +
    ggplot2::theme_minimal()
}

#' Plot gradient heatmap
#' @param gradient_tbl tibble from [sensitivity_diagnostics()] (`$gradients`).
#' @return ggplot object.
#' @export
plot_gradient_heatmap <- function(gradient_tbl) {
  ggplot2::ggplot(gradient_tbl, ggplot2::aes(x = .data$parameter, y = .data$point_id, fill = .data$gradient)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0) +
    ggplot2::labs(x = "Parameter", y = "Sample index", fill = "Gradient") +
    ggplot2::theme_minimal()
}

#' Plot kernel comparison summary
#' @param kernel_tbl tibble from [sensitivity_diagnostics()] (`$kernel_comparison`).
#' @return ggplot object.
#' @export
plot_kernel_comparison <- function(kernel_tbl) {
  ggplot2::ggplot(kernel_tbl, ggplot2::aes(x = .data$kernel, y = .data$loglik, colour = .data$metric, group = .data$metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(x = "Kernel", y = "Log-likelihood", colour = "Metric") +
    ggplot2::theme_minimal()
}

#' Plot acquisition comparison summary
#' @param acquisition_tbl tibble from [sensitivity_diagnostics()] (`$acquisition_comparison`).
#' @return ggplot object.
#' @export
plot_acquisition_comparison <- function(acquisition_tbl) {
  ggplot2::ggplot(acquisition_tbl, ggplot2::aes(x = .data$acquisition, y = .data$mean_score, fill = .data$acquisition)) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_score - .data$sd_score,
                                        ymax = .data$mean_score + .data$sd_score), width = 0.2) +
    ggplot2::labs(x = "Acquisition", y = "Mean acquisition score") +
    ggplot2::theme_minimal()
}

# --- internal utilities -----------------------------------------------------

unit_matrix_to_list <- function(mat) {
  param_names <- colnames(mat)
  purrr::map(seq_len(nrow(mat)), function(i) {
    vec <- mat[i, , drop = TRUE]
    purrr::set_names(as.numeric(vec), param_names)
  })
}

surrogate_gradient <- function(model_list, theta_unit, bounds, param, eps) {
  up <- down <- theta_unit
  up[param] <- min(1, up[param] + eps)
  down[param] <- max(0, down[param] - eps)
  if (abs(up[param] - down[param]) < 1e-8) {
    return(list(gradient = 0, sd = 0))
  }
  preds <- predict_surrogates(model_list, list(up, down))[[1]]
  theta_up <- scale_from_unit(as.list(up), bounds)
  theta_down <- scale_from_unit(as.list(down), bounds)
  delta <- as.numeric(theta_up[[param]] - theta_down[[param]])
  grad <- (preds$mean[1] - preds$mean[2]) / delta
  sd_grad <- sqrt(sum(preds$sd^2)) / delta
  list(gradient = grad, sd = sd_grad)
}

gradient_sampling <- function(surrogates, bounds, outcome, gradient_points, eps) {
  d <- length(bounds)
  param_names <- names(bounds)
  samples <- lhs::randomLHS(gradient_points, d)
  colnames(samples) <- param_names
  # OPTIMIZED: Use lapply + bind_rows instead of imap_dfr
  sample_list <- lapply(seq_len(gradient_points), function(i) {
    theta_unit <- samples[i, , drop = TRUE]
    theta_list <- scale_from_unit(as.list(theta_unit), bounds)
    grads <- sa_gradients(surrogates, theta_list, bounds, outcome = outcome, eps = eps, return_sd = FALSE)
    grads |>
      dplyr::mutate(point_id = i)
  })
  sample_tbl <- dplyr::bind_rows(sample_list)

  sample_tbl
}

kernel_comparison <- function(history, bounds, objective, constraints, kernel_options) {
  constraint_tbl <- parse_constraints(constraints)
  # OPTIMIZED: Use lapply + bind_rows instead of map_dfr
  kernel_results <- lapply(kernel_options, function(kernel) {
    fits <- fit_surrogates(history, objective, constraint_tbl, covtype = convert_kernel(kernel))
    metric_names <- names(fits)
    inner_results <- lapply(seq_along(fits), function(i) {
      model <- fits[[i]]
      metric <- metric_names[i]
      # DiceKriging's km class exposes log-likelihood via the S4 @logLik slot,
      # not an exported logLik() function. Use methods::slot() for portability.
      loglik <- tryCatch(as.numeric(methods::slot(model, "logLik")),
                         error = function(e) NA_real_)
      tibble::tibble(kernel = kernel, metric = metric, loglik = loglik)
    })
    dplyr::bind_rows(inner_results)
  })
  dplyr::bind_rows(kernel_results)
}

acquisition_comparison <- function(surrogates, bounds, objective, constraints, acquisition_options) {
  constraint_tbl <- parse_constraints(constraints)
  best_feasible <- Inf
  # OPTIMIZED: Use lapply + bind_rows instead of map_dfr
  acq_results <- lapply(acquisition_options, function(option) {
    candidates <- lhs_candidate_pool(500, bounds)
    scores <- evaluate_acquisition(option, candidates, surrogates, constraint_tbl, objective, best_feasible)
    tibble::tibble(
      acquisition = option,
      mean_score = mean(scores, na.rm = TRUE),
      sd_score = stats::sd(scores, na.rm = TRUE)
    )
  })
  acquisition_tbl <- dplyr::bind_rows(acq_results)
  acquisition_tbl
}

convert_kernel <- function(kernel_name) {
  switch(kernel_name,
         matern52 = "matern5_2",
         matern32 = "matern3_2",
         sqexp = "gauss",
         kernel_name)
}
