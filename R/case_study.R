#' Summarise calibrated case-study fit
#'
#' @param fit object returned by [bo_calibrate()].
#'
#' @return list with `design` (parameter/value tibble) and
#'   `operating_characteristics` (metric summaries).
#' @importFrom rlang .data
#' @export
summarise_case_study <- function(fit) {
  stopifnot(inherits(fit, "BATON_fit"))
  best_theta <- fit$best_theta
  design_tbl <- tibble::tibble(
    parameter = names(best_theta),
    value = purrr::map_dbl(best_theta, as.numeric)
  )

  posterior_draws <- fit$diagnostics$posterior_draws
  if (!is.null(posterior_draws) && length(posterior_draws) > 0) {
    oc_tbl <- purrr::imap_dfr(posterior_draws, function(draws, metric) {
      tibble::tibble(
        metric = metric,
        estimate = mean(draws, na.rm = TRUE),
        lower = stats::quantile(draws, 0.025, na.rm = TRUE),
        upper = stats::quantile(draws, 0.975, na.rm = TRUE)
      )
    })
  } else {
    best_idx <- best_feasible_index(fit$history, fit$policies$objective %||% names(fit$surrogates)[1])
    metrics <- fit$history$metrics[[best_idx]]
    oc_tbl <- tibble::tibble(
      metric = names(metrics),
      estimate = as.numeric(metrics),
      lower = NA_real_,
      upper = NA_real_
    )
  }

  list(design = design_tbl, operating_characteristics = oc_tbl)
}

#' Case-study diagnostics including Sobol and gradient summaries
#'
#' @param fit object from [bo_calibrate()].
#' @param sobol_samples Monte Carlo samples for Sobol indices.
#' @param gradient_points number of gradient samples.
#' @param eps finite difference step size for gradient computation (default 1e-3).
#'
#' @return list with Sobol, gradient, and covariance diagnostics.
#' @export
case_study_diagnostics <- function(fit,
                                   sobol_samples = 5000,
                                   gradient_points = 200,
                                   eps = 1e-3) {
  stopifnot(inherits(fit, "BATON_fit"))
  bounds <- infer_bounds_from_history(fit$history)
  objective <- fit$policies$objective %||% names(fit$surrogates)[1]
  sobol_tbl <- sa_sobol(fit$surrogates, bounds, outcome = objective, n_mc = sobol_samples)
  gradients_tbl <- gradient_sampling(fit$surrogates, bounds, outcome = objective,
                                     gradient_points = gradient_points, eps = eps)
  cov_matrix <- cov_effects(fit$surrogates, bounds, outcome = objective, n_mc = gradient_points, eps = eps)
  list(
    sobol = sobol_tbl,
    gradients = gradients_tbl,
    covariance = cov_matrix
  )
}

#' Feasible frontier plot highlighting regulatory constraints
#' @param fit object returned by [bo_calibrate()].
#' @return ggplot object.
#' @export
plot_feasible_frontier <- function(fit) {
  history_long <- history_with_metrics(fit$history)
  ggplot2::ggplot(history_long, ggplot2::aes(x = .data$metrics_EN, y = .data$metrics_power, colour = .data$feasible)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(x = "Expected sample size (EN)", y = "Power", colour = "Feasible") +
    ggplot2::theme_minimal()
}

#' Trade-off surface plot across key operating characteristics
#' @param fit object returned by [bo_calibrate()].
#' @return ggplot object.
#' @export
plot_tradeoff_surfaces <- function(fit) {
  history_long <- history_with_metrics(fit$history)
  ggplot2::ggplot(history_long, ggplot2::aes(x = .data$metrics_type1, y = .data$metrics_power, colour = .data$objective)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_colour_viridis_c(option = "C") +
    ggplot2::labs(x = "Type I error", y = "Power", colour = "Objective") +
    ggplot2::theme_minimal()
}

#' Plot Sobol indices for the case study
#' @param sobol_tbl tibble from [case_study_diagnostics()].
#' @return ggplot object.
#' @export
plot_case_sobol <- function(sobol_tbl) {
  plot_sobol_indices(sobol_tbl)
}

#' Plot case-study gradient heatmap
#' @param gradients_tbl tibble from [case_study_diagnostics()].
#' @return ggplot.
#' @export
plot_case_gradient <- function(gradients_tbl) {
  plot_gradient_heatmap(gradients_tbl)
}

#' Plot covariance matrix of gradient effects
#' @param covariance covariance matrix from [case_study_diagnostics()].
#' @return ggplot object.
#' @export
plot_case_covariance <- function(covariance) {
  param_names <- rownames(covariance)
  cov_long <- 
    tibble::tibble(
      parameter_i = rep(param_names, times = length(param_names)),
      parameter_j = rep(param_names, each = length(param_names)),
      value = as.numeric(covariance)
    )
  ggplot2::ggplot(cov_long, ggplot2::aes(x = .data$parameter_i, y = .data$parameter_j, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0) +
    ggplot2::labs(x = "Parameter", y = "Parameter", fill = "Covariance") +
    ggplot2::theme_minimal()
}

# --- internal helpers -------------------------------------------------------

history_with_metrics <- function(history) {
  tidyr::unnest_wider(history, .data$metrics, names_sep = "_")
}

infer_bounds_from_history <- function(history) {
  if (is.null(history) || nrow(history) == 0L) {
    stop("Cannot infer bounds from empty history.", call. = FALSE)
  }
  theta_values <- purrr::map(history$theta, ~ purrr::map_dbl(.x, as.numeric))
  if (length(theta_values) == 0L) {
    stop("Cannot infer bounds from empty history.", call. = FALSE)
  }
  params <- names(theta_values[[1]])
  purrr::set_names(param_bounds(theta_values), params)
}

param_bounds <- function(theta_values) {
  params <- names(theta_values[[1]])
  purrr::map(params, function(param) {
    values <- purrr::map_dbl(theta_values, ~ .x[[param]])
    range(values, na.rm = TRUE)
  })
}
