#' Summarise calibrated case-study fit
#'
#' @param fit object returned by [bo_calibrate()].
#'
#' @return list with `design` (parameter/value tibble) and
#'   `operating_characteristics` (metric summaries).
#' @importFrom dplyr .data
#' @export
summarise_case_study <- function(fit) {
  stopifnot(inherits(fit, "BATON_fit"))
  best_theta <- fit$best_theta
  design_tbl <- tibble::tibble(
    parameter = names(best_theta),
    value = vapply(best_theta, as.numeric, numeric(1))
  )

  posterior_draws <- fit$diagnostics$posterior_draws
  if (!is.null(posterior_draws) && length(posterior_draws) > 0) {
    # Base R lapply + bind_rows instead of purrr::imap_dfr (Task D9); imap
    # passed each element's name (or the integer index for an unnamed list).
    metric_names <- if (is.null(names(posterior_draws))) {
      seq_along(posterior_draws)
    } else {
      names(posterior_draws)
    }
    oc_tbl <- dplyr::bind_rows(lapply(seq_along(posterior_draws), function(i) {
      draws <- posterior_draws[[i]]
      tibble::tibble(
        metric = metric_names[[i]],
        estimate = mean(draws, na.rm = TRUE),
        lower = stats::quantile(draws, 0.025, na.rm = TRUE),
        upper = stats::quantile(draws, 0.975, na.rm = TRUE)
      )
    }))
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
  if (is.null(fit$surrogates) || length(fit$surrogates) == 0) {
    stop("No surrogates in this fit. Slim fits (slim = TRUE) drop them; re-run with slim = FALSE for sensitivity analysis.", call. = FALSE)
  }
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
  require_suggests("ggplot2")
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
  require_suggests("ggplot2")
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
  require_suggests("ggplot2")
  plot_sobol_indices(sobol_tbl)
}

#' Plot case-study gradient heatmap
#' @param gradients_tbl tibble from [case_study_diagnostics()].
#' @return ggplot.
#' @export
plot_case_gradient <- function(gradients_tbl) {
  require_suggests("ggplot2")
  plot_gradient_heatmap(gradients_tbl)
}

#' Plot covariance matrix of gradient effects
#' @param covariance covariance matrix from [case_study_diagnostics()].
#' @return ggplot object.
#' @export
plot_case_covariance <- function(covariance) {
  require_suggests("ggplot2")
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
  # Base R replacement for tidyr::unnest_wider(metrics, names_sep = "_")
  # (Task D9): one numeric column per metric name observed across rows, with
  # NA_real_ where a row lacks that metric; the metrics list-column is dropped.
  metric_names <- unique(unlist(lapply(history$metrics, names)))
  out <- history
  for (nm in metric_names) {
    out[[paste0("metrics_", nm)]] <- vapply(history$metrics, function(m) {
      if (is.null(m) || !(nm %in% names(m))) NA_real_ else as.numeric(m[[nm]])
    }, numeric(1))
  }
  out$metrics <- NULL
  out
}

infer_bounds_from_history <- function(history) {
  if (is.null(history) || nrow(history) == 0L) {
    stop("Cannot infer bounds from empty history.", call. = FALSE)
  }
  # Base R lapply/vapply instead of purrr::map/map_dbl (Task D9).
  theta_values <- lapply(history$theta, function(th) vapply(th, as.numeric, numeric(1)))
  if (length(theta_values) == 0L) {
    stop("Cannot infer bounds from empty history.", call. = FALSE)
  }
  params <- names(theta_values[[1]])
  stats::setNames(param_bounds(theta_values), params)
}

param_bounds <- function(theta_values) {
  params <- names(theta_values[[1]])
  lapply(params, function(param) {
    values <- vapply(theta_values, function(v) v[[param]], numeric(1))
    range(values, na.rm = TRUE)
  })
}
