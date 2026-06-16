#' Estimate constraint reliability for calibrated designs
#'
#' Runs the requested calibration strategies across multiple seeds and
#' re-evaluates the resulting designs with a large Monte Carlo budget to
#' quantify constraint satisfaction probabilities.
#'
#' @inheritParams benchmark_methods
#' @param calibration_seeds integer vector of seeds used for each repeated
#'   calibration.
#' @param validation_reps number of simulator replications used when validating
#'   a calibrated design (passed to `sim_fun` via the `n_rep` argument).
#'
#' @return An object of class `BATON_reliability` with components `summary`
#'   (strategy-level feasibility rates) and `runs` (per-run records).
#' @importFrom rlang .data
#' @export
estimate_constraint_reliability <- function(sim_fun,
                                            bounds,
                                            objective,
                                            constraints,
                                            strategies = c("bo", "grid", "random", "heuristic"),
                                            calibration_seeds = 1:20,
                                            validation_reps = 50000,
                                            bo_args = list(),
                                            grid_args = list(),
                                            random_args = list(),
                                            heuristic_args = list(),
                                            simulators_per_eval = list(),
                                            integer_params = NULL,
                                            progress = TRUE,
                                            ...) {
  validate_bounds(bounds)
  constraint_tbl <- parse_constraints(constraints)
  available <- c("bo", "grid", "random", "heuristic")
  strategies <- intersect(available, strategies)
  if (length(strategies) == 0L) {
    stop("No valid strategies supplied.", call. = FALSE)
  }

  runs <- purrr::imap_dfr(strategies, function(strategy, idx) {
    if (progress) {
      message(sprintf("Reliability: running '%s' across %d seed(s).", strategy, length(calibration_seeds)))
    }
    purrr::map_dfr(calibration_seeds, function(seed) {
      calibration <- run_strategy(
        strategy = strategy,
        seed = seed,
        sim_fun = sim_fun,
        bounds = bounds,
        objective = objective,
        constraint_tbl = constraint_tbl,
        constraints = constraints,
        bo_args = bo_args,
        grid_args = grid_args,
        random_args = random_args,
        heuristic_args = heuristic_args,
        simulators_per_eval = simulators_per_eval,
        integer_params = integer_params,
        progress = progress,
        ...
      )

      validation <- validate_design(
        sim_fun = sim_fun,
        theta = calibration$best_theta,
        objective = objective,
        constraint_tbl = constraint_tbl,
        n_rep = validation_reps,
        seed = seed + 100000L,
        ...
      )

      tibble::tibble(
        strategy = strategy,
        run_id = seed,
        calibration_feasible = calibration$feasible,
        validation_feasible = validation$feasible,
        calibration_objective = calibration$objective,
        validation_objective = validation$metrics[[objective]] %||% NA_real_,
        calibration_metrics = list(calibration$best_metrics),
        validation_metrics = list(validation$metrics),
        theta = list(calibration$best_theta),
        total_sim_calls = calibration$total_sim_calls + validation$n_rep
      )
    })
  })

  summary_tbl <- runs |>
    dplyr::group_by(.data$strategy) |>
    dplyr::summarise(
      runs = dplyr::n(),
      reliability = mean(.data$validation_feasible),
      calibration_feasible_rate = mean(.data$calibration_feasible),
      objective_mean = mean(.data$validation_objective, na.rm = TRUE),
      objective_sd = stats::sd(.data$validation_objective, na.rm = TRUE),
      .groups = "drop"
    )

  structure(
    list(
      summary = summary_tbl,
      runs = runs,
      call = match.call()
    ),
    class = "BATON_reliability"
  )
}

#' Plot constraint reliability comparison
#' @param reliability object returned by [estimate_constraint_reliability()].
#' @return ggplot object visualising feasibility rates.
#' @export
plot_constraint_reliability <- function(reliability) {
  stopifnot(inherits(reliability, "BATON_reliability"))
  ggplot2::ggplot(reliability$summary, ggplot2::aes(x = .data$strategy, y = .data$reliability, fill = .data$strategy)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", .data$reliability * 100)), vjust = -0.4, size = 3) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Strategy", y = "Validation feasibility rate", fill = "Strategy") +
    ggplot2::theme_minimal()
}

# --- internal helper --------------------------------------------------------

validate_design <- function(sim_fun, theta, objective, constraint_tbl, n_rep, seed, ...) {
  res <- sim_fun(theta, fidelity = "high", n_rep = n_rep, seed = seed, ...)
  variance <- attr(res, "variance", exact = TRUE)
  metrics <- if (is.list(res) && !is.null(res$metrics)) res$metrics else res
  metrics <- metrics |> as.list() |> purrr::imap_dbl(function(value, name) as.numeric(value))
  feasible <- is_feasible(metrics, constraint_tbl)
  list(
    metrics = metrics,
    variance = variance,
    feasible = feasible,
    n_rep = n_rep
  )
}
