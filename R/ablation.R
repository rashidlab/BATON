#' Multi-fidelity ablation study for Bayesian optimisation
#'
#' Runs the calibration routine under several fidelity policies to quantify
#' simulation budgets and operating characteristic performance.
#'
#' @inheritParams benchmark_methods
#' @param policies named list mapping policy name to a named numeric vector of
#'   fidelity levels (passed to the `fidelity_levels` argument of
#'   [bo_calibrate()]).
#' @param seeds integer vector of seeds used for each replicate.
#'
#' @return Object of class `BATON_multifidelity` containing `runs` (per-policy
#'   records) and `summary` (policy-level averages).
#' @importFrom rlang .data
#' @export
ablation_multifidelity <- function(sim_fun,
                                   bounds,
                                   objective,
                                   constraints,
                                   policies = list(
                                     low_only = c(low = 500),
                                     med_only = c(med = 2000),
                                     full = c(low = 200, med = 1000, high = 10000)
                                   ),
                                   seeds = 1:20,
                                   bo_args = list(),
                                   integer_params = NULL,
                                   progress = TRUE,
                                   ...) {
  validate_bounds(bounds)
  constraint_tbl <- parse_constraints(constraints)

  runs <- purrr::imap_dfr(policies, function(fidelity_levels, policy_name) {
    if (progress) {
      message(sprintf("Ablation policy '%s' with fidelity levels: %s",
                      policy_name,
                      paste(names(fidelity_levels), fidelity_levels, sep = "=", collapse = ", ")))
    }
    purrr::map_dfr(seeds, function(seed) {
      args <- modify_list(list(
        sim_fun = sim_fun,
        bounds = bounds,
        objective = objective,
        constraints = constraints,
        fidelity_levels = fidelity_levels,
        seed = seed,
        integer_params = integer_params,
        progress = FALSE
      ), bo_args)
      fit <- do.call(bo_calibrate, c(args, list(...)))
      history <- fit$history
      best_idx <- best_feasible_index(history, objective)
      tibble::tibble(
        policy = policy_name,
        seed = seed,
        feasible = history$feasible[[best_idx]],
        objective = history$objective[[best_idx]],
        metrics = list(history$metrics[[best_idx]]),
        theta = list(history$theta[[best_idx]]),
        total_sim_calls = sum(history$n_rep, na.rm = TRUE)
      )
    })
  })

  summary_tbl <- runs |>
    dplyr::group_by(.data$policy) |>
    dplyr::summarise(
      runs = dplyr::n(),
      feasible_rate = mean(.data$feasible),
      objective_mean = mean(.data$objective, na.rm = TRUE),
      objective_sd = stats::sd(.data$objective, na.rm = TRUE),
      sim_calls_mean = mean(.data$total_sim_calls, na.rm = TRUE),
      sim_calls_sd = stats::sd(.data$total_sim_calls, na.rm = TRUE),
      .groups = "drop"
    )

  structure(
    list(
      runs = runs,
      summary = summary_tbl,
      call = match.call()
    ),
    class = "BATON_multifidelity"
  )
}

#' Visualise multi-fidelity budget vs accuracy trade-offs
#' @param ablation object returned by [ablation_multifidelity()].
#' @return ggplot object.
#' @export
plot_multifidelity_tradeoff <- function(ablation) {
  stopifnot(inherits(ablation, "BATON_multifidelity"))
  ggplot2::ggplot(ablation$summary, ggplot2::aes(x = .data$sim_calls_mean, y = .data$objective_mean, colour = .data$policy)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$sim_calls_mean - .data$sim_calls_sd,
                                         xmax = .data$sim_calls_mean + .data$sim_calls_sd), height = 0) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$objective_mean - .data$objective_sd,
                                        ymax = .data$objective_mean + .data$objective_sd), width = 0) +
    ggplot2::labs(x = "Simulation calls (mean +/- sd)", y = "Objective (mean +/- sd)", colour = "Policy") +
    ggplot2::theme_minimal()
}
