#' Multi-fidelity ablation study for Bayesian optimisation
#'
#' Runs the calibration routine under several fidelity policies to quantify
#' simulation budgets and operating characteristic performance.
#'
#' @inheritParams benchmark_methods
#' @param policies named list mapping policy name to a named numeric vector of
#'   fidelity levels. Names must be a subset of `low`, `med`, `high`. Partial
#'   policies are normalized onto the full low/med/high contract expected by
#'   [bo_calibrate()]: each missing level inherits the nearest provided level
#'   at or below it (else the smallest provided). So `c(low = 500)` simulates
#'   500 replications regardless of which fidelity the optimizer nominally
#'   selects, and `c(low = 200, high = 1000)` treats a `med` selection as
#'   another `low`-cost evaluation.
#' @param seeds integer vector of seeds used for each replicate.
#'
#' @return Object of class `BATON_multifidelity` containing `runs` (per-policy
#'   records) and `summary` (policy-level averages).
#' @importFrom dplyr .data
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

  # Base R lapply + bind_rows instead of purrr::imap_dfr/map_dfr (Task D9).
  # imap_dfr passed each policy's name as the second argument (or the integer
  # index for an unnamed list); the explicit policy_names vector preserves that.
  policy_names <- if (is.null(names(policies))) seq_along(policies) else names(policies)
  runs <- dplyr::bind_rows(lapply(seq_along(policies), function(i) {
    policy_name <- policy_names[[i]]
    fidelity_levels <- normalize_fidelity_policy(policies[[i]], policy_name)
    if (progress) {
      message(sprintf("Ablation policy '%s' with fidelity levels: %s",
                      policy_name,
                      paste(names(fidelity_levels), fidelity_levels, sep = "=", collapse = ", ")))
    }
    dplyr::bind_rows(lapply(seeds, function(seed) {
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
    }))
  }))

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

#' Normalize a partial fidelity policy onto the full low/med/high contract
#'
#' `bo_calibrate()` requires `fidelity_levels` named exactly low/med/high
#' (a misnamed level would otherwise surface as a confusing subscript error
#' deep in fidelity selection). Ablation policies legitimately want fewer
#' levels ("always 500 reps"), so each missing level inherits the nearest
#' provided level at or below it (else the smallest provided): the optimizer
#' can nominally select any fidelity, but never spends more than the policy
#' authorized beneath that level.
#'
#' @keywords internal
normalize_fidelity_policy <- function(levels, policy_name) {
  nm <- names(levels)
  if (is.null(nm) || any(!nzchar(nm)) || !all(nm %in% c("low", "med", "high"))) {
    stop(sprintf(
      paste0("Ablation policy '%s': fidelity levels must be named with a ",
             "subset of 'low', 'med', 'high'; got: %s."),
      policy_name,
      if (is.null(nm)) "<unnamed>" else paste(nm, collapse = ", ")
    ), call. = FALSE)
  }
  if (anyDuplicated(nm) > 0) {
    stop(sprintf("Ablation policy '%s': duplicate fidelity level names.",
                 policy_name), call. = FALSE)
  }
  pick <- function(prefs) {
    for (p in prefs) if (p %in% nm) return(as.numeric(levels[[p]]))
  }
  c(low = pick(c("low", "med", "high")),
    med = pick(c("med", "low", "high")),
    high = pick(c("high", "med", "low")))
}

#' Visualise multi-fidelity budget vs accuracy trade-offs
#' @param ablation object returned by [ablation_multifidelity()].
#' @return ggplot object.
#' @export
plot_multifidelity_tradeoff <- function(ablation) {
  require_suggests("ggplot2")
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
