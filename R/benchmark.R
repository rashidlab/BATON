#' Benchmark calibration strategies across repeated runs
#'
#' Executes the Bayesian optimisation routine alongside baseline strategies
#' (grid search, random search, heuristic tuning) under common simulation
#' settings. Returns an object containing run-level results for downstream
#' summarisation and plotting.
#'
#' @param sim_fun simulator callback passed to [bo_calibrate()].
#' @param bounds parameter bounds.
#' @param objective objective metric name.
#' @param constraints constraint specification (see [bo_calibrate()]).
#' @param strategies character vector indicating which strategies to run (subset
#'   of `c("bo", "grid", "random", "heuristic")`).
#' @param bo_args list of arguments forwarded to [bo_calibrate()] (e.g.,
#'   `n_init`, `q`, `budget`, `seeds`).
#' @param grid_args list controlling grid search (`resolution`, `seeds`).
#' @param random_args list controlling random search (`n_samples`, `seeds`).
#' @param heuristic_args named list describing heuristic templates (`seeds`).
#' @param simulators_per_eval optional named list providing fidelity replication
#'   counts per strategy (e.g., `list(bo = c(low = 200, med = 1000, high = 10000))`).
#' @param integer_params parameters constrained to integers.
#' @param progress logical; if `TRUE`, strategy-level progress messages are
#'   emitted.
#' @param ... additional arguments forwarded to `sim_fun`.
#'
#' @return An object of class `BATON_benchmark` with components `results`
#'   (tibble of run-level records) and `call`.
#' @importFrom rlang .data
#' @export
benchmark_methods <- function(sim_fun,
                              bounds,
                              objective,
                              constraints,
                              strategies = c("bo", "grid", "random", "heuristic"),
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
    seeds <- strategy_seeds(strategy, bo_args, grid_args, random_args, heuristic_args)
    if (progress) {
      message(sprintf("Running strategy '%s' with %d replicate(s).", strategy, length(seeds)))
    }
    purrr::map_dfr(seeds, function(seed) {
      result <- run_strategy(
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
      tibble::tibble(
        strategy = strategy,
        run_id = seed,
        feasible = result$feasible,
        objective = result$objective,
        best_theta = list(result$best_theta),
        best_metrics = list(result$best_metrics),
        total_sim_calls = result$total_sim_calls,
        history = list(result$history)
      )
    })
  })

  structure(
    list(
      results = runs,
      call = match.call()
    ),
    class = "BATON_benchmark"
  )
}

#' Summarise benchmark performance across strategies
#' @param benchmark object returned by [benchmark_methods()].
#' @return tibble with strategy-level summaries.
#' @export
summarise_benchmark <- function(benchmark) {
  results <- benchmark$results
  if (is.null(results) || nrow(results) == 0L) {
    return(tibble::tibble())
  }
  expanded <- results |>
    tidyr::unnest_wider(best_metrics, names_sep = "_") |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("best_metrics_"), ~ tidyr::replace_na(.x, NA_real_)))

  metric_cols <- grep("^best_metrics_", names(expanded), value = TRUE)

  expanded |>
    dplyr::group_by(.data$strategy) |>
    dplyr::summarise(
      runs = dplyr::n(),
      feasible_rate = mean(.data$feasible),
      objective_mean = mean(.data$objective, na.rm = TRUE),
      objective_sd = stats::sd(.data$objective, na.rm = TRUE),
      dplyr::across(
        dplyr::all_of(metric_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "oc_{.col}"
      ),
      sims_mean = mean(.data$total_sim_calls, na.rm = TRUE),
      sims_sd = stats::sd(.data$total_sim_calls, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Plot benchmark trajectories (objective vs evaluation)
#' @param benchmark object returned by [benchmark_methods()].
#' @return ggplot object.
#' @export
plot_benchmark_trajectory <- function(benchmark) {
  stopifnot(inherits(benchmark, "BATON_benchmark"))
  history_tbl <- purrr::pmap_dfr(
    list(benchmark$results$strategy, benchmark$results$run_id, benchmark$results$history),
    function(strategy, run_id, history) {
      if (nrow(history) == 0L) {
        return(NULL)
      }
      history |> dplyr::mutate(strategy = strategy, run_id = run_id)
    }
  )

  if (nrow(history_tbl) == 0L) {
    return(ggplot2::ggplot())
  }

  ggplot2::ggplot(history_tbl, ggplot2::aes(x = .data$eval_id, y = .data$objective, colour = .data$strategy, group = interaction(.data$strategy, .data$run_id))) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::labs(
      x = "Simulation evaluations",
      y = "Objective",
      colour = "Strategy"
    ) +
    ggplot2::theme_minimal()
}

#' Plot simulation efficiency comparison
#' @param benchmark object returned by [benchmark_methods()].
#' @return ggplot object.
#' @export
plot_benchmark_efficiency <- function(benchmark) {
  summary_tbl <- summarise_benchmark(benchmark)
  ggplot2::ggplot(summary_tbl, ggplot2::aes(x = .data$sims_mean, y = .data$objective_mean, colour = .data$strategy)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$sims_mean - .data$sims_sd,
                                         xmax = .data$sims_mean + .data$sims_sd), height = 0) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$objective_mean - .data$objective_sd,
                                        ymax = .data$objective_mean + .data$objective_sd), width = 0) +
    ggplot2::labs(
      x = "Simulation calls (mean +/- sd)",
      y = "Objective (mean +/- sd)",
      colour = "Strategy"
    ) +
    ggplot2::theme_minimal()
}

# --- internal helpers -------------------------------------------------------

strategy_seeds <- function(strategy, bo_args, grid_args, random_args, heuristic_args) {
  arg_map <- list(
    bo = bo_args,
    grid = grid_args,
    random = random_args,
    heuristic = heuristic_args
  )
  seeds <- arg_map[[strategy]]$seeds
  if (is.null(seeds)) {
    seeds <- 1L
  }
  seeds
}

run_strategy <- function(strategy,
                         seed,
                         sim_fun,
                         bounds,
                         objective,
                         constraint_tbl,
                         constraints,
                         bo_args,
                         grid_args,
                         random_args,
                         heuristic_args,
                         simulators_per_eval,
                         integer_params,
                         progress,
                         ...) {
  if (strategy == "bo") {
    fidelity <- simulators_per_eval$bo %||% c(low = 200, med = 1000, high = 10000)
    args <- modify_list(list(
      sim_fun = sim_fun,
      bounds = bounds,
      objective = objective,
      constraints = constraints,
      seed = seed,
      fidelity_levels = fidelity,
      progress = FALSE
    ), bo_args)
    fit <- do.call(bo_calibrate, c(args, list(...)))
    history <- fit$history
    best_idx <- best_feasible_index(history, objective)
    best_metrics <- history$metrics[[best_idx]]
    return(list(
      feasible = history$feasible[[best_idx]],
      objective = history$objective[[best_idx]],
      best_theta = history$theta[[best_idx]],
      best_metrics = best_metrics,
      total_sim_calls = sum(history$n_rep, na.rm = TRUE),
      history = history
    ))
  }

  if (strategy == "grid") {
    history <- initialise_history()
    resolution <- grid_args$resolution
    if (is.null(resolution)) {
      resolution <- 10L
    }
    theta_grid <- expand_grid_points(bounds, resolution)
    eval_counter <- 0L
    for (theta in theta_grid) {
      eval_counter <- eval_counter + 1L
      history <- record_evaluation(
        history = history,
        sim_fun = sim_fun,
        theta = coerce_theta_types(theta, integer_params),
        bounds = bounds,
        objective = objective,
        constraint_tbl = constraint_tbl,
        fidelity = "high",
        fidelity_levels = c(high = simulators_per_eval$grid %||% 10000),
        eval_id = eval_counter,
        iter = 0L,
        seed = seed + eval_counter,
        progress = FALSE,
        ...
      )
    }
    best_idx <- best_feasible_index(history, objective)
    return(list(
      feasible = history$feasible[[best_idx]],
      objective = history$objective[[best_idx]],
      best_theta = history$theta[[best_idx]],
      best_metrics = history$metrics[[best_idx]],
      total_sim_calls = sum(history$n_rep, na.rm = TRUE),
      history = history
    ))
  }

  if (strategy == "random") {
    history <- initialise_history()
    n_samples <- random_args$n_samples %||% 500
    set.seed(seed)
    samples <- lhs_design(n_samples, bounds, seed = seed)
    eval_counter <- 0L
    for (theta in samples) {
      eval_counter <- eval_counter + 1L
      history <- record_evaluation(
        history = history,
        sim_fun = sim_fun,
        theta = coerce_theta_types(theta, integer_params),
        bounds = bounds,
        objective = objective,
        constraint_tbl = constraint_tbl,
        fidelity = "high",
        fidelity_levels = c(high = simulators_per_eval$random %||% 10000),
        eval_id = eval_counter,
        iter = 0L,
        seed = seed + eval_counter,
        progress = FALSE,
        ...
      )
    }
    best_idx <- best_feasible_index(history, objective)
    return(list(
      feasible = history$feasible[[best_idx]],
      objective = history$objective[[best_idx]],
      best_theta = history$theta[[best_idx]],
      best_metrics = history$metrics[[best_idx]],
      total_sim_calls = sum(history$n_rep, na.rm = TRUE),
      history = history
    ))
  }

  if (strategy == "heuristic") {
    history <- initialise_history()
    templates <- heuristic_templates(bounds, heuristic_args$template %||% "default")
    eval_counter <- 0L
    for (theta in templates) {
      eval_counter <- eval_counter + 1L
      history <- record_evaluation(
        history = history,
        sim_fun = sim_fun,
        theta = coerce_theta_types(theta, integer_params),
        bounds = bounds,
        objective = objective,
        constraint_tbl = constraint_tbl,
        fidelity = "high",
        fidelity_levels = c(high = simulators_per_eval$heuristic %||% 10000),
        eval_id = eval_counter,
        iter = 0L,
        seed = seed + eval_counter,
        progress = FALSE,
        ...
      )
    }
    best_idx <- best_feasible_index(history, objective)
    return(list(
      feasible = history$feasible[[best_idx]],
      objective = history$objective[[best_idx]],
      best_theta = history$theta[[best_idx]],
      best_metrics = history$metrics[[best_idx]],
      total_sim_calls = sum(history$n_rep, na.rm = TRUE),
      history = history
    ))
  }

  stop(sprintf("Unknown strategy '%s'.", strategy), call. = FALSE)
}

expand_grid_points <- function(bounds, resolution) {
  if (is.null(resolution) || length(resolution) == 0) {
    stop("`resolution` must be provided for grid search.", call. = FALSE)
  }
  params <- names(bounds)
  grids <- purrr::map2(bounds, params, function(range, name) {
    res <- if (is.null(names(resolution))) resolution[[1L]] else resolution[[name]] %||% resolution[[1L]]
    if (!is.numeric(res) || res <= 0) {
      stop(sprintf("Grid resolution for '%s' must be positive.", name), call. = FALSE)
    }
    seq(range[[1]], range[[2]], length.out = res)
  })
  names(grids) <- params
  combos <- do.call(expand.grid, c(grids, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  purrr::map(seq_len(nrow(combos)), function(i) {
    purrr::set_names(as.list(combos[i, , drop = TRUE]), params)
  })
}

heuristic_templates <- function(bounds, template) {
  mid <- purrr::map(bounds, ~ mean(.x))
  templates <- list(mid)
  clip <- function(value, range) {
    min(max(value, range[[1]]), range[[2]])
  }
  if (template == "two_stage") {
    templates <- c(
      templates,
      list(
        purrr::map(bounds, ~ clip(.x[[2]] * 0.9, .x)),
        purrr::map(bounds, ~ clip(.x[[1]] * 1.05, .x))
      )
    )
  }
  purrr::map(templates, purrr::set_names, names(bounds))
}

modify_list <- function(base, overlay) {
  utils::modifyList(base, overlay, keep.null = TRUE)
}
