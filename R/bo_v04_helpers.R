# Copyright (c) 2026. For not-for-profit research and educational use only; all
# other rights reserved. See the LICENSE file for full terms.

# BATON v0.4.0 helpers for cross-philosophy warmstart and multi-seed verification.
#
# Two functions:
#   build_initial_history_from_warmstart(): converts a user-friendly
#     warmstart_from spec (list of designs / paths to *_best_design.csv)
#     into the lower-level initial_history tibble that bo_calibrate already
#     understands.
#   run_multi_seed_verification(): re-evaluates a calibrated design at
#     multiple independent seeds at high fidelity and computes summary
#     statistics + strict-feasibility flag. Used as Stage 4 gate.
#
# Both helpers are package-internal (not exported) but documented for
# maintainers.

#' Build initial_history tibble from a warmstart_from specification.
#'
#' @param warmstart_from List whose elements are one of:
#'   (i) a one-row data frame with parameter columns matching `names(bounds)`,
#'   (ii) a path to a `*_best_design.csv` file (long format with columns
#'     `parameter` and `value`),
#'   (iii) a one-row tibble in the BATON best-design wide format
#'     (column `total_n` etc., with operating-characteristic columns also
#'     present),
#'   (iv) a list with element `theta` (named list of parameter values).
#' @param bounds Named list of length-2 numeric vectors (parameter bounds).
#' @param objective Name of the objective metric (character).
#' @param constraint_tbl Constraint table from `parse_constraints()`.
#' @param progress Logical; if TRUE, emit a message about how many designs
#'   were converted.
#' @return A tibble suitable for use as `initial_history` in `bo_calibrate()`.
#' @keywords internal
build_initial_history_from_warmstart <- function(warmstart_from,
                                                  bounds,
                                                  objective,
                                                  constraint_tbl,
                                                  progress = TRUE) {
  if (!is.list(warmstart_from) || length(warmstart_from) == 0) {
    stop("`warmstart_from` must be a non-empty list of design specifications.",
         call. = FALSE)
  }
  param_names <- names(bounds)
  metric_names <- unique(c(objective, constraint_tbl$metric))

  # Helper: parse a single warmstart spec into a one-row data frame with
  # parameter columns + objective column + constraint metric columns +
  # a "feasible" column.
  parse_one <- function(spec, idx) {
    # Case (ii): path to *_best_design.csv
    if (is.character(spec) && length(spec) == 1) {
      if (!file.exists(spec)) {
        stop(sprintf("warmstart_from[[%d]]: file not found: %s", idx, spec),
             call. = FALSE)
      }
      d <- utils::read.csv(spec, stringsAsFactors = FALSE)
      # best_design.csv has columns: parameter, value, type
      if (!all(c("parameter", "value") %in% names(d))) {
        stop(sprintf(
          "warmstart_from[[%d]] (%s): expected columns 'parameter' and 'value' in best_design.csv",
          idx, spec), call. = FALSE)
      }
      # Pivot to wide
      vals <- as.numeric(d$value)
      names(vals) <- d$parameter
      # Build a one-row data frame with parameter columns and metric columns
      row <- as.data.frame(as.list(vals))
      return(row)
    }
    # Case (iv): list with theta element
    if (is.list(spec) && !is.data.frame(spec) && "theta" %in% names(spec)) {
      theta <- spec$theta
      if (!is.list(theta) && !is.numeric(theta)) {
        stop(sprintf(
          "warmstart_from[[%d]]: 'theta' must be a named list or numeric vector",
          idx), call. = FALSE)
      }
      vals <- unlist(theta)
      row <- as.data.frame(as.list(vals))
      # Pull metrics if present at top level of spec
      for (m in metric_names) {
        if (!is.null(spec[[m]])) {
          row[[m]] <- as.numeric(spec[[m]])
        }
      }
      if (!is.null(spec[["objective"]])) {
        row[["objective"]] <- as.numeric(spec[["objective"]])
      }
      return(row)
    }
    # Case (i) / (iii): one-row data frame
    if (is.data.frame(spec)) {
      if (nrow(spec) != 1) {
        stop(sprintf(
          "warmstart_from[[%d]]: data frame must have exactly one row (got %d)",
          idx, nrow(spec)), call. = FALSE)
      }
      return(as.data.frame(spec, stringsAsFactors = FALSE))
    }
    stop(sprintf(
      "warmstart_from[[%d]]: unrecognized spec type (%s). Expected data frame, file path, or list with 'theta' element.",
      idx, class(spec)[1]), call. = FALSE)
  }

  rows <- lapply(seq_along(warmstart_from), function(i) {
    parse_one(warmstart_from[[i]], i)
  })

  # Bind into a single data frame, taking union of columns
  all_cols <- unique(unlist(lapply(rows, names)))
  rows_padded <- lapply(rows, function(r) {
    missing_cols <- setdiff(all_cols, names(r))
    for (mc in missing_cols) {
      r[[mc]] <- NA_real_
    }
    r[, all_cols, drop = FALSE]
  })
  ws_df <- do.call(rbind, rows_padded)

  # Validate: parameter columns must be present
  missing_params <- setdiff(param_names, names(ws_df))
  if (length(missing_params) > 0) {
    stop(sprintf(
      "warmstart_from is missing parameter columns: %s. Provided columns: %s",
      paste(missing_params, collapse = ", "),
      paste(names(ws_df), collapse = ", ")
    ), call. = FALSE)
  }

  # Build initial_history-shaped tibble. We use the format that
  # bo_calibrate's initial_history handler accepts: parameter columns +
  # objective + fidelity + feasible. The handler will derive theta, unit_x,
  # theta_id, metrics, variance from these columns.
  ih <- data.frame(stringsAsFactors = FALSE)
  for (p in param_names) {
    ih[[p]] <- as.numeric(ws_df[[p]])
  }
  # Objective: use objective column if present, otherwise the metric named
  # like the objective.
  if ("objective" %in% names(ws_df) && !all(is.na(ws_df$objective))) {
    ih$objective <- as.numeric(ws_df$objective)
  } else if (objective %in% names(ws_df)) {
    ih$objective <- as.numeric(ws_df[[objective]])
  } else {
    ih$objective <- NA_real_
  }
  # Constraint metric columns
  for (m in metric_names) {
    if (m %in% names(ws_df)) {
      ih[[m]] <- as.numeric(ws_df[[m]])
    }
  }
  # Fidelity: assume the warmstart designs were calibrated at high fidelity
  # (typical for cross-philosophy warmstart from a previous calibration's
  # best design verified at R = 10,000).
  ih$fidelity <- "high"
  # Feasible: derive from constraint_tbl if all constraint columns present
  feasible_vec <- rep(NA, nrow(ih))
  if (all(constraint_tbl$metric %in% names(ih))) {
    feasible_vec <- vapply(seq_len(nrow(ih)), function(i) {
      metrics_i <- as.list(ih[i, constraint_tbl$metric, drop = FALSE])
      tryCatch(is_feasible(metrics_i, constraint_tbl), error = function(e) NA)
    }, logical(1))
  }
  ih$feasible <- feasible_vec

  if (progress) {
    n_feasible <- sum(ih$feasible, na.rm = TRUE)
    message(sprintf(
      "warmstart_from: built initial_history with %d design(s) (%d feasible)",
      nrow(ih), n_feasible
    ))
  }
  ih
}

#' Run multi-seed verification of a calibrated design.
#'
#' Re-evaluates a single design (`theta`) at multiple independent seeds at
#' high fidelity, returning per-seed operating-characteristic estimates and
#' summary statistics plus a strict-feasibility flag.
#'
#' @param sim_fun Simulator function (same signature as `bo_calibrate`).
#' @param theta Named list / vector of parameter values defining the design.
#' @param seeds Integer vector of seeds to evaluate at.
#' @param n_rep Number of replications per seed.
#' @param objective Name of the objective metric (character).
#' @param constraint_tbl Constraint table from `parse_constraints()`.
#' @param progress Logical; emit per-seed progress messages.
#' @param ... Forwarded to `sim_fun`.
#' @return A list with elements:
#'   \itemize{
#'     \item \code{per_seed}: data frame with one row per seed, columns for
#'       seed, the objective, each constraint metric, feasible flag,
#'       runtime_s.
#'     \item \code{summary}: list with mean/sd of objective and each
#'       constraint metric across seeds, plus strict_feas (1 if all seeds
#'       pass, 0 otherwise) and strict_feas_count (integer, number of
#'       feasible seeds).
#'   }
#' @keywords internal
run_multi_seed_verification <- function(sim_fun,
                                         theta,
                                         seeds,
                                         n_rep,
                                         objective,
                                         constraint_tbl,
                                         progress = TRUE,
                                         ...) {
  metric_names <- unique(c(objective, constraint_tbl$metric))
  per_seed <- vector("list", length(seeds))

  for (i in seq_along(seeds)) {
    seed_i <- as.integer(seeds[i])
    t0 <- Sys.time()
    set.seed(seed_i)
    res <- tryCatch(
      sim_fun(theta = theta, fidelity = "high", n_rep = n_rep,
              seed = seed_i, ...),
      error = function(e) {
        warning(sprintf("Multi-seed eval at seed %d failed: %s",
                        seed_i, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )
    runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(res)) {
      row <- data.frame(seed = seed_i, runtime_s = runtime_s,
                        feasible = NA, stringsAsFactors = FALSE)
      for (m in metric_names) row[[m]] <- NA_real_
      per_seed[[i]] <- row
      next
    }
    metrics_i <- res$metrics
    feasible_i <- tryCatch(
      is_feasible(metrics_i, constraint_tbl),
      error = function(e) NA
    )
    row <- data.frame(seed = seed_i, runtime_s = runtime_s,
                      feasible = feasible_i, stringsAsFactors = FALSE)
    for (m in metric_names) {
      v <- metrics_i[[m]]
      row[[m]] <- if (is.null(v)) NA_real_ else as.numeric(v)[1]
    }
    per_seed[[i]] <- row
    if (progress) {
      message(sprintf(
        "  Stage 4 seed %d: %s = %.4f, feasible = %s (%.1fs)",
        seed_i, objective,
        row[[objective]],
        as.character(feasible_i),
        runtime_s
      ))
    }
  }
  per_seed_df <- do.call(rbind, per_seed)

  # Build summary
  summary_list <- list(
    n_seeds = length(seeds),
    n_rep = as.integer(n_rep),
    seeds = as.integer(seeds),
    objective_mean = mean(per_seed_df[[objective]], na.rm = TRUE),
    objective_sd = stats::sd(per_seed_df[[objective]], na.rm = TRUE),
    strict_feas_count = sum(per_seed_df$feasible == TRUE, na.rm = TRUE),
    strict_feas = if (all(per_seed_df$feasible == TRUE, na.rm = TRUE) &&
                     all(!is.na(per_seed_df$feasible))) 1L else 0L
  )
  for (m in metric_names) {
    summary_list[[paste0(m, "_mean")]] <- mean(per_seed_df[[m]], na.rm = TRUE)
    summary_list[[paste0(m, "_sd")]] <- stats::sd(per_seed_df[[m]], na.rm = TRUE)
  }
  list(per_seed = per_seed_df, summary = summary_list)
}
