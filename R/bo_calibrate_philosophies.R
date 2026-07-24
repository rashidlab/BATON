# BATON v0.4.0 batch wrapper: calibrate multiple design philosophies for a
# single scenario in dependency order, with cross-philosophy warmstart and
# multi-seed verification both enabled by default.

#' Calibrate multiple design philosophies in dependency order with
#' cross-philosophy warmstart and multi-seed verification.
#'
#' This is a batch-level wrapper around `bo_calibrate()` for the common
#' workflow of calibrating multiple design philosophies (Optimal, Minimax,
#' Fleming, Admissible, Alt-Optimal) for a single scenario. It runs
#' philosophies in a dependency-aware order (the "donor" philosophies are
#' calibrated first, then "recipient" philosophies are calibrated with
#' warmstart from the donors' calibrated designs as Stage 0 seeds), and
#' invokes Stage 4 multi-seed verification on every calibrated design.
#'
#' This pipeline addresses two failure modes documented in the BATON
#' bounds-widening campaign (see `inst/doc/bounds_widening_findings.md` if
#' included, or the JASA manuscript Web Appendix G.5):
#'
#' \enumerate{
#'   \item \strong{BO local minima in high-dimensional spaces} — donor
#'     philosophies' calibrated points seed the surrogate's initial belief,
#'     helping it escape local minima that an LHS-only Stage 0 might miss.
#'     Demonstrated by the resolution of the multi-arm Optimal vs Fleming
#'     Pareto inversion.
#'   \item \strong{Stored-vs-actual operating-characteristic gaps} —
#'     calibration-time low-fidelity validation tier values can land
#'     favorably while the seed-averaged high-fidelity multi-seed value
#'     is far from the regulatory boundary. Demonstrated by SYMM H1 with
#'     stored type-I < 0.10 but multi-seed type-I = 0.165.
#' }
#'
#' @param scenario_id Character; identifier for the scenario being
#'   calibrated. Used in messages and the returned manifest.
#' @param sim_fun Simulator function (same as `bo_calibrate()`).
#' @param bounds Named list of length-2 numeric vectors defining the search
#'   box. Same for all philosophies in the batch.
#' @param philosophies Named list, one element per philosophy, each
#'   containing `objective` and `constraints` fields (the philosophy-specific
#'   calibration target). Order does NOT need to be dependency order; the
#'   wrapper sorts internally. Each `objective` (and every constraint name)
#'   must be the name of a metric that `sim_fun` returns; there is no built-in
#'   weighted-loss objective, so compute the weighted loss (e.g.
#'   `weighted_EN_N`) inside the simulator and return it as a named metric.
#' @param dependency_order Optional character vector specifying which
#'   philosophies are donors (calibrated first, no warmstart) and which are
#'   recipients (calibrated second, warmstart from all preceding successful
#'   calibrations). Default: donors are
#'   `c("Fleming", "Minimax", "Admissible")` and recipients are
#'   `c("Optimal", "Alt-Optimal")` (and any others not in donors); within
#'   each group, the order of `philosophies` is preserved.
#' @param multi_seed_verify Logical; pass to `bo_calibrate()`. Default `TRUE`
#'   (more conservative than `bo_calibrate()`'s default; the whole point of
#'   this wrapper is the verification gate).
#' @param multi_seed_n Integer; number of seeds for Stage 4 verification.
#'   Default 5.
#' @param multi_seed_strict Logical; if TRUE, recipient philosophies are
#'   warmstarted only from donors that passed multi-seed verification. If a
#'   recipient's warmstart pool is empty (all donors failed multi-seed),
#'   the recipient is calibrated with LHS-only Stage 0 and a warning is
#'   emitted. Default TRUE.
#' @param ... Additional arguments forwarded to every `bo_calibrate()` call
#'   (e.g., `n_init`, `budget`, `seed`, `fidelity_levels`).
#' @return A list with class `BATON_batch_fit` and elements:
#'   \itemize{
#'     \item \code{scenario_id}: as supplied
#'     \item \code{fits}: named list of `BATON_fit` objects, one per philosophy.
#'       A philosophy whose calibration errors is not fatal to the batch: its
#'       entry is an errored fit with `status = "errored"`, `error_message`
#'       (the message; `error` is kept as a back-compat alias),
#'       `best_theta = NULL`, and an empty `history` tibble.
#'     \item \code{manifest}: data frame with columns
#'       `philosophy`, `dependency_role` ("donor" or "recipient"), `verdict`
#'       (`"MULTI_SEED_PASS"` / `"MULTI_SEED_FAIL"` / `"MULTI_SEED_WARN"` /
#'       `"ERRORED"` for a failed calibration / `NA_character_`), `n_max`,
#'       `objective_value`, `warmstart_donors` (semi-colon-separated names),
#'       runtime_s
#'   }
#' @export
bo_calibrate_philosophies <- function(scenario_id,
                                       sim_fun,
                                       bounds,
                                       philosophies,
                                       dependency_order = NULL,
                                       multi_seed_verify = TRUE,
                                       multi_seed_n = 5L,
                                       multi_seed_strict = TRUE,
                                       ...) {
  if (!is.character(scenario_id) || length(scenario_id) != 1) {
    stop("`scenario_id` must be a single character string.", call. = FALSE)
  }
  if (!is.list(philosophies) || length(philosophies) == 0 ||
      is.null(names(philosophies))) {
    stop("`philosophies` must be a non-empty named list.", call. = FALSE)
  }
  for (name in names(philosophies)) {
    p <- philosophies[[name]]
    if (!is.list(p) || is.null(p$objective) || is.null(p$constraints)) {
      stop(sprintf(
        "philosophies[['%s']] must be a list with `objective` and `constraints` fields.",
        name), call. = FALSE)
    }
  }

  # Determine dependency order. Default: Fleming/Minimax/Admissible as donors;
  # everything else as recipients. The dependency_order argument lets users
  # override this for non-standard philosophy sets.
  default_donors <- c("Fleming", "Minimax", "Admissible", "Admissible_EN",
                      "Admissible_N", "balanced", "minimax", "admissible",
                      "admissible_EN", "admissible_N", "fleming")
  if (is.null(dependency_order)) {
    donors <- intersect(names(philosophies), default_donors)
    recipients <- setdiff(names(philosophies), donors)
  } else {
    if (!all(dependency_order %in% c("donor", "recipient"))) {
      stop("`dependency_order` must be a character vector of 'donor'/'recipient' labels with names matching `philosophies`.",
           call. = FALSE)
    }
    if (is.null(names(dependency_order)) ||
        !setequal(names(dependency_order), names(philosophies))) {
      stop("`dependency_order` must be named with elements matching names(philosophies).",
           call. = FALSE)
    }
    donors <- names(dependency_order)[dependency_order == "donor"]
    recipients <- names(dependency_order)[dependency_order == "recipient"]
  }

  message(sprintf(
    "[bo_calibrate_philosophies] scenario %s: %d donors (%s), %d recipients (%s)",
    scenario_id, length(donors), paste(donors, collapse = ", "),
    length(recipients), paste(recipients, collapse = ", ")
  ))

  fits <- list()
  manifest_rows <- list()

  # Isolate each philosophy's calibration: a failure in one must not discard the
  # other (often multi-hour) completed fits. On error, record an errored fit and
  # continue; downstream (warmstart specs, multi-seed filter, manifest) already
  # treats a fit with NULL best_theta / verdict "ERRORED" gracefully.
  safe_bo_calibrate <- function(name, role, ...) {
    tryCatch(
      bo_calibrate(...),
      error = function(e) {
        warning(sprintf("[%s | %s pass] calibration FAILED: %s",
                        name, role, conditionMessage(e)), call. = FALSE)
        fit <- build_baton_fit(
          history = initialise_history(),
          best_theta = NULL,
          best_objective = NA_real_,
          surrogates = NULL,
          policies = NULL,
          diagnostics = NULL,
          multi_seed_summary = NULL,
          multi_seed_runs = NULL,
          verdict = "ERRORED",
          bounds = NULL,
          constraints = NULL,
          constraint_tbl = NULL,
          status = "errored",
          error_message = conditionMessage(e)
        )
        # Back-compat alias; TODO(v0.8): drop, nothing in-package reads $error
        fit$error <- fit$error_message
        fit
      }
    )
  }

  # Pass 1: calibrate donors (no warmstart)
  for (name in donors) {
    p <- philosophies[[name]]
    t0 <- Sys.time()
    message(sprintf("[%s | donor pass] starting calibration", name))
    fit <- safe_bo_calibrate(
      name, "donor",
      sim_fun = sim_fun,
      bounds = bounds,
      objective = p$objective,
      constraints = p$constraints,
      multi_seed_verify = multi_seed_verify,
      multi_seed_n = multi_seed_n,
      multi_seed_strict = multi_seed_strict,
      ...
    )
    runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fits[[name]] <- fit
    manifest_rows[[name]] <- data.frame(
      philosophy = name,
      dependency_role = "donor",
      verdict = if (is.null(fit$verdict)) NA_character_ else fit$verdict,
      n_max = if (!is.null(fit$best_theta$total_n)) {
        as.integer(fit$best_theta$total_n)
      } else NA_integer_,
      objective_value = fit$best_objective,
      warmstart_donors = "",
      runtime_s = runtime_s,
      stringsAsFactors = FALSE
    )
    message(sprintf("[%s | donor pass] finished in %.1fs (verdict: %s)",
                    name, runtime_s,
                    if (is.null(fit$verdict)) "NA" else fit$verdict))
  }

  # Pass 2: calibrate recipients with warmstart from successful donors.
  # Convert each donor's calibrated design into a warmstart spec (one-row
  # data frame with parameter columns + objective + constraint metrics).
  build_donor_spec <- function(donor_fit, recipient_objective,
                                recipient_constraint_tbl) {
    theta <- donor_fit$best_theta
    if (is.null(theta) || length(theta) == 0) return(NULL)
    spec <- as.data.frame(as.list(unlist(theta)))
    # Attach the OCs of best_theta's OWN evaluation. Match by theta_id (the same
    # hash used to build history) so the seeded (theta, metrics) pair is
    # consistent; only if that lookup fails do we fall back to the last
    # high-fidelity feasible row. (Previously the code always used the last
    # high-fid feasible row, which need not be best_theta's evaluation, seeding
    # the recipient surrogate with a mismatched design/OC pair.)
    h <- donor_fit$history
    if (!is.null(h) && nrow(h) > 0) {
      row_idx <- NA_integer_
      bound_names <- intersect(names(theta), names(bounds))
      if (length(bound_names) == length(bounds) && "theta_id" %in% names(h)) {
        best_id <- theta_to_id(scale_to_unit(theta[bound_names], bounds))
        matched <- which(h$theta_id == best_id)
        if (length(matched) > 0) {
          # best_theta may have been evaluated more than once (e.g. low then high
          # fidelity). Prefer the high-fidelity feasible evaluation, then any
          # feasible one, so we attach the OCs that actually characterize it -
          # not a chronologically-later low-fidelity or infeasible re-eval.
          hf_feas <- matched[which(h$fidelity[matched] == "high" & h$feasible[matched])]
          feas <- matched[which(h$feasible[matched])]
          row_idx <- if (length(hf_feas) > 0) {
            tail(hf_feas, 1)
          } else if (length(feas) > 0) {
            tail(feas, 1)
          } else {
            tail(matched, 1)
          }
        }
      }
      if (is.na(row_idx)) {
        hf_idx <- which(h$fidelity == "high" & h$feasible)
        if (length(hf_idx) > 0) row_idx <- tail(hf_idx, 1)
      }
      if (!is.na(row_idx)) {
        donor_metrics <- h$metrics[[row_idx]]
        for (m in names(donor_metrics)) {
          spec[[m]] <- as.numeric(donor_metrics[[m]])
        }
      }
    }
    # Override objective column with the recipient's objective interpretation
    # (so the surrogate sees the right objective value at this point).
    if (recipient_objective %in% names(spec)) {
      spec$objective <- spec[[recipient_objective]]
    }
    spec
  }

  for (name in recipients) {
    p <- philosophies[[name]]
    # Determine which donors to use as warmstart sources
    eligible_donors <- donors
    if (multi_seed_strict) {
      eligible_donors <- donors[vapply(donors, function(d) {
        v <- fits[[d]]$verdict
        !is.null(v) && !is.na(v) && v == "MULTI_SEED_PASS"
      }, logical(1))]
      if (length(eligible_donors) < length(donors)) {
        skipped <- setdiff(donors, eligible_donors)
        warning(sprintf(
          "[%s | recipient pass] skipping %d donor(s) that did not pass multi-seed verification: %s",
          name, length(skipped), paste(skipped, collapse = ", ")
        ), call. = FALSE)
      }
    }
    warmstart_specs <- lapply(eligible_donors, function(d) {
      build_donor_spec(fits[[d]], p$objective, NULL)
    })
    warmstart_specs <- warmstart_specs[!vapply(warmstart_specs, is.null, logical(1))]
    if (length(warmstart_specs) == 0) {
      warning(sprintf(
        "[%s | recipient pass] no eligible warmstart donors; calibrating with LHS-only Stage 0",
        name), call. = FALSE)
      warmstart_arg <- NULL
    } else {
      warmstart_arg <- warmstart_specs
    }

    t0 <- Sys.time()
    message(sprintf(
      "[%s | recipient pass] starting calibration with %d warmstart donor(s): %s",
      name, length(warmstart_specs),
      paste(eligible_donors[seq_along(warmstart_specs)], collapse = ", ")
    ))
    fit <- safe_bo_calibrate(
      name, "recipient",
      sim_fun = sim_fun,
      bounds = bounds,
      objective = p$objective,
      constraints = p$constraints,
      warmstart_from = warmstart_arg,
      multi_seed_verify = multi_seed_verify,
      multi_seed_n = multi_seed_n,
      multi_seed_strict = multi_seed_strict,
      ...
    )
    runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fits[[name]] <- fit
    manifest_rows[[name]] <- data.frame(
      philosophy = name,
      dependency_role = "recipient",
      verdict = if (is.null(fit$verdict)) NA_character_ else fit$verdict,
      n_max = if (!is.null(fit$best_theta$total_n)) {
        as.integer(fit$best_theta$total_n)
      } else NA_integer_,
      objective_value = fit$best_objective,
      warmstart_donors = paste(eligible_donors[seq_along(warmstart_specs)],
                                collapse = ";"),
      runtime_s = runtime_s,
      stringsAsFactors = FALSE
    )
    message(sprintf("[%s | recipient pass] finished in %.1fs (verdict: %s)",
                    name, runtime_s,
                    if (is.null(fit$verdict)) "NA" else fit$verdict))
  }

  manifest <- do.call(rbind, manifest_rows[c(donors, recipients)])
  rownames(manifest) <- NULL

  structure(
    list(
      scenario_id = scenario_id,
      fits = fits,
      manifest = manifest
    ),
    class = c("BATON_batch_fit", "list")
  )
}
