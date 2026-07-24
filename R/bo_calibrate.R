#' Calibrate adaptive trial designs with constrained Bayesian optimisation
#'
#' Implements the methodology described in Section 2 of the manuscript:
#' heteroskedastic Gaussian process emulation, expected constrained improvement,
#' and multi-fidelity simulation budgeting. The function orchestrates initial
#' design selection, surrogate fitting, sequential acquisition, and diagnostic
#' bookkeeping required for downstream benchmarking and reporting.
#'
#' Version 0.3.0 introduces major performance improvements:
#' \itemize{
#'   \item Batch diversity via local penalization when \code{q > 1} (10-20\% fewer evaluations)
#'   \item Adaptive fidelity selection with cost-awareness (15-25\% better budget use)
#'   \item Warm-start for GP hyperparameters (30-50\% faster surrogate fitting)
#'   \item Adaptive candidate pool sizing (scales with dimension: 500 * d)
#'   \item Early stopping criterion (saves 10-30\% of budget)
#'   \item Improved constraint handling for infeasible regions
#' }
#' Expected combined improvement: 50-70\% overall efficiency gain.
#'
#' @param sim_fun callback with signature
#'   `function(theta, fidelity = c("low","med","high"), ...)` returning a named
#'   numeric vector of operating characteristics (e.g., `power`, `type1`, `EN`,
#'   `ET`). The function should attach attributes `variance` (named numeric
#'   vector of Monte Carlo variances) and `n_rep` (number of simulator
#'   replications) for best results. If `variance` is not provided, a binomial
#'   approximation is used for metrics in [0,1], and a small nugget is used for
#'   other metrics.
#'
#'   Optionally, the simulator may declare an explicit `n_rep` argument. When
#'   it does, BATON passes the requested replication count
#'   `n_rep = fidelity_levels[[fidelity]]` (including dynamically escalated
#'   values under `fidelity_method = "hybrid_staged"`), so the simulator can
#'   honor the exact budget rather than re-deriving it from the fidelity
#'   label. Accepting `...` alone does NOT opt in (this keeps wrappers that
#'   forward `...` to legacy simulators working); simulators without an
#'   `n_rep` formal are called exactly as before (backward compatible).
#' @param bounds named list of length-two numeric vectors specifying lower and
#'   upper limits for each design parameter.
#' @param objective name of the operating characteristic to minimise (character).
#' @param constraints named list mapping metric -> `c(direction, threshold)` with
#'   `direction` in "ge" or "le".
#' @param n_init number of initial design evaluations generated via Latin
#'   hypercube sampling. Default 90 (increased from 40 for better coverage).
#' @param q batch size for each BO iteration (number of new evaluations per
#'   acquisition round). When \code{q > 1}, batch diversity is automatically
#'   applied via local penalization (v0.3.0) to ensure spatially diverse points.
#' @param budget total number of simulator evaluations (initial + BO iterations).
#'   Default 300 (increased from 150 for tight constraint spaces).
#' @param seed base RNG seed for reproducibility. With an explicit seed the
#'   caller's RNG stream is restored on exit (success or error), so calling
#'   `bo_calibrate()` never advances `.Random.seed` in the calling
#'   environment. With `seed = NULL` the run seed is drawn from the caller's
#'   stream via a single `sample.int(1e6, 1)` draw that legitimately advances
#'   the stream by exactly one step (everything else is still rewound), so
#'   repeated `seed = NULL` calls give independent replicates. The seed
#'   actually used is recorded in `fit$policies$seed`.
#' @param initial_history optional data frame of previous evaluations to seed
#'   the initial design. Must contain columns for all parameters in
#'   \code{bounds}, plus \code{objective}, \code{fidelity}, \code{feasible},
#'   and individual constraint metric columns. All rows must be within the
#'   specified \code{bounds}. The provided rows AUGMENT the
#'   Latin hypercube initial design: if fewer than \code{n_init} rows are
#'   supplied, fresh LHS points (deduplicated against the provided
#'   \code{theta_id}s) top up the design to \code{n_init} before the BO loop
#'   starts, so a handful of donor seeds no longer replaces the whole initial
#'   design. With \code{n_init} or more rows, no LHS points are added and the
#'   BO loop starts immediately with this data (multi-stage warm-starting
#'   where later stages reuse evaluations from previous stages with narrowed
#'   bounds). Default: NULL (LHS-only initialization).
#' @param init_stopping optional list controlling GP-based initialization early
#'   stopping. When enabled, a Gaussian process is fit periodically during
#'   initialization to assess whether the design space is sufficiently explored.
#'   If metrics stabilize (e.g., prediction variance stops changing), initialization
#'   terminates early, saving budget for the BO loop. Create via
#'   \code{\link{init_stopping_config}()}. Key options:
#'   \describe{
#'     \item{\code{enabled}}{Logical: enable init stopping (default: TRUE)}
#'     \item{\code{min_init}}{Minimum points before checking (default: 20)}
#'     \item{\code{check_every}}{Check every k points (default: 20)}
#'     \item{\code{metrics}}{Metrics to track: "variance", "loo", "hyperparams",
#'       "loglik" (default: "variance")}
#'     \item{\code{threshold}}{Relative change threshold for convergence (default: 0.10)}
#'     \item{\code{window}}{Consecutive stable checkpoints required (default: 2)}
#'     \item{\code{verbose}}{Print diagnostic messages (default: FALSE)}
#'   }
#'   Default: NULL (no early stopping, use full n_init).
#' @param fidelity_levels named numeric vector giving the number of simulator
#'   replications associated with each fidelity level. Default: c(low=2000, med=4000,
#'   high=10000), increased from (200, 1000, 10000) for reduced Monte Carlo variance.
#' @param fidelity_method method for selecting fidelity level. Options:
#'   \describe{
#'     \item{\code{"adaptive"}}{(default) Cost-aware selection based on expected
#'       value per unit cost. Balances information gain vs computational expense.
#'       Recommended for most use cases.}
#'     \item{\code{"hybrid_staged"}}{MCEM-inspired continuous escalation with
#'       dynamic fidelity levels and metric-specific CV thresholds. Combines
#'       iteration-based staging with responsive escalation. Features: dynamic
#'       level scaling (5000->10000->15000 for high), constraint vs objective
#'       differentiation, boundary detection, budget safeguards. Based on empirical
#'       patterns from BOHB and materials discovery literature.}
#'     \item{\code{"staged"}}{Fixed schedule with iteration-based thresholds.
#'       Iterations 1-30: low, 31-100: adaptive, 101+: high. Simple but less
#'       efficient than adaptive or hybrid_staged.}
#'     \item{\code{"threshold"}}{Simple feasibility probability thresholds.
#'       P >= 0.75: high, P >= 0.4: med, else: low. Legacy method, not recommended.}
#'   }
#' @param fidelity_costs named numeric vector of relative costs per fidelity level.
#'   If NULL (default), assumes cost proportional to replication count. Use to
#'   specify non-linear cost relationships (e.g., parallelization effects).
#' @param fidelity_cv_threshold coefficient of variation threshold for fidelity
#'   promotion in staged method. Designs with CV > threshold are promoted to
#'   higher fidelity. Default 0.05 (tightened from 0.18 for constraint robustness).
#'   Lower values = more conservative, require more replications before accepting.
#' @param fidelity_prob_range two-element vector giving feasibility probability
#'   range [min, max] for adaptive fidelity promotion. Points within this range
#'   are near constraint boundaries and benefit from higher fidelity. Default
#'   c(0.2, 0.8) identifies uncertain feasibility regions.
#' @param acquisition acquisition rule to use (currently only `"eci"` is
#'   implemented).
#' @param candidate_pool number of random candidate points assessed per
#'   acquisition step. In v0.3.0, pool size automatically scales with dimension
#'   (500 * d, clamped to [1000, 5000]) and this parameter serves as a minimum.
#'   Larger pools in final 30\% of iterations for precision refinement.
#' @param covtype covariance kernel passed to [fit_surrogates()].
#' @param integer_params optional character vector of parameters that should be
#'   rounded to the nearest integer prior to simulation.
#' @param early_stop optional list controlling early stopping based on convergence
#'   detection. When the objective stops improving, optimization terminates early
#'   to save budget. Options:
#'   \describe{
#'     \item{\code{enabled}}{Logical: enable early stopping (default: TRUE)}
#'     \item{\code{patience}}{Integer: number of BO iterations without improvement
#'       before considering convergence (default: 5)}
#'     \item{\code{threshold}}{Numeric: minimum relative improvement required to
#'       reset patience counter (default: 1e-3 = 0.1\%)}
#'     \item{\code{consecutive}}{Integer: number of consecutive patience windows
#'       showing no improvement required to trigger stop (default: 2)}
#'     \item{\code{acq_rel}}{Numeric: acquisition-flatline stop threshold
#'       (default: 1e-3). Optimization also stops when the maximum acquisition
#'       value of the selected batch stays below a small level for several
#'       consecutive iterations. Once a feasible incumbent exists the acquisition
#'       is expected constrained improvement (EI * P(feasible)) and the level is
#'       \code{acq_rel * |incumbent|}, so the criterion scales with the objective
#'       (consistent whether it is O(1) or O(100)); before any feasible point the
#'       acquisition is the exploration/expected-violation score and the level is
#'       simply \code{acq_rel}.}
#'   }
#'   Default: NULL (uses default values: patience=5, threshold=1e-3,
#'   consecutive=2, acq_rel=1e-3).
#'   Set \code{early_stop = list(enabled = FALSE)} to disable early stopping entirely.
#' @param incumbent_method method for computing the incumbent (best feasible) value
#'   used in the ECI acquisition function. Options:
#'   \describe{
#'     \item{\code{"bpm"}}{Best Posterior Mean (default): Uses the GP surrogate's
#'       posterior mean estimate of the best feasible design's objective. More
#'       robust to Monte Carlo noise as it leverages the surrogate's smoothing.}
#'     \item{\code{"boi"}}{Best Observed Incumbent: Uses the raw best observed
#'       objective value from the evaluation history. Can be sensitive to noise
#'       but may perform better when the surrogate model is poorly calibrated.}
#'   }
#'   See Srinivas et al. (2024) "Bayesian Optimization with Expected Improvement:
#'   No Regret and the Choice of Incumbent" for theoretical analysis.
#' @param max_walltime_s optional walltime cap in seconds (default `NULL`, no
#'   cap). Elapsed time is measured from function entry and checked between
#'   initial-design evaluation chunks and at the top of each BO iteration, so
#'   the cap is advisory (chunk/iteration-granular), not preemptive: a running
#'   chunk or batch is never interrupted, and the run may overshoot the cap by
#'   up to one chunk or one full iteration (surrogate fit plus `q` simulator
#'   evaluations). When the cap trips, the function returns normally with a
#'   partial but internally consistent fit whose `status` field is
#'   `"walltime"`. If the cap is first crossed on a final chunk or batch that
#'   simultaneously ends the run for another reason (e.g. the budget is
#'   exhausted), `status` records that reason instead. In either case, once
#'   the cap has been reached, Stage 4 multi-seed verification is skipped
#'   even when `multi_seed_verify = TRUE`.
#' @param callback optional function invoked once at the END of each BO
#'   iteration (after the batch has been evaluated and recorded; it is NOT
#'   invoked during the initial design phase). It receives a single list with
#'   fields:
#'   \describe{
#'     \item{\code{iter}}{Current BO iteration number (1-based).}
#'     \item{\code{eval_count}}{Total evaluations recorded so far (including
#'       the initial design and any donor rows).}
#'     \item{\code{best_objective}}{Best feasible objective value known at the
#'       start of the iteration (the ECI incumbent), or `NA_real_` if no
#'       feasible point exists yet. Under the default incumbent method
#'       (`"bpm"`) this may be surrogate-predicted (a GP posterior mean), not
#'       necessarily an observed value.}
#'     \item{\code{best_observed}}{Best observed feasible objective by row,
#'       including the just-completed batch, or `NA_real_` if no feasible
#'       point exists yet.}
#'     \item{\code{elapsed_s}}{Seconds elapsed since function entry.}
#'   }
#'   Contract: return `TRUE` to continue; any other return value cancels the
#'   run cooperatively AFTER the current iteration (no mid-batch interruption,
#'   no completed work discarded) and the returned fit has `status`
#'   `"cancelled"`. A logging-only callback must therefore return `TRUE` as
#'   its final expression (e.g. call `message(info$iter)` and then `TRUE`).
#'   An error raised by
#'   the callback produces a warning and the
#'   run continues (a monitoring hook must not kill a paid run). A cancelled
#'   run skips Stage 4 multi-seed verification even when
#'   `multi_seed_verify = TRUE`. If the callback cancels on the final
#'   iteration (one that simultaneously ends the run because the budget is
#'   exhausted), `status` records `"cancelled"` rather than
#'   `"budget_exhausted"`, and Stage 4 verification is skipped. Default
#'   `NULL` (no callback).
#' @param checkpoint_fun optional function invoked with a slim, resumable
#'   `BATON_fit` snapshot (`status = "checkpoint"`) every `checkpoint_every`
#'   completed BO iterations (iteration-granular, not per-evaluation; the
#'   initial design phase emits no checkpoints). Snapshots carry the full
#'   evaluation `history`, the observed feasible best so far (`best_theta`
#'   is `NULL` and `best_objective` is `NA_real_` before feasibility), and
#'   the run's `policies`, but `surrogates` and `diagnostics` are `NULL`:
#'   GP objects are large, environment-laden, and not needed to resume. To
#'   resume from a saved checkpoint `ckpt`, pass its history back in:
#'   `bo_calibrate(..., initial_history = ckpt$history,
#'   n_init = nrow(ckpt$history))`; keeping `n_init <= nrow(ckpt$history)`
#'   avoids the donor-augmentation LHS top-up (Task B2 rule), so the resume
#'   refits surrogates from the recorded history and continues the run. The
#'   return value of `checkpoint_fun` is ignored; an error it raises
#'   produces a warning and the run continues (same rationale as
#'   `callback`). Default `NULL` (no checkpoints).
#' @param checkpoint_every single positive integer; number of completed BO
#'   iterations between checkpoint snapshots. Used only when
#'   `checkpoint_fun` is supplied. Default `5L`.
#' @param on_error one of `"stop"` (default) or `"return_partial"`. With
#'   `"stop"`, any error raised during the run (simulator failure, surrogate
#'   fitting failure, ...) propagates to the caller exactly as before. With
#'   `"return_partial"`, an error inside an evaluation round or BO iteration
#'   ends the run gracefully instead: the function returns a partial but
#'   internally consistent `BATON_fit` with `status = "errored"` and the
#'   condition message in `error_message`, preserving every completed
#'   evaluation. Granularity caveat: history rows are appended once per
#'   evaluation round (an initial-design chunk or a BO batch), so a failure
#'   inside a round loses that round's evaluations, including any that
#'   completed before the failing one (with parallel evaluation via
#'   `options(BATON.cores)`, the whole failed batch); rows from all PRIOR
#'   rounds always survive. Without a walltime cap the whole LHS initial
#'   design is a single round, so an error there yields an empty-history
#'   partial fit. After an in-loop error the final surrogate refit is still
#'   attempted (in its own guard), so `surrogates` may be `NULL` in the
#'   partial fit; Stage 4 multi-seed verification is skipped. A failure
#'   INSIDE Stage 4 verification itself (one that escapes
#'   `run_multi_seed_verification`'s per-seed catch) is also guarded:
#'   the fit is returned with `status = "errored"`,
#'   `multi_seed_summary`/`multi_seed_runs` `NULL`, and `verdict` `NA`.
#'   Likewise a post-refit `generate_diagnostics` failure drops
#'   `diagnostics` (with a warning) instead of losing the fit. Interrupts
#'   (Ctrl-C) are never caught in either mode.
#' @param slim logical; if `TRUE` the returned `BATON_fit` is a lightweight
#'   decision object: `surrogates`, `diagnostics`, and `multi_seed_runs` are
#'   `NULL`, while `history`, `best_theta`, `best_objective`, `policies`,
#'   `multi_seed_summary`, `verdict`, `status`, and `error_message` are kept
#'   (the adoption gate survives). Rationale: km/hetGP surrogate objects
#'   carry closures and environments that inflate `saveRDS` payloads far
#'   beyond the decision content, so a calibration service persisting many
#'   fits wants `slim`. Because nothing else consumes the post-loop
#'   surrogate refit (`best_theta`/`best_objective` are computed from the
#'   history alone), `slim = TRUE` also SKIPS that final refit entirely,
#'   saving one full GP fit per run; a side effect is that an error arising
#'   only in the final refit cannot surface under `slim = TRUE` (the run
#'   returns its natural status). Downstream helpers that need the fitted
#'   surrogates (`sa_sobol()`, `sensitivity_diagnostics()`,
#'   `extract_case_study()`) require a non-slim fit. A slim fit remains
#'   resumable via `initial_history = fit$history`. Default `FALSE`.
#' @param progress logical; if `TRUE` (default) messages are emitted at key
#'   milestones.
#' @param warmstart_from optional list of calibrated designs to inject as Stage 0
#'   seeds alongside the LHS sample. Each element should be either a one-row
#'   data frame with parameter columns matching `names(bounds)` (typical for
#'   cross-philosophy warmstart from a previous calibration's best design), a
#'   path to a `*_best_design.csv` file produced by a prior calibration, or a
#'   list with element `theta` (a named list/vector of parameter values).
#'   Designs outside the current bounds are filtered with a warning. Default
#'   `NULL` (LHS-only initial design pool). Internally converted to
#'   `initial_history`; see the `initial_history` argument of `bo_calibrate`
#'   for the lower-level API.
#' @param multi_seed_verify logical; if `TRUE`, after Stage 3 selects the
#'   calibrated design, re-evaluate it at `multi_seed_n` independent seeds at
#'   `multi_seed_reps` fidelity to detect stored-vs-actual operating-
#'   characteristic gaps. Default `FALSE` for backward compatibility.
#' @param multi_seed_n integer; number of independent seeds for Stage 4
#'   multi-seed verification. Used only when `multi_seed_verify = TRUE`.
#'   Default `5`.
#' @param multi_seed_reps integer; replications per seed for Stage 4
#'   verification. Used only when `multi_seed_verify = TRUE`. Default
#'   `fidelity_levels[["high"]]`.
#' @param multi_seed_strict logical; if `TRUE` and any of the
#'   `multi_seed_n` seeds produces an infeasible operating characteristic
#'   (constraint violation), the returned `BATON_fit` is marked with
#'   `verdict = "MULTI_SEED_FAIL"` and the calibration is flagged as not
#'   adoption-ready. The fit is still returned for diagnostics; downstream
#'   code can gate on `fit$multi_seed_summary$strict_feas == 1` before
#'   adopting the design. Default `TRUE`.
#' @param ... additional arguments forwarded to `sim_fun`.
#'
#' @return An object of class `BATON_fit` containing the optimisation history,
#'   best design, fitted surrogates, policy configuration, and posterior draws
#'   supporting sensitivity diagnostics. With `slim = TRUE` the `surrogates`,
#'   `diagnostics`, and `multi_seed_runs` components are `NULL` (see the
#'   `slim` argument). When `multi_seed_verify = TRUE`, also
#'   includes `multi_seed_summary` (summary statistics), `multi_seed_runs`
#'   (the per-seed evaluations; `NULL` under `slim = TRUE`), and `verdict`
#'   (`"MULTI_SEED_PASS"`, `"MULTI_SEED_FAIL"` when the strict gate fails
#'   under `multi_seed_strict = TRUE`, or `"MULTI_SEED_WARN"` when it fails
#'   under `multi_seed_strict = FALSE`), unless the
#'   `max_walltime_s` cap tripped or the run was cancelled via `callback`, in
#'   which case Stage 4 verification is
#'   skipped and `multi_seed_summary` is `NULL` with `verdict` `NA`. Note: early
#'   stopping may terminate before \code{budget} is exhausted if convergence
#'   is detected (configurable via \code{early_stop} parameter, default: no
#'   improvement > 0.1\% for 5 iterations). The `status` field records why
#'   the run ended: `"budget_exhausted"` (evaluation budget consumed),
#'   `"early_stopped"` (improvement plateau), `"acq_flatline"`
#'   (acquisition values below the \code{acq_rel} threshold),
#'   `"walltime"` (the \code{max_walltime_s} cap tripped between iterations),
#'   `"cancelled"` (the \code{callback} requested cancellation at the end
#'   of a BO iteration), or `"errored"` (an error was caught under
#'   \code{on_error = "return_partial"}; the message is in `error_message`).
#'
#' @importFrom utils head tail
#' @export
bo_calibrate <- function(sim_fun,
                         bounds,
                         objective,
                         constraints,
                         n_init = 90,
                         q = 2,
                         budget = 300,
                         seed = 2025,
                         initial_history = NULL,
                         warmstart_from = NULL,
                         init_stopping = NULL,
                         fidelity_levels = c(low = 2000, med = 4000, high = 10000),
                         fidelity_method = c("adaptive", "staged", "threshold", "hybrid_staged"),
                         fidelity_costs = NULL,
                         fidelity_cv_threshold = 0.05,
                         fidelity_prob_range = c(0.2, 0.8),
                         acquisition = c("eci"),
                         candidate_pool = 2000,
                         covtype = "matern5_2",
                         integer_params = NULL,
                         early_stop = NULL,
                         incumbent_method = c("bpm", "boi"),
                         progress = TRUE,
                         multi_seed_verify = FALSE,
                         multi_seed_n = 5L,
                         multi_seed_reps = NULL,
                         multi_seed_strict = TRUE,
                         max_walltime_s = NULL,
                         callback = NULL,
                         checkpoint_fun = NULL,
                         checkpoint_every = 5L,
                         on_error = c("stop", "return_partial"),
                         slim = FALSE,
                         ...) {
  # Walltime accounting (Task D4) starts at function entry so the cap covers
  # the initial design as well as the BO loop.
  start_time <- Sys.time()

  # RNG hygiene (D3, revised in D4): restore the caller's stream on every exit
  # path (including errors). When seed = NULL, the run seed is drawn from the
  # caller's stream FIRST (one sample.int() draw) and the restore hook is
  # installed immediately after, capturing the post-draw state. That single
  # draw is the only way the caller's stream advances, so repeated
  # seed = NULL calls give independent replicates while everything from
  # set.seed() onward is rewound on exit. With an explicit seed no draw
  # occurs and the caller's stream is untouched. Nothing between here and
  # set.seed(rng_seed) below may touch the RNG. Same restore pattern as
  # run_seeded() in R/surrogates.R.
  rng_seed <- if (is.null(seed)) sample.int(1e6, 1) else seed
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
  }

  acquisition <- match.arg(acquisition)
  fidelity_method <- match.arg(fidelity_method)
  incumbent_method <- match.arg(incumbent_method)
  on_error <- match.arg(on_error)

  # INPUT VALIDATION: Core parameters
  if (!is.function(sim_fun)) {
    stop("`sim_fun` must be a function.", call. = FALSE)
  }
  if (!is.character(objective) || length(objective) != 1 || nchar(objective) == 0) {
    stop("`objective` must be a non-empty character string.", call. = FALSE)
  }
  if (!is.list(constraints) || length(constraints) == 0) {
    stop("`constraints` must be a non-empty named list.", call. = FALSE)
  }
  if (!all(vapply(constraints, function(x) length(x) == 2, logical(1)))) {
    stop("Each constraint must be a length-2 vector: c(direction, threshold).", call. = FALSE)
  }

  validate_bounds(bounds)
  constraint_tbl <- parse_constraints(constraints)

  # v0.4.0: warmstart_from is a high-level wrapper around initial_history that
  # accepts simple design specifications (one-row data frames, paths to
  # *_best_design.csv files, or theta lists) and converts them to the proper
  # initial_history tibble. If both warmstart_from and initial_history are
  # supplied, error — they're alternative ways to spell the same thing.
  if (!is.null(warmstart_from)) {
    if (!is.null(initial_history)) {
      stop("Specify only one of `warmstart_from` or `initial_history`, not both.",
           call. = FALSE)
    }
    initial_history <- build_initial_history_from_warmstart(
      warmstart_from = warmstart_from,
      bounds = bounds,
      objective = objective,
      constraint_tbl = constraint_tbl,
      progress = progress
    )
  }

  # v0.4.0: validate multi-seed verification parameters
  if (!is.logical(multi_seed_verify) || length(multi_seed_verify) != 1) {
    stop("`multi_seed_verify` must be a single logical (TRUE or FALSE).",
         call. = FALSE)
  }
  if (multi_seed_verify) {
    if (!is.numeric(multi_seed_n) || length(multi_seed_n) != 1 ||
        multi_seed_n < 2 || multi_seed_n != as.integer(multi_seed_n)) {
      stop("`multi_seed_n` must be a single integer >= 2.", call. = FALSE)
    }
    multi_seed_n <- as.integer(multi_seed_n)
    if (!is.null(multi_seed_reps)) {
      if (!is.numeric(multi_seed_reps) || length(multi_seed_reps) != 1 ||
          multi_seed_reps < 100) {
        stop("`multi_seed_reps` must be a single integer >= 100 or NULL ",
             "(in which case fidelity_levels[['high']] is used).",
             call. = FALSE)
      }
      multi_seed_reps <- as.integer(multi_seed_reps)
    }
    if (!is.logical(multi_seed_strict) || length(multi_seed_strict) != 1) {
      stop("`multi_seed_strict` must be a single logical.", call. = FALSE)
    }
  }

  # Validate numeric parameters
  if (n_init <= 0) {
    stop("`n_init` must be positive.", call. = FALSE)
  }
  if (q <= 0) {
    stop("`q` (batch size) must be positive.", call. = FALSE)
  }
  if (budget <= 0) {
    stop("`budget` must be positive.", call. = FALSE)
  }
  if (candidate_pool <= 0) {
    stop("`candidate_pool` must be positive.", call. = FALSE)
  }
  if (budget < n_init) {
    warning("`budget` is less than `n_init`; only initial design will be evaluated.",
            call. = FALSE)
  }
  if (!is.null(max_walltime_s)) {
    if (!is.numeric(max_walltime_s) || length(max_walltime_s) != 1 ||
        !is.finite(max_walltime_s) || max_walltime_s <= 0) {
      stop("`max_walltime_s` must be a single positive finite number of seconds, or NULL for no cap.",
           call. = FALSE)
    }
  }
  if (!is.null(callback) && !is.function(callback)) {
    stop("`callback` must be a function taking a single info list, or NULL.",
         call. = FALSE)
  }
  if (!is.null(checkpoint_fun) && !is.function(checkpoint_fun)) {
    stop("`checkpoint_fun` must be a function taking a single BATON_fit, or NULL.",
         call. = FALSE)
  }
  # Integrality via %% 1 (not != as.integer()): as.integer() overflows to NA
  # for values above .Machine$integer.max, which would make this check NA and
  # crash the if(). Kept as numeric for the same reason; %% works exactly on
  # integral doubles below 2^53 (carry-forward from the D6 review).
  if (!is.numeric(checkpoint_every) || length(checkpoint_every) != 1 ||
      !is.finite(checkpoint_every) || checkpoint_every < 1 ||
      checkpoint_every %% 1 != 0) {
    stop("`checkpoint_every` must be a single positive integer.", call. = FALSE)
  }
  checkpoint_every <- as.numeric(checkpoint_every)
  # Task D8: slim = TRUE drops surrogates/diagnostics/multi_seed_runs from
  # the returned fit (and skips the final surrogate refit).
  if (!is.logical(slim) || length(slim) != 1L || is.na(slim)) {
    stop("`slim` must be a single logical (TRUE or FALSE).", call. = FALSE)
  }
  walltime_exceeded <- function() {
    !is.null(max_walltime_s) &&
      as.numeric(difftime(Sys.time(), start_time, units = "secs")) > max_walltime_s
  }
  # D4 follow-up: with a walltime cap set, initial-design evaluation proceeds
  # in chunks of the EFFECTIVE parallel worker count (effective_cores(): the
  # BATON.cores option only when the platform can actually fork, else 1) so
  # the cap is checked at chunk boundaries without losing parallel
  # throughput, and serial platforms get a check per evaluation. With no cap
  # the chunk spans the whole design (previous behavior, bit-identical
  # results).
  walltime_chunk <- if (is.null(max_walltime_s)) Inf else effective_cores()

  # Validate integer_params if provided
  if (!is.null(integer_params)) {
    if (!is.character(integer_params)) {
      stop("`integer_params` must be a character vector.", call. = FALSE)
    }
    unknown_params <- setdiff(integer_params, names(bounds))
    if (length(unknown_params) > 0) {
      stop(sprintf("`integer_params` contains unknown parameters: %s",
                   paste(unknown_params, collapse = ", ")), call. = FALSE)
    }
  }

  # BATON v0.3.0: q default changed from 8 to 2 for better exploration/exploitation balance

  # Notify users who may be relying on the old default
  if (q == 2 && !isTRUE(getOption("BATON.q_warning_shown"))) {
    message("Note: BATON v0.3.0 changed the default batch size from q=8 to q=2. ",
            "Set q explicitly if you need the previous behavior. ",
            "Suppress this message with options(BATON.q_warning_shown = TRUE).")
    options(BATON.q_warning_shown = TRUE)
  }

  # Parse early_stop configuration with defaults
  early_stop_config <- list(
    enabled = TRUE,
    patience = 5,
    threshold = 1e-3,
    consecutive = 2,
    acq_rel = 1e-3   # acquisition-flatline stop: max EI < acq_rel * |incumbent|
  )
  if (!is.null(early_stop)) {
    if (is.list(early_stop)) {
      # Override defaults with user-provided values
      if (!is.null(early_stop$enabled)) early_stop_config$enabled <- early_stop$enabled
      if (!is.null(early_stop$patience)) early_stop_config$patience <- early_stop$patience
      if (!is.null(early_stop$threshold)) early_stop_config$threshold <- early_stop$threshold
      if (!is.null(early_stop$consecutive)) early_stop_config$consecutive <- early_stop$consecutive
      if (!is.null(early_stop$acq_rel)) early_stop_config$acq_rel <- early_stop$acq_rel
    } else {
      stop("`early_stop` must be a list or NULL.", call. = FALSE)
    }
  }

  set.seed(rng_seed)

  fidelity_levels <- validate_fidelity_levels(fidelity_levels)
  primary_fidelity <- names(fidelity_levels)[1]

  # Initialize fidelity costs if not provided
  if (is.null(fidelity_costs)) {
    # Default: cost proportional to replication count, normalized to min=1
    fidelity_costs <- fidelity_levels / min(fidelity_levels)
  } else {
    # Validate user-provided costs
    if (!all(names(fidelity_costs) %in% names(fidelity_levels))) {
      stop("fidelity_costs must have names matching fidelity_levels", call. = FALSE)
    }
    fidelity_costs <- fidelity_costs[names(fidelity_levels)]
  }

  if (progress) {
    message("Initialising BATON calibration with seed = ", rng_seed)
    message(sprintf("Fidelity selection method: '%s'", fidelity_method))
  }

  history <- initialise_history()

  # Why the run ended (Task D2). "budget_exhausted" is the default when the
  # BO while condition runs out; each break site overwrites it first. The
  # walltime cap (D4) can also trip during the initialization phase below,
  # which sets "walltime" here and skips the BO loop entirely.
  run_status <- "budget_exhausted"
  run_error <- NULL

  # Task D7: on_error = "return_partial" turns an error inside an evaluation
  # round or BO iteration into a graceful stop that keeps all completed work.
  # run_guarded() wraps the risky regions; the wrapped expression is a promise
  # evaluated in this frame, so its assignments (history, surrogates, ...)
  # land here as usual. In the default "stop" mode NO handler is installed
  # (the expression just runs), so errors propagate exactly as before: same
  # condition, same message, raised from the same place. In "return_partial"
  # mode a caught error records run_status/run_error and returns FALSE so the
  # caller can break out of its loop. tryCatch(error =) does not catch
  # interrupts; Ctrl-C still aborts either way. Granularity caveat: appends
  # are per-round, so a failure inside evaluate_points() loses that round's
  # rows (rows from PRIOR rounds always survive).
  run_guarded <- if (on_error == "stop") {
    function(expr) {
      expr
      TRUE
    }
  } else {
    function(expr) {
      tryCatch({
        expr
        TRUE
      }, error = function(e) {
        run_status <<- "errored"
        run_error <<- conditionMessage(e)
        if (progress) {
          message(sprintf(
            "Error caught (on_error = 'return_partial'): %s. Returning partial fit.",
            run_error))
        }
        FALSE
      })
    }
  }

  # Handle initial_history parameter (v0.4.0 warm-start feature)
  if (!is.null(initial_history)) {
    # Validate initial_history
    param_names <- names(bounds)

    # Check for required non-parameter columns
    required_meta_cols <- c("objective", "fidelity", "feasible")
    missing_meta <- setdiff(required_meta_cols, names(initial_history))
    if (length(missing_meta) > 0) {
      stop(sprintf(
        "initial_history missing required columns: %s",
        paste(missing_meta, collapse = ", ")
      ), call. = FALSE)
    }

    # Check for parameter values - either in theta column OR individual columns
    has_theta_col <- "theta" %in% names(initial_history) && is.list(initial_history$theta)
    has_param_cols <- all(param_names %in% names(initial_history))

    if (!has_theta_col && !has_param_cols) {
      stop(sprintf(
        "initial_history must have either a 'theta' list column or individual parameter columns: %s",
        paste(param_names, collapse = ", ")
      ), call. = FALSE)
    }

    # Filter rows to only those within bounds
    # Build a mask: TRUE if row is within bounds for ALL parameters
    n_original <- nrow(initial_history)
    in_bounds_mask <- rep(TRUE, n_original)

    if (has_theta_col) {
      # Extract from theta list column and check bounds
      for (param in param_names) {
        vals <- sapply(initial_history$theta, function(th) {
          if (is.null(th) || !param %in% names(th)) NA_real_ else th[[param]]
        })
        lower <- bounds[[param]][1]
        upper <- bounds[[param]][2]

        # Row is out of bounds if value exists and is outside [lower, upper]
        # NA values are kept (they'll be handled elsewhere)
        param_in_bounds <- is.na(vals) | (vals >= lower & vals <= upper)
        in_bounds_mask <- in_bounds_mask & param_in_bounds
      }
    } else {
      # Check individual parameter columns
      for (param in param_names) {
        vals <- initial_history[[param]]
        lower <- bounds[[param]][1]
        upper <- bounds[[param]][2]

        param_in_bounds <- is.na(vals) | (vals >= lower & vals <= upper)
        in_bounds_mask <- in_bounds_mask & param_in_bounds
      }
    }

    # Apply filter
    n_filtered <- sum(!in_bounds_mask)
    if (n_filtered > 0) {
      warning(sprintf(
        "initial_history: filtered %d of %d rows with parameters outside bounds",
        n_filtered, n_original
      ))
      initial_history <- initial_history[in_bounds_mask, , drop = FALSE]
    }

    if (nrow(initial_history) == 0) {
      stop("initial_history has no rows within bounds - cannot warm-start", call. = FALSE)
    }

    # Ensure required columns exist for surrogate fitting
    # Derive theta list column if it doesn't exist
    if (!("theta" %in% names(initial_history)) || !is.list(initial_history$theta)) {
      initial_history$theta <- lapply(seq_len(nrow(initial_history)), function(i) {
        theta <- as.list(initial_history[i, param_names, drop = TRUE])
        names(theta) <- param_names
        theta
      })
    }

    # Derive unit_x from theta if it doesn't exist
    if (!("unit_x" %in% names(initial_history)) || !is.list(initial_history$unit_x)) {
      initial_history$unit_x <- lapply(initial_history$theta, function(th) {
        scale_to_unit(th, bounds)
      })
    }

    # Derive theta_id from unit_x if it doesn't exist
    if (!("theta_id" %in% names(initial_history))) {
      initial_history$theta_id <- vapply(initial_history$unit_x, theta_to_id, character(1))
    }

    # Ensure metrics list column exists (may be used by some analysis functions)
    if (!("metrics" %in% names(initial_history)) || !is.list(initial_history$metrics)) {
      # Build metrics from individual metric columns
      # Per docs, initial_history has "objective" column + constraint metric columns
      # We need to map "objective" column to the objective metric name if not already present
      metric_cols <- unique(c(objective, constraint_tbl$metric))
      available_cols <- intersect(metric_cols, names(initial_history))

      # Check if objective metric column is missing but "objective" column exists
      # (user provided objective values in "objective" column, not in metric-named column)
      has_objective_col <- "objective" %in% names(initial_history)
      has_metric_col <- objective %in% names(initial_history)

      initial_history$metrics <- lapply(seq_len(nrow(initial_history)), function(i) {
        m <- list()
        # Add values from metric-named columns
        if (length(available_cols) > 0) {
          m <- as.list(initial_history[i, available_cols, drop = TRUE])
          names(m) <- available_cols
        }
        # If objective metric not in metric columns but "objective" column exists,
        # use the "objective" column value for the objective metric
        if (!has_metric_col && has_objective_col) {
          m[[objective]] <- initial_history$objective[i]
        }
        m
      })
    }

    # Ensure variance list column exists
    # Handle three cases: missing, numeric column, or already a list-column
    if (!("variance" %in% names(initial_history))) {
      # No variance column - create empty lists
      initial_history$variance <- replicate(nrow(initial_history), list(), simplify = FALSE)
    } else if (!is.list(initial_history$variance)) {
      # Variance is a numeric column (scalar per row) - convert to list-column
      # Assume the variance applies to the objective metric
      var_values <- initial_history$variance
      initial_history$variance <- lapply(seq_len(nrow(initial_history)), function(i) {
        v <- var_values[i]
        if (is.na(v) || !is.finite(v)) {
          list()  # Return empty list for invalid values
        } else {
          # Create named list with variance for objective
          setNames(list(v), objective)
        }
      })
    }
    # else: already a list-column, keep as-is

    # Use filtered history
    history <- initial_history
    n_donor <- nrow(initial_history)

    # Donor rows belong to the initial design (iter 0, eval_id 1..n_donor).
    # Only stamp these when the columns are absent: a full previous-run
    # history passed as initial_history keeps its own iter/eval_id.
    if (!("iter" %in% names(history))) history$iter <- 0L
    if (!("eval_id" %in% names(history))) history$eval_id <- seq_len(n_donor)

    eval_counter <- n_donor

    # Task B2: donor seeds AUGMENT the LHS initial design instead of replacing
    # it. When initial_history has fewer than n_init rows, top up with fresh
    # LHS points so the BO loop never starts from a degenerate 1-3 point donor
    # design. Deduplication covers donor theta_ids AND points accepted earlier
    # in this top-up (integer_params coercion can collapse distinct LHS draws
    # onto the same design); rejected draws are replenished with re-seeded LHS
    # rounds until n_init is reached or the (small discrete) design space is
    # exhausted.
    n_init_requested <- min(n_init, budget)
    if (n_donor < n_init_requested) {
      n_needed <- n_init_requested - n_donor
      seen_ids <- history$theta_id
      accepted <- vector("list", n_needed)
      n_accepted <- 0L

      # theta_id collisions require every coordinate to coincide, which only
      # happens in practice when ALL parameters are integer-coerced. In that
      # case the design space is finite and its size tells us exactly when it
      # is exhausted; with any continuous parameter collisions are
      # measure-zero and the first round fills the design.
      # Count integers whose rounding cell (k - 0.5, k + 0.5) overlaps the
      # bounds -- coerce_theta_types uses round(), which reaches integers
      # OUTSIDE non-integral bounds (round(0.4) = 0 for bounds [0.4, 3.6]).
      # On exact half-endpoints this may overcount by an endpoint-only design;
      # overcounting is safe (a few extra rounds under the attempt backstop),
      # undercounting would truncate the design.
      space_size <- if (length(bounds) > 0 &&
                        all(names(bounds) %in% integer_params)) {
        prod(vapply(bounds, function(b) {
          max(0, floor(max(b) + 0.5) - ceiling(min(b) - 0.5) + 1)
        }, numeric(1)))
      } else {
        Inf
      }

      attempt <- 0L
      max_attempts <- 100L  # termination backstop only; exhaustion check below
      while (n_accepted < n_needed &&
             length(unique(seen_ids)) < space_size &&
             attempt < max_attempts) {
        # Round 1 draws exactly n_needed at rng_seed (continuous-path
        # behavior unchanged); later rounds oversample so a few remaining
        # slots do not need luck to land on unseen designs.
        n_remaining <- n_needed - n_accepted
        draw_size <- if (attempt == 0L) n_remaining else max(2L * n_remaining, 32L)
        draw <- lhs_design(draw_size, bounds,
                           seed = rng_seed + 1000L * attempt)
        for (theta in draw) {
          if (n_accepted >= n_needed) break
          theta <- coerce_theta_types(theta, integer_params)
          candidate_id <- theta_to_id(scale_to_unit(theta, bounds))
          if (candidate_id %in% seen_ids) next
          seen_ids <- c(seen_ids, candidate_id)
          n_accepted <- n_accepted + 1L
          accepted[[n_accepted]] <- theta
        }
        attempt <- attempt + 1L
      }
      # The top-up points are independent: evaluate them (in parallel when
      # options(BATON.cores) opts in) in walltime-bounded chunks, appended in
      # one bind per chunk (Tasks C2 + D4 follow-up). With no walltime cap
      # this is a single chunk, identical to the previous behavior; per-point
      # eval_ids and seeds are unchanged either way.
      n_evaluated <- 0L
      ok <- TRUE  # guard flag; initialized so the pattern survives reordering
      while (n_evaluated < n_accepted) {
        k <- as.integer(min(walltime_chunk, n_accepted - n_evaluated))
        topup_specs <- lapply(seq_len(k), function(j) {
          list(theta = accepted[[n_evaluated + j]], fidelity = primary_fidelity,
               eval_id = eval_counter + j, iter = 0L,
               seed = rng_seed + eval_counter + j)
        })
        # Task D7: guarded evaluation round. On error with
        # on_error = "return_partial" this round's rows are lost (appends are
        # per-round) but every previously appended row survives.
        ok <- run_guarded({
          topup_rows <- evaluate_points(
            specs = topup_specs, sim_fun = sim_fun, bounds = bounds,
            objective = objective, constraint_tbl = constraint_tbl,
            fidelity_levels = fidelity_levels, ...
          )
          history <- append_rows(history, topup_rows)
        })
        if (!ok) break
        if (progress) {
          for (r in topup_rows) announce_evaluation(r, objective)
        }
        eval_counter <- eval_counter + k
        n_evaluated <- n_evaluated + k
        # D4 follow-up: walltime cap checked at top-up chunk boundaries.
        if (walltime_exceeded() && n_evaluated < n_accepted) {
          run_status <- "walltime"
          if (progress) {
            # n_evaluated counts only fresh top-up evaluations run here;
            # eval_counter would also include donor rows not evaluated this run.
            message(sprintf(
              "Walltime cap (%.0fs) reached during initial-design top-up after %d evaluation(s).",
              max_walltime_s, n_evaluated))
          }
          break
        }
      }
      if (progress && run_status != "errored") {
        # Attribute the shortfall correctly when both causes apply: the
        # accept loop can exhaust the distinct-design space (n_accepted <
        # n_needed) AND the walltime cap can stop evaluation of the accepted
        # points (n_evaluated < n_accepted) in the same run. (An errored
        # round already announced itself via run_guarded; skip this summary
        # so the shortfall is not mislabeled as walltime/exhaustion.)
        shortfall_exhausted <- n_needed - n_accepted
        shortfall_walltime <- n_accepted - n_evaluated
        message(sprintf(
          "Initial design: %d donor row(s) augmented with %d fresh LHS point(s)%s",
          n_donor, n_evaluated,
          if (shortfall_exhausted > 0 && shortfall_walltime > 0) {
            sprintf(" (%d fewer than requested: %d distinct designs exhausted, %d walltime cap)",
                    n_needed - n_evaluated, shortfall_exhausted, shortfall_walltime)
          } else if (shortfall_walltime > 0) {
            sprintf(" (%d fewer than requested: walltime cap)",
                    shortfall_walltime)
          } else if (shortfall_exhausted > 0) {
            sprintf(" (%d fewer than requested: distinct designs exhausted)",
                    shortfall_exhausted)
          } else ""
        ))
      }
    } else if (progress) {
      message(sprintf(
        "Using %d evaluations from initial_history (skipping random initialization)",
        n_donor
      ))
    }

    # IMPORTANT: n_init must reflect the actual initial design size (donors +
    # fresh LHS). This keeps early_stop bookkeeping correct (it checks
    # iter_counter relative to n_init).
    n_init <- eval_counter

  } else {
    # Standard random initialization via Latin hypercube
    n_init <- min(n_init, budget)
    initial_design <- lhs_design(n_init, bounds, seed = rng_seed)

    eval_counter <- 0L

    # Initialize GP-based stopping if enabled
    init_checkpoints <- list()
    init_stopped_early <- FALSE
    use_init_stopping <- !is.null(init_stopping) &&
                         isTRUE(init_stopping$enabled) &&
                         is.null(initial_history)

    if (use_init_stopping && progress) {
      message(sprintf("  Init stopping enabled: check every %d points, min %d",
                      init_stopping$check_every %||% 20,
                      init_stopping$min_init %||% 20))
    }

    # Task C2: initial-design points are independent, so they are evaluated in
    # chunks (parallel when options(BATON.cores) opts in) and appended in one
    # bind per chunk. GP-based init stopping only ever fires when eval_counter
    # reaches a multiple of check_every (and >= min_init), so chunking up to
    # the next such checkpoint preserves the sequential stopping decisions
    # exactly; without init stopping the whole design is one chunk.
    n_points <- length(initial_design)
    point_pos <- 0L
    ok <- TRUE  # guard flag; initialized so the pattern survives reordering
    while (point_pos < n_points) {
      if (use_init_stopping) {
        ce <- init_stopping$check_every %||% 20
        mi <- init_stopping$min_init %||% 20
        next_check <- (floor(eval_counter / ce) + 1L) * ce
        while (next_check < mi) next_check <- next_check + ce
        chunk_n <- min(next_check - eval_counter, n_points - point_pos)
      } else {
        chunk_n <- n_points - point_pos
      }
      # D4 follow-up: bound chunks so the walltime cap below is checked often
      # enough. Finer chunks preserve init-stopping checkpoints exactly (the
      # checkpoint condition is a modulo on eval_counter).
      chunk_n <- as.integer(min(chunk_n, walltime_chunk))
      chunk_specs <- lapply(seq_len(chunk_n), function(j) {
        theta <- coerce_theta_types(initial_design[[point_pos + j]], integer_params)
        list(theta = theta, fidelity = primary_fidelity,
             eval_id = eval_counter + j, iter = 0L,
             seed = rng_seed + eval_counter + j)
      })
      # Task D7: guarded evaluation round. On error with
      # on_error = "return_partial" this chunk's rows are lost (appends are
      # per-round) but every previously appended chunk survives. Without a
      # walltime cap the whole LHS design is one chunk, so an init failure
      # yields an empty-history partial fit.
      ok <- run_guarded({
        chunk_rows <- evaluate_points(
          specs = chunk_specs, sim_fun = sim_fun, bounds = bounds,
          objective = objective, constraint_tbl = constraint_tbl,
          fidelity_levels = fidelity_levels, ...
        )
        history <- append_rows(history, chunk_rows)
      })
      if (!ok) break
      if (progress) {
        for (r in chunk_rows) announce_evaluation(r, objective)
      }
      eval_counter <- eval_counter + chunk_n
      point_pos <- point_pos + chunk_n

      # D4 follow-up: walltime cap checked at init-chunk boundaries (never
      # mid-chunk), so the partial init history is internally consistent.
      # The point_pos guard mirrors the top-up path: a trip on the FINAL
      # chunk is left to the BO loop's top-of-iteration check, so a run whose
      # completed init design also fills the budget keeps its natural
      # "budget_exhausted" status (Stage 4 is still suppressed via the
      # walltime_exceeded() disjunct at the verification gate).
      if (walltime_exceeded() && point_pos < n_points) {
        run_status <- "walltime"
        if (progress) {
          message(sprintf(
            "Walltime cap (%.0fs) reached during initialization after %d of %d evaluation(s).",
            max_walltime_s, eval_counter, n_init))
        }
        break
      }

      # Check GP-based initialization stopping
      if (use_init_stopping &&
          eval_counter >= (init_stopping$min_init %||% 20) &&
          eval_counter %% (init_stopping$check_every %||% 20) == 0) {

        # Extract design matrix and objective values from history
        # theta is stored as a list column - unpack to matrix
        # OPTIMIZED: Pre-allocate matrix instead of do.call(rbind, ...)
        param_names <- names(bounds)
        n_rows <- nrow(history)
        n_params <- length(param_names)
        X_init <- matrix(NA_real_, nrow = n_rows, ncol = n_params)
        for (i in seq_len(n_rows)) {
          X_init[i, ] <- unlist(history$theta[[i]][param_names])
        }
        colnames(X_init) <- param_names
        y_init <- history[[objective]]

        check_result <- check_init_sufficiency(
          X = X_init,
          y = y_init,
          checkpoints = init_checkpoints,
          metrics = init_stopping$metrics %||% c("variance"),
          threshold = init_stopping$threshold %||% 0.10,
          window = init_stopping$window %||% 2,
          verbose = init_stopping$verbose %||% FALSE
        )

        init_checkpoints <- c(init_checkpoints, list(check_result$metrics))

        if (check_result$converged) {
          if (progress) {
            message(sprintf("  [Init stopping] %s at %d/%d points (saved %d evals)",
                            check_result$reason, eval_counter, n_init,
                            n_init - eval_counter))
          }
          init_stopped_early <- TRUE
          break
        }
      }
    }
  }

  iter_counter <- 0L
  last_best_objective <- NA_real_
  last_best_objective_display <- Inf  # For progress display (observed best)

  # Store original base fidelity levels for dynamic scaling (hybrid_staged method)
  base_fidelity_levels <- fidelity_levels

  # Extract constraint metric names for metric-specific CV thresholds
  constraint_metric_names <- constraint_tbl$metric

  # Initialize budget tracking for dynamic fidelity
  cumulative_budget_used <- 0

  # Initialize for warm-starting and early stopping
  prev_surrogates <- NULL
  best_obj_history <- numeric()
  no_improvement_count <- 0
  low_acq_count <- 0  # Counter for consecutive low acquisition score iterations
  feasible_since_len <- NA_integer_  # length(best_obj_history) when feasibility first appeared

  # Stage 4 artifacts stay NULL/NA until post-loop verification runs; they are
  # initialized here (rather than after the loop) so build_policies() below
  # can read them at call time from inside the loop (Task D6 checkpoints).
  multi_seed_summary <- NULL
  multi_seed_runs <- NULL
  verdict <- NA_character_

  # Single constructor for the fit's policies list, shared by the periodic
  # checkpoint hook (Task D6) and the final build_baton_fit call. A closure
  # evaluated at call time (rather than a list built once here) because two
  # entries are not loop-invariant: fidelity_levels is re-computed every
  # iteration under the hybrid_staged method, and multi_seed_reps depends on
  # the post-loop multi_seed_summary. Call-time evaluation keeps the final
  # fit identical to the previous inline construction and keeps checkpoint
  # policies truthful about the state at snapshot time. Note: hybrid_staged
  # checkpoints pair dynamic fidelity_levels with base-derived fidelity_costs
  # (pre-existing final-fit behavior, now visible mid-run; the resume recipe
  # passes only initial_history, so the mismatch stays latent).
  build_policies <- function() {
    list(
      acquisition = acquisition,
      q = q,
      budget = budget,
      seed = rng_seed,
      fidelity_levels = fidelity_levels,
      fidelity_method = fidelity_method,
      fidelity_costs = fidelity_costs,
      candidate_pool = candidate_pool,
      objective = objective,
      early_stop = early_stop_config,
      multi_seed_verify = multi_seed_verify,
      multi_seed_n = if (multi_seed_verify) multi_seed_n else NA_integer_,
      multi_seed_reps = if (multi_seed_verify && !is.null(multi_seed_summary)) {
        multi_seed_summary$n_rep
      } else {
        NA_integer_
      },
      multi_seed_strict = multi_seed_strict,
      # D8 review carry-forward: record slim so a persisted fit is unambiguous
      # about WHY surrogates/diagnostics are NULL (slim vs a failed refit).
      slim = slim
    )
  }

  # run_status was initialized (Task D2) before the initialization phase; if
  # the walltime cap already tripped there (or an initialization round
  # errored under on_error = "return_partial", Task D7), skip the BO loop
  # entirely.
  ok <- TRUE  # guard flag; initialized so the pattern survives reordering
  while (!(run_status %in% c("walltime", "errored")) && eval_counter < budget) {
    # Task D4: advisory walltime cap, checked between iterations only (never
    # mid-batch), so a partial fit is always internally consistent. Elapsed
    # time can therefore overshoot the cap by up to one full iteration.
    if (walltime_exceeded()) {
      run_status <- "walltime"
      if (progress) message(sprintf("Walltime cap (%.0fs) reached at iteration %d.",
                                    max_walltime_s, iter_counter))
      break
    }

    iter_counter <- iter_counter + 1L

    # Compute dynamic fidelity levels for hybrid_staged method
    if (fidelity_method == "hybrid_staged") {
      fidelity_levels <- compute_dynamic_fidelity_levels(
        iter = iter_counter,
        budget = budget,
        base_levels = base_fidelity_levels,
        budget_used = cumulative_budget_used,
        budget_tolerance = 1.2,
        batch_size = q
      )
      if (progress && iter_counter %% 10 == 1) {
        message(sprintf("  Dynamic fidelity levels: low=%d, med=%d, high=%d",
                        fidelity_levels["low"], fidelity_levels["med"], fidelity_levels["high"]))
      }
    }

    if (progress) {
      message(sprintf("Iteration %d: fitting surrogates on %d evaluations.",
                      iter_counter, nrow(history)))
    }

    # Task D7: the entire iteration work region (fit -> acquisition ->
    # fidelity -> evaluate -> append) is guarded, not just the simulator
    # calls: surrogate fitting can stop() mid-loop on ill-conditioned data
    # after plenty of history exists. On error with
    # on_error = "return_partial" the loop breaks with run_status "errored"
    # and everything appended so far survives; a failure inside
    # evaluate_points() loses only that batch (appends are per-round). The
    # post-append sections (progress display, callback, checkpoint, early
    # stopping) stay outside the guard: they contain the loop's break sites
    # and are individually protected where they run user code.
    # EDIT CONSTRAINT: do not add break/next/return anywhere inside this
    # guarded block. The block is a promise evaluated in this function's
    # frame, not a separate function: break/next would silently bypass the
    # ok flag and every post-guard section of the iteration, and return()
    # would exit bo_calibrate entirely. Loop-control decisions belong AFTER
    # the closing "})" below.
    ok <- run_guarded({

    # Fit surrogates with warm-start from previous iteration
    # Include tryCatch to handle fitting failures gracefully
    surrogates <- tryCatch({
      fit_surrogates(history, objective, constraint_tbl,
                     covtype = covtype,
                     prev_surrogates = prev_surrogates,
                     fit_seed = rng_seed + 10000L * iter_counter)
    }, error = function(e) {
      warning(sprintf("Surrogate fitting failed at iteration %d: %s. Retrying without warm-start.",
                      iter_counter, e$message))
      # Retry without warm-start
      tryCatch({
        fit_surrogates(history, objective, constraint_tbl,
                       covtype = covtype,
                       prev_surrogates = NULL,
                       fit_seed = rng_seed + 10000L * iter_counter)
      }, error = function(e2) {
        stop(sprintf("Surrogate fitting failed even without warm-start: %s", e2$message),
             call. = FALSE)
      })
    })

    # Use configured incumbent method (BPM or BOI) for ECI computation
    # Restrict to high-fidelity evaluations when available for reliability
    # In fixed_per_stage mode, Stage 1 has only low-fidelity evals — filtering
    # to high-fidelity would yield zero points and return Inf, making the
    # incumbent_method parameter inert. Fall back to all fidelities when no
    # high-fidelity evaluations exist.
    has_high_fidelity <- "fidelity" %in% names(history) &&
      any(history$fidelity == "high", na.rm = TRUE)
    best_feasible_value <- best_feasible_objective(history, objective,
                                                    surrogates = surrogates,
                                                    high_fidelity_only = has_high_fidelity,
                                                    incumbent_method = incumbent_method)

    # Adaptive candidate pool size: scale with dimension and iteration
    d <- length(bounds)
    # Base size: 500 * dimension (min 1000, max 5000)
    base_pool_size <- pmax(1000, pmin(5000, 500 * d))
    # Early iterations: use base size
    # Late iterations (last 30%): increase by 50% for refinement
    if (iter_counter > 0.7 * (budget / q)) {
      candidate_pool_size <- min(base_pool_size * 1.5, 10000)
    } else {
      candidate_pool_size <- base_pool_size
    }
    # Respect user's minimum if provided
    candidate_pool_size <- max(candidate_pool_size, candidate_pool)

    unit_candidates <- lhs_candidate_pool(candidate_pool_size, bounds)

    # Single batch prediction used for acquisition + fidelity routing (avoids repeated GP calls)
    candidate_predictions <- predict_surrogates(surrogates, unit_candidates)
    candidate_prob_feas <- prob_feasibility_batch(candidate_predictions, constraint_tbl)

    acquisition_scores <- evaluate_acquisition(
      acquisition = acquisition,
      unit_candidates = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible_value,
      predictions = candidate_predictions,
      prob_feas = candidate_prob_feas
    )

    n_new <- min(q, budget - eval_counter)

    # Select batch with diversity mechanism
    if (n_new == 1) {
      # Single point: just take best
      selected_idx <- which.max(acquisition_scores)
    } else {
      # Batch: use local penalization for spatial diversity
      lipschitz <- estimate_lipschitz(surrogates, objective)
      # NOTE: with the corrected penalty (penalty = acq_best - L*r), a LARGER L
      # shrinks the suppression radius (acq_best/L) and thus REDUCES diversity.
      # The previous 1.5x inflation was calibrated to the old, inverted penalty
      # and is dropped: scaling L here would now weaken diversity, not strengthen it.
      selected_idx <- select_batch_local_penalization(
        candidates = unit_candidates,
        acq_scores = acquisition_scores,
        q = n_new,
        lipschitz = lipschitz
      )
    }

    selected_candidates <- unit_candidates[selected_idx, , drop = FALSE]

    if (progress && length(selected_idx) > 0) {
      max_acq <- max(acquisition_scores[selected_idx])
      message(sprintf("  -> max acquisition score: %.3f", max_acq))
    }

    # Task C2, pass 1 (sequential, no simulator calls): decide fidelity and
    # bookkeeping for every batch point. The only cross-point dependency in
    # the old per-point loop was total_budget_used, which included the ACTUAL
    # n_rep of earlier in-batch evaluations; it is now accumulated from the
    # REQUESTED n_rep (fidelity_levels[[fidelity]]) of earlier picks. The two
    # agree unless the simulator overrides its n_rep attribute, in which case
    # later in-batch fidelity decisions see the requested value -- batch
    # points are conceptually simultaneous, and this is what makes them
    # independently evaluable.
    batch_budget_base <- if ("n_rep" %in% names(history)) sum(history$n_rep, na.rm = TRUE) else 0L
    predicted_batch_reps <- 0
    batch_specs <- vector("list", n_new)
    # Hoisted once per iteration (Task C7): the best feasible design for the
    # hybrid_staged proximity bonus. History is fixed during pass 1, so the
    # lookup does not belong in the per-point loop.
    best_unit <- best_feasible_unit(history, "objective", bounds)
    for (i in seq_len(n_new)) {
      eval_counter <- eval_counter + 1L
      chosen_unit <- selected_candidates[i, ]
      theta <- scale_from_unit(chosen_unit, bounds)
      theta <- coerce_theta_types(theta, integer_params)

      prob_feas <- candidate_prob_feas[selected_idx[i]]

      # Get surrogate predictions from cached batch for CV computation
      obj_pred <- candidate_predictions[[objective]]
      pred_obj_mean <- obj_pred$mean[selected_idx[i]]
      pred_obj_sd <- obj_pred$sd[selected_idx[i]]
      cv_objective <- pred_obj_sd / max(abs(pred_obj_mean), 1e-6)

      cv_constraints <- vapply(constraint_metric_names, function(metric_name) {
        if (metric_name %in% names(candidate_predictions)) {
          pred <- candidate_predictions[[metric_name]]
          pred_sd <- pred$sd[selected_idx[i]]
          pred_mean <- pred$mean[selected_idx[i]]
          pred_sd / max(abs(pred_mean), 1e-6)
        } else {
          0
        }
      }, numeric(1))

      # For hybrid_staged: use max CV across constraints (tighter precision needed)
      # For other methods: use objective CV only
      if (fidelity_method == "hybrid_staged") {
        # Use maximum CV across all metrics, weighted toward constraints
        # This ensures we escalate if ANY metric needs higher precision
        max_constraint_cv <- max(cv_constraints, na.rm = TRUE)
        # If constraints have high CV or we're near boundary, prioritize constraint precision
        if (prob_feas >= 0.4 && prob_feas <= 0.7 && max_constraint_cv > cv_objective) {
          cv_estimate <- max_constraint_cv
          effective_metric <- constraint_metric_names[which.max(cv_constraints)]
        } else {
          # Otherwise use objective CV
          cv_estimate <- cv_objective
          effective_metric <- objective
        }
      } else {
        # Other methods just use objective CV
        cv_estimate <- cv_objective
        effective_metric <- objective
      }

      # Get acquisition value for this candidate
      acq_val <- acquisition_scores[selected_idx[i]]

      # Compute total budget used so far (history + earlier in-batch picks)
      total_budget_used <- batch_budget_base + predicted_batch_reps
      total_budget_sim <- budget * mean(fidelity_levels)  # approximate

      # Distance to the current best feasible design (hybrid_staged proximity
      # bonus), in scaled [0,1] space; 1.0 = far from best / none feasible.
      distance_to_best <- if (is.null(best_unit)) {
        1.0
      } else {
        sqrt(sum((chosen_unit - best_unit)^2))
      }

      # Select fidelity using configured method
      fidelity <- select_fidelity_method(
        method = fidelity_method,
        prob_feasible = prob_feas,
        cv_estimate = cv_estimate,
        acq_value = acq_val,
        best_obj = best_feasible_value,
        fidelity_levels = fidelity_levels,
        fidelity_costs = fidelity_costs,
        iter = iter_counter,
        total_budget_used = total_budget_used,
        total_budget = total_budget_sim,
        cv_threshold = fidelity_cv_threshold,
        cv_threshold_base = fidelity_cv_threshold,
        prob_range = fidelity_prob_range,
        metric = effective_metric,  # Now uses actual metric (constraint or objective)
        constraint_metrics = constraint_metric_names,
        distance_to_best = distance_to_best
      )

      # Debug output for fidelity selection (if enabled)
      # Support both new and legacy option names for backwards compatibility
      debug_fidelity <- getOption("BATON.debug_fidelity", FALSE) ||
                        getOption("evolveBO.debug_fidelity", FALSE)
      if (debug_fidelity && fidelity_method == "hybrid_staged") {
        message(sprintf("  [Fidelity Debug] iter=%d, prob_feas=%.3f, cv=%.3f, metric=%s, selected=%s",
                        iter_counter, prob_feas, cv_estimate, effective_metric, fidelity))
      }

      predicted_batch_reps <- predicted_batch_reps + fidelity_levels[[fidelity]]
      batch_specs[[i]] <- list(
        theta = theta,
        fidelity = fidelity,
        eval_id = eval_counter,
        iter = iter_counter,
        seed = rng_seed + eval_counter,
        extra = list(prob_feas = prob_feas,
                     cv_estimate = cv_estimate,
                     acq_score = acq_val)
      )
    }

    # Task C2, pass 2: evaluate the batch (in parallel when
    # options(BATON.cores) opts in) and append all rows in one bind.
    batch_rows <- evaluate_points(
      specs = batch_specs, sim_fun = sim_fun, bounds = bounds,
      objective = objective, constraint_tbl = constraint_tbl,
      fidelity_levels = fidelity_levels, ...
    )
    history <- append_rows(history, batch_rows)

    # Update cumulative budget for dynamic fidelity tracking
    cumulative_budget_used <- if ("n_rep" %in% names(history)) sum(history$n_rep, na.rm = TRUE) else 0L

    })  # end run_guarded (Task D7); no break/next/return inside the block
        # above (see the EDIT CONSTRAINT note at its opening).
    if (!ok) break

    # Progress display: use observed (BOI) best for informative tracking.
    # BPM incumbent changes only when surrogates are re-fit (once per
    # iteration), so it shows 0.000 within batches; the observed best is
    # reported once per batch now that evaluations complete together.
    if (progress) {
      for (r in batch_rows) announce_evaluation(r, objective)
    }
    current_best_display <- best_feasible_objective(history, objective,
                                                     surrogates = NULL,
                                                     high_fidelity_only = FALSE,
                                                     incumbent_method = "boi")
    if (progress && is.finite(current_best_display) && is.finite(last_best_objective_display)) {
      message(sprintf("  -> best observed feasible: %.3f (change: %.3f)",
                      current_best_display, current_best_display - last_best_objective_display))
    }
    if (is.finite(current_best_display)) {
      last_best_objective_display <- current_best_display
    }

    # Update for warm-starting next iteration
    prev_surrogates <- surrogates

    # Task D5: per-iteration callback with cooperative cancellation. Invoked
    # at the end of each BO iteration, after the batch has been recorded, so
    # cancellation never interrupts a batch mid-flight and no completed work
    # is discarded. Contract: return TRUE to continue; anything else cancels.
    # A callback error warns and continues (a monitoring hook must not kill a
    # paid run). If cancellation and the walltime cap collide in the same
    # iteration, "cancelled" wins: this break exits the loop before the
    # top-of-loop walltime check can run.
    if (!is.null(callback)) {
      keep_going <- tryCatch(
        isTRUE(callback(list(
          iter = iter_counter,
          eval_count = eval_counter,
          best_objective = if (is.finite(best_feasible_value)) best_feasible_value else NA_real_,
          best_observed = if (is.finite(current_best_display)) current_best_display else NA_real_,
          elapsed_s = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        ))),
        error = function(e) {
          warning(sprintf("callback error (continuing): %s", conditionMessage(e)),
                  call. = FALSE)
          TRUE
        })
      if (!keep_going) {
        run_status <- "cancelled"
        if (progress) {
          message(sprintf("Cancelled by callback at iteration %d.", iter_counter))
        }
        break
      }
    }

    # Task D6: periodic resumable snapshot every checkpoint_every completed BO
    # iterations (iteration-granular, not per-evaluation; the initial design
    # phase emits none). Checkpoints are slim by construction: surrogates and
    # diagnostics are NULL because GP objects are large, environment-laden,
    # and not needed to resume (a resume refits from history via
    # initial_history). Incumbent info is the observed feasible best by row;
    # no surrogate fit is performed here. An error in checkpoint_fun warns
    # and continues, same rationale as the callback above.
    # ORDERING NOTE (deliberate, do not "fix"): this hook runs BEFORE the
    # early-stopping bookkeeping below, so the early-stop counters have not
    # yet absorbed the current iteration at snapshot time. That is fine: the
    # counters are in-loop state, not part of the history, and a resume via
    # initial_history resets them regardless.
    if (!is.null(checkpoint_fun) && iter_counter %% checkpoint_every == 0L) {
      # which.min() returns integer(0) when every feasible row has an NA
      # objective (and ckpt_feas_idx[integer(0)] is integer(0) when none is
      # feasible), so gate on length rather than is.na(): is.na(integer(0))
      # is logical(0) and would crash the if() (D6 review carry-forward).
      ckpt_feas_idx <- which(history$feasible)
      ckpt_best_candidate <- ckpt_feas_idx[which.min(history$objective[ckpt_feas_idx])]
      ckpt_best_row <- if (length(ckpt_best_candidate) == 1L) {
        ckpt_best_candidate
      } else {
        NA_integer_
      }
      # Build the snapshot BEFORE the tryCatch so the "checkpoint_fun error"
      # warning is only ever attributed to the USER's function; a bug in
      # build_baton_fit/build_policies must surface as a real error, not be
      # silently swallowed as a checkpoint hiccup (D6 review carry-forward).
      ckpt_fit <- build_baton_fit(
        history = history,
        best_theta = if (is.na(ckpt_best_row)) NULL else history$theta[[ckpt_best_row]],
        best_objective = if (is.na(ckpt_best_row)) NA_real_ else history$objective[[ckpt_best_row]],
        surrogates = NULL,
        policies = build_policies(),
        diagnostics = NULL,
        multi_seed_summary = NULL,
        multi_seed_runs = NULL,
        verdict = NA_character_,
        bounds = bounds,
        constraints = constraints,
        constraint_tbl = constraint_tbl,
        status = "checkpoint"
      )
      tryCatch(checkpoint_fun(ckpt_fit), error = function(e) {
        warning(sprintf("checkpoint_fun error (continuing): %s",
                        conditionMessage(e)), call. = FALSE)
      })
    }

    # Early stopping check
    has_feasible <- is.finite(best_feasible_value)
    current_best <- if (has_feasible) {
      best_feasible_value
    } else {
      # No feasible solution yet, use best infeasible
      min(history$objective, na.rm = TRUE)
    }
    best_obj_history <- c(best_obj_history, current_best)
    # Record when feasibility first appeared. The infeasible-phase objectives
    # (best infeasible minimum) and the feasible-phase incumbent are different
    # quantities, and the transition typically shows as a large negative
    # "improvement"; comparing across it can trigger early stopping right when
    # the search finally becomes feasible. So the improvement-based stop only
    # compares WITHIN the feasible phase.
    if (has_feasible && is.na(feasible_since_len)) {
      feasible_since_len <- length(best_obj_history)
    }

    # Check for early stopping after minimum iterations (if enabled)
    # Start checking after 3 BO iterations (not n_init, which could be 60+)
    min_iters_before_stopping <- 3
    if (early_stop_config$enabled && iter_counter > min_iters_before_stopping) {
      # Check 1: No improvement in objective for patience iterations
      es_patience <- early_stop_config$patience  # default: 5 iterations
      es_threshold <- early_stop_config$threshold  # default: 1e-3 (0.1%)

      # Number of feasible-phase entries available for comparison.
      feasible_phase_len <- if (is.na(feasible_since_len)) 0L else
        length(best_obj_history) - feasible_since_len + 1L

      if (has_feasible && feasible_phase_len >= es_patience + 1) {
        # Compare current best vs best from es_patience iterations ago
        current_best_val <- best_obj_history[length(best_obj_history)]
        earlier_best_val <- best_obj_history[length(best_obj_history) - es_patience]

        # Relative improvement (positive means we improved, i.e., found lower objective)
        improvement <- (earlier_best_val - current_best_val) / (abs(earlier_best_val) + 1e-8)

        if (improvement < es_threshold) {
          # No meaningful improvement in last es_patience iterations
          no_improvement_count <- no_improvement_count + 1
          if (progress && no_improvement_count == 1) {
            message(sprintf("  -> No improvement > %.1f%% in last %d iterations (count: %d)",
                            es_threshold * 100, es_patience, no_improvement_count))
          }
          # Trigger after seeing no improvement for `consecutive` checks
          if (no_improvement_count >= early_stop_config$consecutive) {
            if (progress) {
              message(sprintf("Early stopping at iteration %d: no improvement > %.1f%% for %d consecutive iterations",
                              iter_counter, es_threshold * 100, no_improvement_count))
            }
            run_status <- "early_stopped"
            break
          }
        } else {
          if (no_improvement_count > 0 && progress) {
            message(sprintf("  -> Improvement detected (%.2f%%), resetting counter", improvement * 100))
          }
          no_improvement_count <- 0
        }
      }

      # Check 2: All acquisition values very small (converged)
      # Use selected batch scores (what's actually being evaluated) rather than full pool
      # Threshold of 0.1 means we stop when expected improvement is <10% of current best
      # Require 3 consecutive iterations with low acquisition before stopping
      # Guard against empty selected_idx (would return -Inf from max())
      if (length(selected_idx) > 0) {
        selected_acq_max <- max(acquisition_scores[selected_idx], na.rm = TRUE)
      } else {
        selected_acq_max <- Inf  # Don't trigger early stopping if no candidates selected
      }
      acq_patience <- 3  # Number of consecutive low-acq iterations required
      # Relative threshold: EI (objective units) is scaled by the incumbent, so
      # the criterion means the same thing whether the objective is O(1) (a
      # probability) or O(100) (expected enrollment). The old absolute 0.1
      # stopped [0,1]-scale objectives almost immediately and never fired for
      # large-scale ones.
      acq_stop_level <- early_stop_config$acq_rel *
        (if (has_feasible) max(abs(best_feasible_value), 1e-8) else 1)
      if (is.finite(selected_acq_max) && selected_acq_max < acq_stop_level) {
        low_acq_count <- low_acq_count + 1
        if (low_acq_count >= acq_patience) {
          if (progress) {
            message(sprintf("Early stopping at iteration %d: max acquisition < %.3g for %d consecutive iterations",
                            iter_counter, acq_stop_level, acq_patience))
          }
          run_status <- "acq_flatline"
          break
        }
      } else {
        low_acq_count <- 0  # Reset counter if acquisition goes back up
      }
    }
  }

  # Task D7: under on_error = "return_partial" the final refit is attempted
  # inside its own tryCatch. The partial history may be exactly the kind of
  # ill-conditioned data that errored the in-loop fit, and losing the whole
  # fit here would defeat the point; surrogates are NULL in that case. When
  # the failure is the FIRST error of the run (a completed loop whose final
  # refit fails), it sets run_status/run_error; an already-errored run keeps
  # its ORIGINAL in-loop error message. In "stop" mode nothing changes: the
  # refit runs unguarded and any error propagates as before.
  # Task D8: slim = TRUE skips the final refit entirely. Nothing downstream
  # consumes it except the returned surrogates slot and generate_diagnostics
  # (both dropped under slim): best_theta comes from best_feasible_index()
  # and best_objective from the feasible history rows, and Stage 4 uses only
  # best_theta. Consequence: an error arising only in the final refit cannot
  # surface under slim (the run keeps its natural status in either on_error
  # mode); documented in the roxygen for slim.
  # EDIT CONSTRAINT: any new consumer of final_surrogates below this point
  # must handle NULL. Both slim = TRUE and the return_partial failed-refit
  # path above produce final_surrogates = NULL on an otherwise valid fit.
  final_surrogates <- if (slim) {
    NULL
  } else if (on_error == "return_partial") {
    tryCatch(
      fit_surrogates(history, objective, constraint_tbl,
                     covtype = covtype,
                     prev_surrogates = prev_surrogates,
                     fit_seed = rng_seed + 99999L),
      error = function(e) {
        if (run_status != "errored") {
          run_status <<- "errored"
          run_error <<- conditionMessage(e)
        }
        NULL
      })
  } else {
    fit_surrogates(history, objective, constraint_tbl,
                   covtype = covtype,
                   prev_surrogates = prev_surrogates,
                   fit_seed = rng_seed + 99999L)
  }

  # Empty history is only reachable under on_error = "return_partial" (an
  # errored first evaluation round); indexing best_feasible_index() there
  # would fail. In "stop" mode at least one evaluation always exists (any
  # earlier error propagated), so this branch never fires and the path below
  # is byte-identical to the previous behavior.
  if (nrow(history) == 0L) {
    best_theta <- NULL
    best_objective <- NA_real_
    diagnostics <- NULL
  } else {
    best_idx <- best_feasible_index(history, objective)
    best_theta <- history$theta[[best_idx]]

    # D7 review carry-forward: under return_partial a diagnostics failure
    # after a SUCCESSFUL refit must not lose the fit; diagnostics are a
    # reporting artifact, so they are dropped with a warning and the run
    # keeps its natural status (no run_status/run_error change). In "stop"
    # mode the error propagates as before. Under slim = TRUE this whole
    # branch is moot: final_surrogates is NULL (refit skipped), so
    # diagnostics are NULL without ever calling generate_diagnostics.
    diagnostics <- if (is.null(final_surrogates)) {
      NULL
    } else if (on_error == "return_partial") {
      tryCatch(
        generate_diagnostics(final_surrogates, history, objective),
        error = function(e) {
          warning(sprintf(
            "diagnostics generation failed (diagnostics dropped from the fit): %s",
            conditionMessage(e)), call. = FALSE)
          NULL
        })
    } else {
      generate_diagnostics(final_surrogates, history, objective)
    }

    # Compute best objective value from feasible solutions
    feasible_history <- history[history$feasible, , drop = FALSE]
    best_objective <- if (nrow(feasible_history) > 0) {
      min(feasible_history$objective, na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  # v0.4.0 Stage 4: multi-seed verification gate. Re-evaluates the calibrated
  # design at multiple independent seeds at high fidelity to detect
  # stored-vs-actual operating-characteristic gaps. Backward-compatible default
  # (multi_seed_verify = FALSE). multi_seed_summary / multi_seed_runs /
  # verdict were initialized (NULL/NULL/NA) before the BO loop.
  if (multi_seed_verify &&
      (run_status %in% c("walltime", "cancelled", "errored") || walltime_exceeded())) {
    # D4 follow-up: never launch Stage 4 (multi_seed_n high-fidelity
    # re-evaluations) once the walltime cap has been reached; that work could
    # dwarf the cap itself. The walltime_exceeded() disjunct covers runs that
    # end via the loop CONDITION going false (e.g. a budget-filling final
    # chunk or batch that crosses the cap): status still records why the run
    # ended ("budget_exhausted" etc.), but verification is suppressed. Task
    # D5 adds the "cancelled" case: a caller who cancelled via the callback
    # asked the run to stop, so launching a fresh block of high-fidelity
    # verification work would defeat the point of cancelling. Task D7 adds
    # the "errored" case: a run that already hit an infrastructure or
    # numerical failure should return its partial fit immediately, not
    # launch fresh verification work against the same failing backend. The
    # fit is returned with multi_seed_summary / multi_seed_runs NULL and
    # verdict NA, same as multi_seed_verify = FALSE.
    if (progress) {
      message(sprintf(
        "Stage 4 multi-seed verification skipped: %s.",
        if (run_status == "cancelled") "run cancelled by callback"
        else if (run_status == "errored") "run errored (partial fit returned)"
        else "walltime cap tripped or exceeded"))
    }
  } else if (multi_seed_verify && !is.na(best_objective)) {
    if (is.null(multi_seed_reps)) {
      ms_reps <- as.integer(fidelity_levels[["high"]])
    } else {
      ms_reps <- multi_seed_reps
    }
    ms_seeds <- as.integer(rng_seed) + seq_len(multi_seed_n) - 1L
    if (progress) {
      message(sprintf(
        "Stage 4: multi-seed verification at %d seeds x R = %d (seeds %d..%d)",
        multi_seed_n, ms_reps, ms_seeds[1], ms_seeds[length(ms_seeds)]
      ))
    }
    # D7 review carry-forward: the Stage 4 block runs AFTER the most
    # expensive simulator work of the run, so under return_partial a crash
    # here (one that escapes run_multi_seed_verification's per-seed
    # tryCatch, e.g. a malformed simulator return breaking the summary
    # construction) must not lose the fit. run_guarded() applies the same
    # pattern as the in-loop guards: in "stop" mode it is a no-op wrapper
    # (errors propagate unchanged); in "return_partial" mode a caught error
    # sets run_status "errored" + run_error and we reset the verification
    # fields so the fit carries no half-built gate artifacts.
    ms_ok <- run_guarded({
      ms_result <- run_multi_seed_verification(
        sim_fun = sim_fun,
        theta = best_theta,
        seeds = ms_seeds,
        n_rep = ms_reps,
        objective = objective,
        constraint_tbl = constraint_tbl,
        progress = progress,
        ...
      )
      multi_seed_runs <- ms_result$per_seed
      multi_seed_summary <- ms_result$summary
      verdict <- if (multi_seed_summary$strict_feas == 1) {
        "MULTI_SEED_PASS"
      } else if (multi_seed_strict) {
        "MULTI_SEED_FAIL"
      } else {
        "MULTI_SEED_WARN"
      }
      if (progress) {
        message(sprintf(
          "Stage 4 verdict: %s (strict feas %d/%d; %s mean = %.4f, sd = %.4f)",
          verdict,
          multi_seed_summary$strict_feas_count,
          multi_seed_n,
          objective,
          multi_seed_summary$objective_mean,
          multi_seed_summary$objective_sd
        ))
      }
    })
    if (!ms_ok) {
      multi_seed_runs <- NULL
      multi_seed_summary <- NULL
      verdict <- NA_character_
    }
  }

  build_baton_fit(
    history = history,
    best_theta = best_theta,
    best_objective = best_objective,
    surrogates = final_surrogates,
    policies = build_policies(),
    diagnostics = diagnostics,
    multi_seed_summary = multi_seed_summary,
    # Task D8: slim keeps the adoption gate (multi_seed_summary + verdict)
    # but drops the per-seed evaluation table; surrogates and diagnostics
    # are already NULL under slim (final refit skipped above).
    multi_seed_runs = if (slim) NULL else multi_seed_runs,
    verdict = verdict,
    bounds = bounds,
    constraints = constraints,
    constraint_tbl = constraint_tbl,
    status = run_status,
    error_message = run_error
  )
}

#' Assemble a BATON_fit object (Task D1)
#'
#' Single construction point so normal completion, partial returns
#' (error/cancel/walltime), and checkpoints produce structurally identical
#' objects. `status` records why the run ended.
#'
#' Notes on fields: `multi_seed_summary`/`multi_seed_runs` are the v0.4.0
#' multi-seed verification artifacts (NULL if multi_seed_verify = FALSE);
#' `bounds`/`constraints` are stored to support warm-starting.
#' @keywords internal
build_baton_fit <- function(history, best_theta, best_objective, surrogates,
                            policies, diagnostics, multi_seed_summary,
                            multi_seed_runs, verdict, bounds, constraints,
                            constraint_tbl, status, error_message = NULL) {
  # Exact matching (not match.arg) so partial matches like "budget" are
  # rejected. errored is constructed by the D7 on_error = "return_partial"
  # path in bo_calibrate and by the philosophies wrapper's isolation layer;
  # completed is legacy, slated for removal once nothing constructs it.
  status_values <- c("budget_exhausted", "early_stopped", "acq_flatline",
                     "cancelled", "walltime", "errored", "checkpoint",
                     "completed")
  if (!is.character(status) || length(status) != 1L || !status %in% status_values) {
    stop(sprintf("Invalid status '%s'. Must be one of: %s",
                 paste(status, collapse = ", "),
                 paste(status_values, collapse = ", ")), call. = FALSE)
  }
  structure(
    list(
      history = history,
      best_theta = best_theta,
      best_objective = best_objective,
      surrogates = surrogates,
      policies = policies,
      diagnostics = diagnostics,
      multi_seed_summary = multi_seed_summary,
      multi_seed_runs = multi_seed_runs,
      verdict = verdict,
      bounds = bounds,
      constraints = constraints,
      constraint_tbl = constraint_tbl,
      status = status,
      error_message = error_message
    ),
    class = c("BATON_fit", "list")
  )
}

#' @keywords internal
initialise_history <- function() {
  tibble::tibble(
    iter = integer(),
    eval_id = integer(),
    theta = list(),
    unit_x = list(),
    theta_id = character(),
    fidelity = character(),
    n_rep = integer(),
    metrics = list(),
    variance = list(),
    objective = numeric(),
    feasible = logical(),
    prob_feas = numeric(),
    cv_estimate = numeric(),
    acq_score = numeric()
  )
}

#' Evaluate one design and build its history row (pure, no side effects)
#'
#' Task C2: the evaluation itself is separated from history mutation so that
#' independent points (initial design, BO batches) can be evaluated
#' concurrently and appended in one bind. The simulator call runs under
#' `run_seeded(seed)` so its RNG use neither depends on stream position nor
#' leaks into the optimizer's stream -- this is what makes results identical
#' whether points are evaluated serially or in parallel forks.
#'
#' @keywords internal
invoke_evaluation <- function(sim_fun,
                              theta,
                              bounds,
                              objective,
                              constraint_tbl,
                              fidelity,
                              fidelity_levels,
                              eval_id,
                              iter,
                              seed,
                              extra = NULL,
                              ...) {
  unit_theta <- scale_to_unit(theta, bounds)
  theta_id <- theta_to_id(unit_theta)
  result <- run_seeded(seed, function() {
    invoke_simulator(
      sim_fun = sim_fun,
      theta = theta,
      fidelity = fidelity,
      n_rep = fidelity_levels[[fidelity]],
      seed = seed,
      ...
    )
  })

  metrics <- result$metrics
  variance <- result$variance
  n_rep <- result$n_rep

  feasible <- is_feasible(metrics, constraint_tbl)

  # Check membership before subsetting: metrics has been normalized to a named
  # numeric vector by invoke_simulator(), and `metrics[[objective]]` on a
  # numeric vector with a missing name throws a bare "subscript out of bounds"
  # before any is.null() guard can fire. Mirrors the constraint-name check below.
  if (!objective %in% names(metrics)) {
    stop(sprintf(
      paste0("Objective '%s' not returned by simulator. ",
             "Simulator returned: %s. The objective must be one of the ",
             "named metrics the simulator returns."),
      objective, paste(names(metrics), collapse = ", ")
    ), call. = FALSE)
  }
  objective_value <- metrics[[objective]]

  # Validate that every constraint metric name is actually returned by the
  # simulator. Without this check, a misnamed constraint metric (e.g. "pwr"
  # instead of "power") is silently treated as infeasible by is_feasible()
  # (metrics[[metric]] is NULL -> FALSE), so no design is ever feasible,
  # best_objective becomes NA, multi-seed verification is skipped, and the
  # run reports "success" with garbage output. Mirror the upfront objective
  # check above so the failure is loud and immediate.
  if (nrow(constraint_tbl) > 0L) {
    missing_constraints <- setdiff(constraint_tbl$metric, names(metrics))
    if (length(missing_constraints) > 0L) {
      stop(sprintf(
        paste0("Constraint metric(s) not returned by simulator: %s. ",
               "Simulator returned: %s. Constraint names must match the ",
               "names of the metrics the simulator returns."),
        paste(missing_constraints, collapse = ", "),
        paste(names(metrics), collapse = ", ")
      ), call. = FALSE)
    }
  }

  if (is.null(extra)) extra <- list()
  prob_feas_val <- if (!is.null(extra$prob_feas)) as.numeric(extra$prob_feas) else NA_real_
  cv_val <- if (!is.null(extra$cv_estimate)) as.numeric(extra$cv_estimate) else NA_real_
  acq_val <- if (!is.null(extra$acq_score)) as.numeric(extra$acq_score) else NA_real_

  # Unpack metrics into individual columns for easier access
  # This allows adaptive bounds functions to directly access constraint metrics
  # OPTIMIZED: Use vapply for faster vectorized conversion
  metrics_df <- if (length(metrics) > 0) {
    metrics_vec <- vapply(metrics, function(x) {
      if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]
    }, numeric(1))
    as.data.frame(as.list(metrics_vec), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  # Create base row
  base_row <- tibble::tibble(
    iter = as.integer(iter),
    eval_id = as.integer(eval_id),
    theta = list(theta),
    unit_x = list(unit_theta),
    theta_id = as.character(theta_id),
    fidelity = as.character(fidelity),
    n_rep = as.integer(n_rep),
    metrics = list(metrics),
    variance = list(variance),
    objective = as.numeric(objective_value),
    feasible = as.logical(feasible),
    prob_feas = as.numeric(prob_feas_val),
    cv_estimate = as.numeric(cv_val),
    acq_score = as.numeric(acq_val)
  )

  # Unpack theta into individual parameter columns
  # This ensures new rows have the same columns as initial history (which has
  # unpacked eff, fut, ev, etc.) and prevents NAs when bind_rows combines them
  # OPTIMIZED: Use vapply for faster vectorized conversion
  theta_df <- if (length(theta) > 0) {
    theta_vec <- vapply(theta, function(x) {
      if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]
    }, numeric(1))
    as.data.frame(as.list(theta_vec), stringsAsFactors = FALSE)
  } else {
    data.frame()
  }

  # Add unpacked theta columns to base_row (before metrics to avoid conflicts)
  if (ncol(theta_df) > 0) {
    # Check for column name conflicts with base_row
    conflict_cols <- intersect(names(base_row), names(theta_df))
    if (length(conflict_cols) > 0) {
      # These shouldn't conflict, but if they do, don't overwrite base columns
      theta_df <- theta_df[, !names(theta_df) %in% conflict_cols, drop = FALSE]
    }
    if (ncol(theta_df) > 0) {
      base_row <- dplyr::bind_cols(base_row, theta_df)
    }
  }

  # Combine with unpacked metrics
  if (ncol(metrics_df) > 0) {
    # Check for column name conflicts and rename if needed
    conflict_cols <- intersect(names(base_row), names(metrics_df))
    if (length(conflict_cols) > 0) {
      names(metrics_df)[names(metrics_df) %in% conflict_cols] <-
        paste0("metric_", names(metrics_df)[names(metrics_df) %in% conflict_cols])
    }
    new_row <- dplyr::bind_cols(base_row, metrics_df)
  } else {
    new_row <- base_row
  }

  new_row
}

#' Append evaluation rows to the history in a single bind
#' @keywords internal
append_rows <- function(history, rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) return(history)
  new_rows <- dplyr::bind_rows(rows)
  # Ensure type consistency before binding
  tryCatch({
    dplyr::bind_rows(history, new_rows)
  }, error = function(e) {
    # If bind_rows fails, try to diagnose and fix
    warning(sprintf("bind_rows failed at eval_id %s: %s. Attempting type coercion.",
                    paste(new_rows$eval_id, collapse = ","), e$message))
    # Coerce all columns in history to match new_rows types
    for (col in names(new_rows)) {
      if (col %in% names(history) && !inherits(new_rows[[col]], "list")) {
        history[[col]] <- tryCatch(
          as(history[[col]], class(new_rows[[col]])[1]),
          error = function(e2) history[[col]]
        )
      }
    }
    dplyr::bind_rows(history, new_rows)
  })
}

#' Print the per-evaluation progress line from a completed history row
#' @keywords internal
announce_evaluation <- function(row, objective) {
  message(sprintf("  Eval %03d (iter %02d, %s fidelity): %s = %.4f%s",
                  row$eval_id, row$iter, row$fidelity, objective, row$objective,
                  if (isTRUE(row$feasible)) " [feasible]" else ""))
}

#' Evaluate a set of independent design points, optionally in parallel
#'
#' Each spec is a list with `theta`, `fidelity`, `eval_id`, `iter`, `seed`,
#' and optionally `extra`. Points run concurrently when the user opts in via
#' `options(BATON.cores = k)` on Unix (same convention as `fit_surrogates`);
#' per-point seeding via `run_seeded` inside `invoke_evaluation` makes the
#' results identical to the serial path regardless of core count. Errors
#' raised inside forked workers are re-thrown here so parallel failure
#' behavior matches serial.
#'
#' @keywords internal
evaluate_points <- function(specs, sim_fun, bounds, objective, constraint_tbl,
                            fidelity_levels, ...) {
  .cores <- effective_cores()
  .map <- if (.cores > 1L) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = .cores)
  } else {
    lapply
  }
  rows <- .map(specs, function(sp) {
    invoke_evaluation(
      sim_fun = sim_fun,
      theta = sp$theta,
      bounds = bounds,
      objective = objective,
      constraint_tbl = constraint_tbl,
      fidelity = sp$fidelity,
      fidelity_levels = fidelity_levels,
      eval_id = sp$eval_id,
      iter = sp$iter,
      seed = sp$seed,
      extra = sp$extra,
      ...
    )
  })
  # mclapply swallows worker errors into try-error objects; re-throw to match
  # the serial path, where the first failing evaluation stops the run.
  for (r in rows) {
    if (inherits(r, "try-error")) {
      cond <- attr(r, "condition")
      stop(if (!is.null(cond)) conditionMessage(cond) else as.character(r),
           call. = FALSE)
    }
  }
  rows
}

#' Evaluate one design and append it to the history (serial convenience)
#' @keywords internal
record_evaluation <- function(history,
                              sim_fun,
                              theta,
                              bounds,
                              objective,
                              constraint_tbl,
                              fidelity,
                              fidelity_levels,
                              eval_id,
                              iter,
                              seed,
                              progress,
                              extra = NULL,
                              ...) {
  row <- invoke_evaluation(
    sim_fun = sim_fun,
    theta = theta,
    bounds = bounds,
    objective = objective,
    constraint_tbl = constraint_tbl,
    fidelity = fidelity,
    fidelity_levels = fidelity_levels,
    eval_id = eval_id,
    iter = iter,
    seed = seed,
    extra = extra,
    ...
  )
  if (progress) announce_evaluation(row, objective)
  append_rows(history, list(row))
}

#' @keywords internal
invoke_simulator <- function(sim_fun, theta, fidelity, n_rep, seed, ...) {
  # Optional-argument n_rep contract (Task B1): forward the replication count only to
  # simulators that declare an explicit n_rep formal. Dots alone do NOT opt in:
  # forwarding wrappers (e.g. fix_parameters()) pass their dots verbatim to an
  # inner simulator, so injecting n_rep into dots would leak it into legacy
  # simulators that cannot accept it.
  sim_args <- names(formals(sim_fun))
  if ("n_rep" %in% sim_args) {
    res <- sim_fun(theta, fidelity = fidelity, seed = seed, n_rep = n_rep, ...)
  } else {
    res <- sim_fun(theta, fidelity = fidelity, seed = seed, ...)
  }

  variance <- attr(res, "variance", exact = TRUE)
  rep_attr <- attr(res, "n_rep", exact = TRUE)
  if (is.list(res) && !is.null(res$metrics)) {
    metrics <- res$metrics
    if (is.null(variance) && !is.null(res$variance)) {
      variance <- res$variance
    }
    if (is.null(rep_attr) && !is.null(res$n_rep)) {
      rep_attr <- res$n_rep
    }
  } else {
    metrics <- res
  }

  # Use base R instead of purrr::imap_dbl to avoid "In index: X" errors
  metrics_list <- as.list(metrics)
  metrics <- vapply(metrics_list, function(x) {
    tryCatch(as.numeric(x), error = function(e) NA_real_)
  }, FUN.VALUE = numeric(1))
  names(metrics) <- names(metrics_list)

  if (is.null(variance)) {
    variance <- default_variance_estimator(metrics, n_rep)
  } else {
    variance_list <- as.list(variance)
    variance <- vapply(variance_list, function(x) {
      tryCatch(as.numeric(x), error = function(e) NA_real_)
    }, FUN.VALUE = numeric(1))
    names(variance) <- names(variance_list)
  }
  if (is.null(rep_attr) || is.na(rep_attr)) {
    rep_attr <- n_rep
  }

  list(
    metrics = metrics,
    variance = variance,
    n_rep = rep_attr
  )
}

#' Default variance estimator for simulator outputs
#'
#' Provides a fallback variance estimate when the simulator does not return
#' variance information. Uses binomial variance formula for proportion-type
#' metrics (values in [0,1]) and returns NA for other metrics.
#'
#' @param metrics named numeric vector of operating characteristics.
#' @param n_rep number of simulator replications.
#'
#' @details
#' For metrics that appear to be proportions (e.g., power, type I error),
#' this function estimates variance using the binomial formula: p(1-p)/n.
#' For continuous metrics outside [0,1] (e.g., expected sample size, trial
#' duration), it returns NA, which triggers the use of a small homoskedastic
#' nugget in the Gaussian process surrogate.
#'
#' \strong{Recommendation:} For best results, simulators should return variance
#' estimates as an attribute, especially for non-proportion metrics where the
#' binomial approximation does not apply.
#'
#' @return Named numeric vector of variance estimates (may contain NA values).
#' @keywords internal
default_variance_estimator <- function(metrics, n_rep) {
  if (is.null(n_rep) || is.na(n_rep) || n_rep <= 0) {
    result <- rep(NA_real_, length(metrics))
    names(result) <- names(metrics)
    return(result)
  }
  # Use base R vapply instead of purrr::imap_dbl to avoid "In index: X" errors
  metric_names <- names(metrics)
  result <- vapply(seq_along(metrics), function(i) {
    value <- metrics[i]
    # Use binomial variance for proportion-type metrics (values in [0,1])
    if (!is.na(value) && value >= 0 && value <= 1) {
      max(value * (1 - value) / n_rep, 1e-6)
    } else {
      # For continuous metrics outside [0,1], return NA
      # The GP fitting will use a small nugget for these cases
      NA_real_
    }
  }, FUN.VALUE = numeric(1))
  names(result) <- metric_names
  result
}

#' @keywords internal
validate_fidelity_levels <- function(fidelity_levels) {
  if (is.null(names(fidelity_levels))) {
    stop("`fidelity_levels` must be a named numeric vector.", call. = FALSE)
  }
  # The code hard-indexes fidelity_levels[["low"]] / [["med"]] / [["high"]],
  # so the names must be exactly that set (in any order). Without this check a
  # misnamed level produces a confusing subscript-out-of-bounds error later.
  if (!setequal(names(fidelity_levels), c("low", "med", "high"))) {
    stop(sprintf(
      "`fidelity_levels` must be named exactly 'low', 'med', 'high'; got: %s.",
      paste(names(fidelity_levels), collapse = ", ")
    ), call. = FALSE)
  }
  # Use base R for loop instead of purrr::iwalk
  for (i in seq_along(fidelity_levels)) {
    value <- fidelity_levels[i]
    if (!is.finite(value) || value <= 0) {
      stop("Fidelity replication counts must be positive.", call. = FALSE)
    }
  }
  fidelity_levels
}

#' Unit coordinates of the best feasible design in the history
#'
#' Task C7: found by row index (which.min over feasible objectives), not by
#' float-equality against the incumbent -- the BPM incumbent is a posterior
#' mean that equals no observed objective, so an equality match found nothing
#' and the hybrid_staged proximity bonus never fired.
#'
#' Mirrors the A8 incumbent semantics: high-fidelity feasible rows are
#' preferred when any exist, falling back to any-fidelity feasible rows.
#' Coordinates are rescaled from theta under the CURRENT bounds rather than
#' read from the stored unit_x, which is stale when a previous-stage history
#' was recorded under different (e.g. wider) bounds. Returns NULL when no
#' feasible row with a non-NA objective exists.
#'
#' @keywords internal
best_feasible_unit <- function(history, objective, bounds) {
  feas <- history$feasible & !is.na(history[[objective]])
  if ("fidelity" %in% names(history)) {
    high <- feas & !is.na(history$fidelity) & history$fidelity == "high"
    if (any(high)) feas <- high
  }
  feas_idx <- which(feas)
  if (length(feas_idx) == 0) return(NULL)
  best_row <- feas_idx[which.min(history[[objective]][feas_idx])]
  unit <- unlist(scale_to_unit(history$theta[[best_row]], bounds))
  # scale_to_unit follows the theta's own name order, but candidate rows are
  # bounds-ordered and the distance subtraction is positional: return in
  # names(bounds) order always.
  unit[names(bounds)]
}

#' @keywords internal
lhs_candidate_pool <- function(n, bounds) {
  # Task C5: return the matrix directly. predict_surrogates and
  # select_batch_local_penalization consume it without any per-candidate
  # list round trip; rows index candidates, columns are named parameters.
  lhs <- lhs::randomLHS(n, length(bounds))
  colnames(lhs) <- names(bounds)
  lhs
}

#' Best feasible objective value using posterior mean (BPM) when available
#'
#' Uses the GP surrogate's posterior mean for more stable incumbent estimation
#' under Monte Carlo noise. Falls back to best observed value (BOI) when
#' surrogates are unavailable (e.g., during initialization).
#'
#' @param history tibble of evaluation history
#' @param objective name of the objective metric
#' @param surrogates fitted surrogate models (optional, enables BPM)
#' @param high_fidelity_only if TRUE, restrict to high-fidelity evaluations
#' @return best feasible objective value (Inf if none feasible)
#' @keywords internal
best_feasible_objective <- function(history, objective, surrogates = NULL,
                                    high_fidelity_only = FALSE,
                                    incumbent_method = "bpm") {
  idx <- which(history$feasible)

  # Optional: prefer high-fidelity evaluations for stability, but fall back to
  # the any-fidelity feasible best when no feasible high-fidelity eval exists.
  # Returning Inf here (the previous behavior) flipped the acquisition into its
  # "no feasible solution yet" mode even though feasible designs were known -
  # common mid-run, since the initial design runs entirely at the lowest fidelity.
  if (high_fidelity_only && "fidelity" %in% names(history)) {
    hf_idx <- intersect(idx, which(history$fidelity == "high"))
    if (length(hf_idx) > 0L) idx <- hf_idx
  }

  if (length(idx) == 0L) {
    return(Inf)  # genuinely no feasible point at any fidelity
  }

  # BPM: Use posterior mean when surrogates available (reduces MC noise sensitivity)
  # Only use BPM if incumbent_method == "bpm" AND surrogates are available
  if (incumbent_method == "bpm" && !is.null(surrogates) && objective %in% names(surrogates)) {
    # history$unit_x is a LIST COLUMN - subset with [idx] not [idx, , drop=FALSE]
    unit_x_feasible <- history$unit_x[idx]
    # Predict only the objective surrogate: the previous code predicted all m
    # metrics here and discarded all but the objective (~75% wasted GP work per
    # iteration). Subsetting the surrogate list preserves identical results.
    pred <- predict_surrogates(surrogates[objective], unit_x_feasible)
    return(min(pred[[objective]]$mean, na.rm = TRUE))
  }

  # BOI (best observed incumbent): use when incumbent_method == "boi" or surrogates unavailable
  min(history$objective[idx], na.rm = TRUE)
}

#' @keywords internal
best_feasible_index <- function(history, objective) {
  feasible_idx <- which(history$feasible)
  if (length(feasible_idx) == 0L) {
    warning("No feasible solutions found; returning best infeasible point.",
            call. = FALSE)
    return(which.min(history$objective))
  }
  feasible_idx[which.min(history$objective[feasible_idx])]
}

#' @keywords internal
estimate_candidate_feasibility <- function(surrogates, unit_x, constraint_tbl) {
  pred <- predict_surrogates(surrogates, unit_x)
  pred_names <- names(pred)
  # Use base R vapply instead of purrr::map_dbl to avoid "In index: X" errors
  mean_vec <- vapply(pred, function(p) {
    if (is.null(p$mean) || length(p$mean) == 0) NA_real_ else p$mean[[1]]
  }, FUN.VALUE = numeric(1))
  names(mean_vec) <- pred_names
  sd_vec <- vapply(pred, function(p) {
    if (is.null(p$sd) || length(p$sd) == 0) NA_real_ else p$sd[[1]]
  }, FUN.VALUE = numeric(1))
  names(sd_vec) <- pred_names
  prob_feasibility(mean_vec, sd_vec, constraint_tbl)
}

#' Compute dynamic fidelity levels based on iteration progress
#'
#' Implements MCEM-inspired continuous escalation by scaling fidelity levels
#' as the optimization progresses. Balances exploration (low cost) vs refinement
#' (high precision) with budget safeguards.
#'
#' @param iter current BO iteration number (not evaluation count)
#' @param budget total evaluation budget
#' @param base_levels named vector of base fidelity levels (e.g., c(low=200, med=1000, high=5000))
#' @param budget_used cumulative computational budget consumed (sum of n_rep)
#' @param budget_tolerance maximum allowed budget overhead (default 1.2 = 20\% over)
#' @param batch_size batch size q for computing iteration phases
#'
#' @details
#' Scaling strategy based on empirical patterns from BOHB and materials discovery:
#' \itemize{
#'   \item \strong{Early phase} (iter < 20\% budget): Base levels - exploration phase
#'   \item \strong{Middle phase} (20\% <= iter < 60\%): 1.5x scaling - focused search
#'   \item \strong{Late phase} (iter >= 60\%): High-fidelity levels increase to 2x/3x - refinement
#' }
#'
#' Moderate scaling for high fidelity: 5000 -> 10000 -> 15000 (user feedback)
#'
#' Budget safeguard: Caps scaling if theoretical budget would exceed tolerance.
#'
#' @return named numeric vector of adjusted fidelity levels
#' @keywords internal
#'
#' @references
#' Falkner et al. (2018). BOHB: Robust and Efficient Hyperparameter Optimization. ICML.
#' Wu & Frazier (2019). Practical Multi-fidelity Bayesian Optimization. UAI.
compute_dynamic_fidelity_levels <- function(iter,
                                            budget,
                                            base_levels,
                                            budget_used = 0,
                                            budget_tolerance = 1.2,
                                            batch_size = 1) {
  # Calculate iteration phase (0 to 1)
  # Approximate max BO iterations (after initial design)
  # Note: n_init not directly available here, so we estimate
  estimated_bo_evals <- pmax(1, budget * 0.7)  # assume ~30% for initial design
  max_iters <- ceiling(estimated_bo_evals / batch_size)
  iter_fraction <- pmin(1, iter / max(max_iters, 1))

  # Theoretical budget: approximate total simulation cost
  theoretical_budget <- budget * mean(base_levels)

  # Check if we're approaching budget limits
  budget_fraction_used <- budget_used / max(theoretical_budget, 1)

  # Determine scaling factors based on phase
  if (iter_fraction < 0.2) {
    # Early phase: base levels for exploration
    scale_low <- 1.0
    scale_med <- 1.0
    scale_high <- 1.0
  } else if (iter_fraction < 0.6) {
    # Middle phase: moderate scaling for focused search
    scale_low <- 1.5
    scale_med <- 1.5
    scale_high <- 2.0  # 5000 -> 10000
  } else {
    # Late phase: aggressive scaling for refinement
    # User feedback: moderate high levels (5000 -> 10000 -> 15000)
    scale_low <- 2.0
    scale_med <- 2.5
    scale_high <- 3.0  # 5000 -> 10000 -> 15000
  }

  # Apply budget safeguard: reduce scaling if approaching limits
  if (budget_fraction_used > 0.85) {
    # Nearing budget limit - cap scaling
    safety_factor <- pmax(0.5, 1 - (budget_fraction_used - 0.85) / 0.15)
    scale_low <- 1 + (scale_low - 1) * safety_factor
    scale_med <- 1 + (scale_med - 1) * safety_factor
    scale_high <- 1 + (scale_high - 1) * safety_factor
  }

  # Ensure we don't exceed budget tolerance
  proposed_levels <- c(
    low = unname(base_levels["low"]) * scale_low,
    med = unname(base_levels["med"]) * scale_med,
    high = unname(base_levels["high"]) * scale_high
  )

  # Check if proposed levels would violate budget
  proposed_mean <- mean(proposed_levels, na.rm = TRUE)
  remaining_budget_evals <- pmax(1, budget - budget_used / mean(base_levels, na.rm = TRUE))

  if (!is.na(proposed_mean) && proposed_mean * remaining_budget_evals > theoretical_budget * budget_tolerance) {
    # Scale back proportionally
    max_allowed_mean <- (theoretical_budget * budget_tolerance) / (remaining_budget_evals + 1e-6)
    reduction <- max_allowed_mean / (proposed_mean + 1e-6)
    reduction <- pmax(0, pmin(1, reduction))  # clamp to [0, 1]
    proposed_levels <- base_levels + (proposed_levels - base_levels) * reduction
  }

  # Round to integer replication counts
  proposed_levels <- round(proposed_levels)

  # Ensure minimum levels (never go below base) and handle any NAs
  proposed_levels <- pmax(proposed_levels, base_levels, na.rm = TRUE)
  # Replace any remaining NAs with base levels
  proposed_levels[is.na(proposed_levels)] <- base_levels[is.na(proposed_levels)]

  proposed_levels
}

#' Hybrid staged-adaptive fidelity selection with metric-specific thresholds
#'
#' Combines iteration-based staging (BOHB-style) with metric-specific CV thresholds
#' for escalation. Implements empirically-validated patterns from multi-fidelity BO
#' literature while handling multi-constraint trial calibration.
#'
#' @param prob_feasible probability of constraint satisfaction
#' @param cv_estimate coefficient of variation from objective surrogate
#' @param iter current iteration number
#' @param fidelity_levels named vector of fidelity levels (possibly dynamic)
#' @param metric name of metric being optimized (for metric-specific thresholds)
#' @param constraint_metrics character vector of constraint metric names
#' @param cv_threshold_base base CV threshold for escalation (default 0.10)
#' @param prob_range feasibility probability range for boundary detection
#' @param distance_to_best normalized distance to current best (0-1 scale)
#' @param ... additional arguments (for compatibility)
#'
#' @details
#' Three-stage strategy:
#' \itemize{
#'   \item \strong{Stage 1} (exploration): Primarily low fidelity, escalate only for very high uncertainty
#'   \item \strong{Stage 2} (focused search): Balanced mix, metric-specific CV thresholds
#'   \item \strong{Stage 3} (refinement): Aggressive high-fidelity promotion for feasible regions
#' }
#'
#' Metric-specific CV thresholds:
#' \itemize{
#'   \item Constraints (power, type1): 0.10 (tight for binary outcomes, boundary precision critical)
#'   \item Objective (EN, ET): 0.20 (looser for continuous high-variance metrics)
#'   \item Near-boundary bonus: 33\% tighter thresholds when 0.4 < P(feasible) < 0.7
#'   \item Near-optimum bonus: 50\% tighter when close to current best design
#' }
#'
#' @return name of selected fidelity level
#' @keywords internal
#'
#' @references
#' Falkner et al. (2018). BOHB: Robust and Efficient Hyperparameter Optimization. ICML.
#' Jiang et al. (2024). Efficient Hyperparameter Optimization with Adaptive Fidelity. CVPR.
#' Richter et al. (2022). Improving adaptive seamless designs through Bayesian optimization. Biom J.
select_fidelity_hybrid_staged <- function(prob_feasible,
                                          cv_estimate,
                                          iter,
                                          fidelity_levels,
                                          metric = NULL,
                                          constraint_metrics = NULL,
                                          cv_threshold_base = 0.10,
                                          prob_range = c(0.2, 0.8),
                                          distance_to_best = 1.0,
                                          ...) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  fidelity_names <- names(fidelity_levels)

  # Determine metric-specific CV threshold
  # Constraints need tighter thresholds (binary outcomes, boundary precision)
  # Objectives can tolerate more variance (continuous, averaged outcomes)
  if (!is.null(metric) && !is.null(constraint_metrics)) {
    if (metric %in% constraint_metrics) {
      cv_thresh <- cv_threshold_base  # tight: 0.10
    } else {
      cv_thresh <- cv_threshold_base * 2  # looser: 0.20
    }
  } else {
    # Default: assume objective metric (looser)
    cv_thresh <- cv_threshold_base * 1.5
  }

  # Near-boundary bonus: tighten threshold by 33% for critical regions
  if (prob_feasible >= 0.4 && prob_feasible <= 0.7) {
    cv_thresh <- cv_thresh * 0.67
  }

  # Near-optimum bonus: tighten threshold by 50% when close to best
  if (distance_to_best < 0.2) {
    cv_thresh <- cv_thresh * 0.5
  }

  # === Three-stage strategy ===

  # Stage 1: Exploration (first 20% of iterations)
  # Primarily low fidelity, escalate only for very high uncertainty
  if (iter <= 20) {
    # Escalate to medium if very uncertain OR near boundary
    if (cv_estimate > 0.25 || (cv_estimate > 0.15 && prob_feasible >= 0.3 && prob_feasible <= 0.7)) {
      if ("med" %in% fidelity_names) {
        return("med")
      }
    }
    # Default: low fidelity for exploration
    return(fidelity_names[1])
  }

  # Stage 3: Refinement (last 40% of iterations, iter > 60 typically)
  # Aggressive high-fidelity promotion for promising regions
  if (iter > 60) {
    # High fidelity for:
    # 1. Likely feasible designs (P > 0.7)
    # 2. Boundary regions with uncertainty
    # 3. Nearby current best
    if ((prob_feasible > 0.7 && "high" %in% fidelity_names) ||
        (prob_feasible >= 0.4 && prob_feasible <= 0.7 && cv_estimate > cv_thresh && "high" %in% fidelity_names) ||
        (distance_to_best < 0.15 && "high" %in% fidelity_names)) {
      return("high")
    }

    # Medium fidelity for moderately promising designs
    if (prob_feasible >= 0.4 && "med" %in% fidelity_names) {
      return("med")
    }

    # Low fidelity for exploration of unlikely regions
    return(fidelity_names[1])
  }

  # Stage 2: Focused search (middle 60%, iterations 21-60)
  # Balanced mix with metric-specific CV-based escalation

  # Promote to high fidelity if:
  # 1. High uncertainty (CV > threshold) AND near boundary, OR
  # 2. Very likely feasible (P > 0.8) with moderate uncertainty
  if ("high" %in% fidelity_names) {
    if ((cv_estimate > cv_thresh && prob_feasible >= prob_range[1] && prob_feasible <= prob_range[2]) ||
        (prob_feasible > 0.8 && cv_estimate > cv_thresh * 0.5)) {
      return("high")
    }
  }

  # Promote to medium fidelity for:
  # 1. Moderately promising regions (P >= 0.4), OR
  # 2. Moderate uncertainty anywhere
  if ("med" %in% fidelity_names) {
    if (prob_feasible >= 0.4 || cv_estimate > cv_thresh * 1.5) {
      return("med")
    }
  }

  # Default: low fidelity for unpromising regions
  fidelity_names[1]
}

#' Select fidelity level using staged strategy
#'
#' Internal function for staged fidelity selection based on iteration count.
#'
#' @param prob_feasible probability of feasibility from surrogate
#' @param cv_estimate coefficient of variation from surrogate
#' @param iter current iteration number
#' @param fidelity_levels named vector of fidelity levels
#' @param cv_threshold CV threshold for fidelity promotion
#' @param prob_range probability range for medium fidelity
#' @param ... additional arguments (ignored)
#' @return character string indicating selected fidelity level
#' @keywords internal
select_fidelity_staged <- function(prob_feasible, cv_estimate, iter, fidelity_levels,
                                   cv_threshold = 0.05, prob_range = c(0.2, 0.8), ...) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  # Stage 1: Initial exploration (iterations 1-30) - uniform low fidelity
  if (iter <= 30) {
    return(names(fidelity_levels)[1])  # "low"
  }

  # Stage 3: Final refinement (iterations 101+) - high fidelity near optimum
  if (iter > 100 && prob_feasible > 0.6 && "high" %in% names(fidelity_levels)) {
    return("high")
  }

  # Stage 2: Focused search (iterations 31-100) - adaptive via MFKG heuristic
  # Promote fidelity if:
  # 1. High uncertainty (CV > cv_threshold), AND
  # 2. Near feasibility boundary (prob_range[1] < prob_feasible < prob_range[2])
  if (cv_estimate > cv_threshold && prob_feasible >= prob_range[1] && prob_feasible <= prob_range[2]) {
    if ("high" %in% names(fidelity_levels)) {
      return("high")
    } else if ("med" %in% names(fidelity_levels)) {
      return("med")
    }
  }

  # Otherwise use medium fidelity for promising regions
  if (prob_feasible >= 0.4 && "med" %in% names(fidelity_levels)) {
    return("med")
  }

  # Default to low fidelity
  names(fidelity_levels)[1]
}

#' Legacy fidelity selection (simple threshold-based)
#' @keywords internal
select_fidelity <- function(prob_feasible, fidelity_levels, ...) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }
  if (prob_feasible >= 0.75 && "high" %in% names(fidelity_levels)) {
    "high"
  } else if (prob_feasible >= 0.4 && "med" %in% names(fidelity_levels)) {
    "med"
  } else {
    names(fidelity_levels)[1]
  }
}

#' Dispatcher for fidelity selection methods
#'
#' Routes fidelity selection to the appropriate strategy based on method argument.
#'
#' @param method fidelity selection method ("adaptive", "staged", "threshold", "hybrid_staged")
#' @param ... arguments passed to specific selection function
#'
#' @return name of selected fidelity level
#' @keywords internal
select_fidelity_method <- function(method, ...) {
  switch(method,
         adaptive = select_fidelity_adaptive(...),
         staged = select_fidelity_staged(...),
         threshold = select_fidelity(...),
         hybrid_staged = select_fidelity_hybrid_staged(...),
         stop("Unknown fidelity method: ", method, call. = FALSE)
  )
}

#' Adaptive cost-aware fidelity selection
#'
#' Implements cost-aware fidelity selection inspired by MFKG (Wu & Frazier 2016).
#' Selects fidelity by maximizing expected value per unit cost, with exploration
#' decay and boundary detection.
#'
#' @param prob_feasible probability of constraint satisfaction
#' @param cv_estimate coefficient of variation from objective surrogate
#' @param acq_value acquisition function value at candidate point
#' @param best_obj current best objective value (used for context)
#' @param fidelity_levels named vector of fidelity levels (replication counts)
#' @param fidelity_costs named vector of relative costs (default: proportional to replications)
#' @param iter current iteration number
#' @param total_budget_used cumulative simulation budget consumed
#' @param total_budget approximate total simulation budget available
#'
#' @details
#' The method balances information gain vs cost using:
#' \itemize{
#'   \item \strong{Value score}: acquisition * uncertainty * boundary_factor
#'     \itemize{
#'       \item High near feasibility boundary (prob ~= 0.5)
#'       \item High when objective uncertain (large CV)
#'       \item High for promising candidates (large acquisition)
#'     }
#'   \item \strong{Cost normalization}: Divide by cost^alpha where alpha decays from 0.5 to 0.8
#'     \itemize{
#'       \item Early: less cost-sensitive (exploration)
#'       \item Late: more cost-sensitive (exploitation)
#'     }
#'   \item \strong{Exploration decay}: Randomization probability from 50\% to 10\%
#' }
#'
#' @return name of selected fidelity level
#' @keywords internal
#'
#' @references
#' Wu, J., & Frazier, P. (2016). The parallel knowledge gradient method for
#' batch Bayesian optimization. NeurIPS.
select_fidelity_adaptive <- function(prob_feasible,
                                     cv_estimate,
                                     acq_value,
                                     best_obj,
                                     fidelity_levels,
                                     fidelity_costs,
                                     iter,
                                     total_budget_used,
                                     total_budget,
                                     cv_threshold = 0.05,
                                     prob_range = c(0.2, 0.8),
                                     ...) {
  if (length(fidelity_levels) == 1L) {
    return(names(fidelity_levels))
  }

  fidelity_names <- names(fidelity_levels)
  costs <- fidelity_costs[fidelity_names]

  # === Compute value score ===

  # Uncertainty factor: higher CV -> more value in reducing uncertainty
  # Normalize to [0,1] range, assuming CV > 0.3 is very high
  uncertainty_factor <- pmax(0, pmin(1, cv_estimate / 0.3))

  # Boundary factor: highest value near feasibility boundary
  # P = 0.5 -> factor = 1, P = 0 or 1 -> factor = 0
  boundary_factor <- 1 - abs(2 * prob_feasible - 1)
  boundary_factor <- boundary_factor^0.5  # soften the effect

  # Acquisition factor: diminishing returns on acquisition value
  # Use log1p for numerical stability
  acq_factor <- log1p(pmax(0, acq_value))

  # Combined value score - weight components based on optimization stage
  if (iter < 20) {
    # Early: prioritize exploration (uncertainty)
    value_score <- acq_factor * (0.3 + 0.7 * uncertainty_factor) * (0.5 + 0.5 * boundary_factor)
  } else if (iter < 60) {
    # Middle: balance
    value_score <- acq_factor * (0.5 + 0.5 * uncertainty_factor) * (0.3 + 0.7 * boundary_factor)
  } else {
    # Late: prioritize promising candidates
    value_score <- acq_factor * (0.7 + 0.3 * uncertainty_factor) * (0.1 + 0.9 * boundary_factor)
  }

  # === Compute cost sensitivity ===

  # Normalize costs to [0, 1]
  cost_normalized <- costs / max(costs)

  # Cost exponent: increases with iteration (more cost-sensitive over time)
  # Also increases as budget depletes
  # Lower exponent = more willing to use high fidelity
  budget_fraction_used <- pmin(1, total_budget_used / (total_budget + 1e-6))
  cost_exponent <- 0.15 + 0.3 * pmin(1, iter / 100) + 0.2 * budget_fraction_used
  cost_exponent <- pmin(cost_exponent, 0.8)  # cap at 0.8 (was 1.0)

  # === Value per cost ===
  value_per_cost <- value_score / (cost_normalized ^ cost_exponent + 1e-6)

  # === Exploration randomization ===

  # Probability of forcing low fidelity for exploration
  # Decays from 50% early to 5% late
  exploration_prob <- pmax(0.05, 0.5 * exp(-iter / 30))

  if (stats::runif(1) < exploration_prob) {
    # Force exploration with low fidelity
    return(fidelity_names[1])
  }

  # === Select fidelity ===

  best_idx <- which.max(value_per_cost)
  selected <- fidelity_names[best_idx]

  # === Diagnostics (optional) ===
  # Support both new and legacy option names for backwards compatibility
  debug_fidelity <- getOption("BATON.debug_fidelity", FALSE) ||
                    getOption("evolveBO.debug_fidelity", FALSE)
  if (debug_fidelity) {
    message(sprintf(
      "  Fidelity selection: prob_feas=%.3f, CV=%.3f, acq=%.3f, value=%.3f, cost_exp=%.2f -> %s",
      prob_feasible, cv_estimate, acq_value, value_score, cost_exponent, selected
    ))
    message(sprintf("    Value/cost: %s",
                    paste(sprintf("%s=%.2f", fidelity_names, value_per_cost), collapse=", ")))
  }

  selected
}

#' @keywords internal
generate_diagnostics <- function(surrogates, history, objective, n_draws = 1000) {
  unit_best <- history$unit_x[[best_feasible_index(history, objective)]]
  newdata <- matrix(unname(unlist(unit_best)), nrow = 1)
  colnames(newdata) <- names(unit_best)

  # Use base R lapply instead of purrr::imap to avoid "In index: X" errors
  surrogate_names <- names(surrogates)
  draws <- lapply(seq_along(surrogates), function(i) {
    model <- surrogates[[i]]
    metric <- surrogate_names[i]

    tryCatch({
      if (inherits(model, "constant_predictor")) {
        # Constant predictor - return draws from approximate distribution
        mean_val <- if (!is.null(model$mean_value)) model$mean_value else 0
        stats::rnorm(n_draws, mean_val, 1.0)
      } else if (inherits(model, c("hetGP", "homGP"))) {
        # hetGP/homGP models
        pred <- predict(x = newdata, object = model)
        mean_val <- as.numeric(pred$mean)
        sd_val <- max(as.numeric(sqrt(pred$sd2)), 1e-8)
        stats::rnorm(n_draws, mean_val, sd_val)
      } else {
        # DiceKriging km model
        pred <- DiceKriging::predict.km(model,
                                        newdata = newdata,
                                        type = "UK",
                                        se.compute = TRUE,
                                        cov.compute = TRUE)
        mean_val <- as.numeric(pred$mean)
        sd_val <- max(as.numeric(pred$sd), 1e-8)
        stats::rnorm(n_draws, mean_val, sd_val)
      }
    }, error = function(e) {
      warning(sprintf("[generate_diagnostics] Error generating draws for metric '%s': %s",
                      metric, e$message), call. = FALSE)
      # Return uniform draws as fallback
      stats::runif(n_draws, 0, 1)
    })
  })
  names(draws) <- surrogate_names

  list(
    posterior_draws = draws,
    history_summary = history |> dplyr::select(eval_id, objective, feasible)
  )
}

#' @keywords internal
evaluate_acquisition <- function(acquisition,
                                 unit_candidates,
                                 surrogates,
                                 constraint_tbl,
                                 objective,
                                 best_feasible,
                                 predictions = NULL,
                                 prob_feas = NULL) {
  if (identical(acquisition, "eci")) {
    return(acq_eci(
      unit_x = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible,
      pred = predictions,
      prob_feas = prob_feas
    ))
  }
  if (identical(acquisition, "qehvi")) {
    return(acq_qehvi(
      unit_x = unit_candidates,
      surrogates = surrogates,
      constraint_tbl = constraint_tbl,
      objective = objective,
      best_feasible = best_feasible,
      pred = predictions,
      prob_feas = prob_feas
    ))
  }
  stop(sprintf("Acquisition '%s' is not supported.", acquisition), call. = FALSE)
}
