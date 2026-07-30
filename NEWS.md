# BATON Changelog

## BATON 0.7.1 (2026-07-30)

Documentation release; no code changes beyond two small guards shipped
with the earlier README overhaul.

- **BATON-introduction rewritten as a self-contained, discursive Getting
  Started guide** (1000+ lines): the evaluator contract and minimal
  template with per-line narration, a fully narrated worked example with
  real captured output, per-knob settings guidance with failure modes,
  the quick/balanced/thorough run profiles, live demonstrations of
  infeasibility, verification, and integer rounding, and a trust
  checklist. No longer defers to the GitHub README.
- **New vignette `BATON-methods`**: how the optimizer works. GP
  surrogates under Monte Carlo noise, expected constrained improvement,
  local-penalization batching, multi-fidelity budgeting, the iteration
  loop, and the seeding architecture behind core-count-independent
  reproducibility, with small executed demonstrations. Paraphrases and
  cites the companion manuscript; where manuscript and code differ, the
  vignette describes the code.
- **New vignette `BATON-case-study`**: an end-to-end calibration of a
  group-sequential design carried from problem statement through
  quick-profile exploration, a precomputed balanced run, the trust
  pipeline (including an instructive multi-seed FAIL and its remedy),
  sensitivity analysis, benchmarking against random and grid search,
  constraint reliability, multi-philosophy calibration with warm-start
  donors, and long-run service controls. Heavy results ship precomputed
  (12 KB) with the generating script included.

## BATON 0.7.0 (2026-07-24)

Service-readiness release (Phase D of the 2026-07 review plan): run controls
for calibrations driven by a service or scheduler, RNG hygiene, and a lighter
dependency footprint. Numerical output is unchanged: the golden-master
calibration snapshot is byte-identical to v0.6.0.

### Service controls

- **`status` field.** Every `BATON_fit` records why the run ended:
  `"budget_exhausted"`, `"early_stopped"`, `"acq_flatline"`, `"walltime"`,
  `"cancelled"`, or `"errored"`.
- **`max_walltime_s` graceful stop.** Optional walltime cap checked between
  initial-design evaluation chunks and at the top of each BO iteration
  (advisory, never preemptive: a running chunk or batch is not interrupted).
  On trip the run returns a partial but internally consistent fit with
  `status = "walltime"`, and Stage 4 multi-seed verification is skipped even
  when `multi_seed_verify = TRUE` (it could dwarf the cap itself). With a
  cap set, initial-design evaluation proceeds in worker-count chunks so the
  check also fires during the init phase; with no cap, evaluation order and
  seeds are bit-identical to v0.6.0.
- **`callback` cooperative cancellation.** Optional hook invoked at the end
  of each BO iteration with `iter`, `eval_count`, `best_objective` (the ECI
  incumbent), `best_observed` (best observed feasible objective including
  the just-completed batch), and `elapsed_s`. Return `TRUE` to continue;
  any other value cancels after the current iteration with
  `status = "cancelled"` and no completed work discarded. A callback error
  warns and continues: a monitoring hook must not kill a paid run.
- **`checkpoint_fun` resumable snapshots.** Every `checkpoint_every`
  completed BO iterations the hook receives a slim `BATON_fit` snapshot
  (`status = "checkpoint"`, full evaluation history, no GP objects). Resume
  a saved checkpoint with `bo_calibrate(..., initial_history =
  ckpt$history, n_init = nrow(ckpt$history))`.
- **`on_error = "return_partial"`.** An error inside an evaluation round or
  BO iteration (simulator failure, surrogate fitting failure) ends the run
  gracefully: the fit comes back with `status = "errored"`, the condition
  message in `error_message`, and every completed prior round preserved. A
  failure inside a round loses only that round's evaluations. The default
  `on_error = "stop"` propagates errors exactly as before.
- **`slim = TRUE` lightweight return.** Drops `surrogates`, `diagnostics`,
  and `multi_seed_runs` (km/hetGP objects carry closures and environments
  that inflate `saveRDS` payloads) while keeping the decision content:
  `history`, `best_theta`, `best_objective`, `policies`,
  `multi_seed_summary`, `verdict`, and `status`. The final surrogate refit
  is skipped entirely (nothing else consumes it), saving one full GP fit
  per run. Slim fits stay resumable via `initial_history = fit$history`;
  `sa_sobol()`, `sensitivity_diagnostics()`, and `extract_case_study()`
  need a non-slim fit.

### RNG hygiene

- `bo_calibrate()` no longer clobbers the caller's RNG stream. With an
  explicit `seed`, `.Random.seed` is restored on exit (success or error)
  and the caller's stream is untouched. With `seed = NULL`, the stream
  advances by exactly one draw (the internal seed draw), so back-to-back
  `seed = NULL` runs remain independent replicates.

### Dependency slimming

- `ggplot2` and `randtoolbox` moved to Suggests; `purrr`, `tidyr`, and
  `rlang` dropped entirely (all call sites replaced with base R).
  Installation is lighter: the calibration core needs only DiceKriging,
  hetGP, lhs, dplyr, and tibble. Plotting (`plot_*`) requires ggplot2 to be
  installed and fails with an actionable install hint naming the calling
  function when it is absent. The init-stopping GP test grid uses Sobol
  points when randtoolbox is installed and falls back to random points
  (with a message) otherwise; the LHS initial design is unaffected.

## BATON 0.6.0 (2026-07-23)

Completes the 2026-07 review plan: the Phase B contract changes and the
Phase C performance work. All performance changes are gated by the
golden-master calibration snapshot (byte-identical results before and after)
and a core-count-independence test.

### Performance (Phase C)

- **Parallel simulator evaluation.** The initial design, warm-start top-up,
  and BO batches evaluate through `evaluate_points()`, opt-in parallel via
  `options(BATON.cores = k)` on Unix. Every simulator call runs under
  `run_seeded()`, so results are identical serial vs parallel even for
  simulators drawing from the ambient RNG stream. History rows are appended
  in one bind per round instead of one per evaluation.
- **Matrix candidate scoring.** `predict_surrogates()` accepts the candidate
  matrix directly and the acquisition hot loop is list-free end to end
  (measured 58% faster scoring of a 2000-candidate pool). Invalid candidates
  now error instead of being silently dropped, which had misaligned
  predictions with candidate indices.
- **Matrix history aggregation.** Per-theta tibble construction in surrogate
  fitting is replaced by a single `rowsum`-based pass with exact
  `mean(na.rm = TRUE)` semantics (all-NA groups stay NaN).
- **Batched sensitivity predicts.** `sa_sobol()`, `sa_gradients()`,
  `cov_effects()`, and gradient sampling issue one predict per metric
  instead of one per (point, parameter) -- measured 36x faster gradient
  covariance estimation.
- **hybrid_staged proximity bonus now fires.** The distance-to-best lookup
  matched rows by float equality against the posterior-mean incumbent and
  never matched; it now finds the best feasible row by index (preferring
  high-fidelity rows, rescaled under current bounds, bounds-ordered).

### Fixes

- Ablation policies with a subset of fidelity levels (e.g.
  `low_only = c(low = 500)`) are normalized onto the low/med/high contract
  instead of erroring.
- Histories without a `variance` column keep their homoskedastic km fit
  (a regression during the aggregation rewrite would have silently
  downgraded them to constant predictors).
- Warm-start top-up replenishes dedup rejections, detects true design-space
  exhaustion for integer grids (including `round()` reach outside
  non-integral bounds), and never underfills `n_init` while distinct
  designs remain.

### Contract changes (Phase B of the 2026-07 review plan)

- **Optional `n_rep` simulator argument.** `bo_calibrate()` now passes the
  requested replication count (`fidelity_levels[[fidelity]]`, including
  hybrid_staged escalated values) to simulators that declare an explicit
  `n_rep` formal. Accepting `...` alone does not opt in, so wrappers that
  forward dots to legacy simulators keep working; simulators without the
  argument are called exactly as before. `fix_parameters()` wrappers pass
  `n_rep` through when both caller and wrapped simulator support it.
- **Warm-start seeds augment the LHS initial design.** `warmstart_from` /
  `initial_history` rows no longer replace the initial design: when fewer than
  `n_init` rows are provided, fresh LHS points (deduplicated against donor
  designs and against integer-coerced duplicates, with re-seeded replenishment
  rounds) top up the design to `n_init`. Histories with `n_init` or more rows
  are used as-is, preserving multi-stage warm-starting with narrowed bounds.

## BATON 0.5.0 (2026-07-23)

Correctness and performance release from a full-package code review. All
changes are covered by new regression tests; a golden-master calibration
snapshot (with a core-count-independence check) guards the numerical output.

### Performance

- **Warm-start now functions.** `parinit` was passed inside `km()`'s `control`
  list, which DiceKriging ignores, so the advertised v0.3.0 warm-start was a
  no-op and every iteration paid a full cold-start MLE. It is now a top-level
  `km()` argument (and `mleHetGP`/`mleHomGP` are seeded from the previous
  `theta`), giving ~30-50% faster surrogate fitting (measured ~34% at n=120,
  d=6).
- **Opt-in parallel GP fitting.** The `m` per-metric fits can run concurrently
  via `options(BATON.cores = k)` on Unix. Each fit runs under a deterministic
  seed (`run_seeded()`), so results are identical regardless of core count -
  reproducibility does not depend on the host.
- **Objective-only incumbent.** `best_feasible_objective()` predicted all `m`
  surrogates and used only the objective; it now predicts just the objective.

### Correctness

- **Local penalization sign corrected.** The batch-diversity penalty grew with
  distance, so `q > 1` batches clustered around the first pick (opposite of
  Gonzalez et al. 2016); batches now spread out.
- **Expected constraint violation** now uses the correct Gaussian expected
  shortfall `sd*(z*Phi(z)+phi(z))` instead of a form that collapsed to
  `sd*phi(z)` (which peaked at the boundary and vanished for deeply infeasible
  points).
- **Welford utilities:** `welford_mean_var()` matches metrics by name (was
  positional); `pool_welford_results()` handles `n=1` chunks (was NA-poisoning).
- **`theta_to_id()`** uses fixed-notation per-coordinate formatting, stable at
  bound edges and across session options (duplicate detection no longer breaks).
- **Multi-seed verification** accepts the documented named-vector simulator
  contract (previously crashed with `$ operator is invalid for atomic vectors`
  after the full budget was spent).
- **Incumbent fidelity fallback:** returns the any-fidelity feasible best when
  no feasible high-fidelity evaluation exists, instead of `Inf`.
- **Early stopping:** the improvement check no longer compares across the
  infeasible-to-feasible transition, and the acquisition-flatline threshold is
  now relative to the incumbent (new `early_stop$acq_rel`, default 1e-3).
- **`bo_calibrate_philosophies()`** isolates per-philosophy failures (one
  failure no longer discards completed fits; errored fits get
  `status = "errored"` / verdict `"ERRORED"`), and donor warm-start attaches the
  operating characteristics of `best_theta`'s own evaluation.
- **`warmstart_from` path fixed.** `build_initial_history_from_warmstart()`
  initialised a 0-row data frame and threw `replacement has 1 row, data has 0`
  for every warmstart call; it now seeds the correct row count.
- **`is_feasible()`** is robust to atomic metric vectors with missing names
  (was an opaque `subscript out of bounds`; now surfaces the clear
  name-contract error).
- **`parse_constraints()`** validates input (rejects unnamed lists, non-numeric
  or non-finite thresholds, malformed entries) instead of silently producing
  `NA` feasibility.
- Removed the always-erroring `compute_global_sensitivity()` export (it called
  a function that does not exist); use `sa_sobol()` or
  `compute_parameter_importance()`.
- Benchmark grid/random/heuristic arms no longer crash on the documented
  named-vector fidelity input.
- Unknown (`NA`) noise is imputed with the maximum observed noise, not the
  minimum (unknown-variance points are now trusted least, not most).

## BATON 0.4.0 (2026-05-09)

### Hardening release: cross-philosophy warmstart + multi-seed verification

This release adds two safeguards to `bo_calibrate()` and a new batch
wrapper `bo_calibrate_philosophies()`. Both safeguards are motivated by
empirical failure modes documented in the bounds-widening campaign for
the JASA 2026 manuscript on "Constrained Bayesian Optimization for
Calibration of Bayesian Adaptive Clinical Trials" (Web Appendix G.5).

### New features

- **`warmstart_from` parameter (`bo_calibrate()`).** Accepts a list of
  previously calibrated designs (one-row data frames with parameter
  columns matching `names(bounds)`, paths to `*_best_design.csv` files,
  or lists with a `theta` element). Inject these designs as initial
  Stage 0 seeds alongside the LHS sample, helping the GP surrogate
  escape local minima that an LHS-only Stage 0 might miss. Designs
  outside the current bounds are filtered with a warning. Default
  `NULL` (LHS-only initial design pool, backward compatible).

  *Why this matters:* In the manuscript campaign, the multi-arm
  Optimal philosophy calibration settled at a local minimum with
  $E_0[N] = 67.19$ that Fleming's calibrated point at $E_0[N] = 66.18$
  Pareto-dominated under Optimal's own loss function. Seeding Optimal's
  BO Stage 0 with Fleming's top-3 high-fidelity designs recovered an
  Optimal design at $E_0[N] = 65.66$, strictly better than both the
  original Optimal calibration and Fleming under Optimal's loss. This
  resolution is now built into the package.

- **`multi_seed_verify` / `multi_seed_n` / `multi_seed_reps` /
  `multi_seed_strict` parameters (`bo_calibrate()`).** When
  `multi_seed_verify = TRUE`, after Stage 3 selects the calibrated
  design, BATON re-evaluates that design at `multi_seed_n` independent
  seeds at `multi_seed_reps` fidelity (defaults: 5 seeds at
  `fidelity_levels[["high"]]` reps each). The returned `BATON_fit`
  object gains a `multi_seed_summary` element (mean/SD of objective and
  each constraint metric across seeds), a `multi_seed_runs` element
  (per-seed evaluation rows), and a `verdict` field
  (`"MULTI_SEED_PASS"`, `"MULTI_SEED_FAIL"`, or `"MULTI_SEED_WARN"`).
  When `multi_seed_strict = TRUE` (default) and any of the seeds
  produces an infeasible operating characteristic, the verdict is
  `"MULTI_SEED_FAIL"` and downstream code can gate on
  `fit$multi_seed_summary$strict_feas == 1` before adopting the
  design. Default `multi_seed_verify = FALSE` for backward
  compatibility — opt-in.

  *Why this matters:* In the campaign, three SYMM cohortA bidirectional
  calibrated designs had Stage 3 stored type-I error below the 0.10
  regulatory floor but multi-seed type-I error in [0.122, 0.170] —
  stored-vs-actual gaps of approximately 0.075, more than 20 Monte
  Carlo standard errors. Without multi-seed verification, those
  designs would have shipped type-I-inflated. Multi-seed verification
  is now built into the package.

- **`bo_calibrate_philosophies()` batch wrapper.** New function that
  calibrates multiple design philosophies for a single scenario in
  dependency order. Donor philosophies (Fleming, Minimax, Admissible
  variants) are calibrated first without warmstart; recipient
  philosophies (Optimal, Alt-Optimal) are calibrated second with
  `warmstart_from` populated by the donors' calibrated designs. Every
  calibration runs with `multi_seed_verify = TRUE` by default. The
  function returns a manifest documenting per-philosophy verdicts and
  warmstart provenance. See `?bo_calibrate_philosophies`.

  *Why this matters:* The two failure modes above are most relevant
  when calibrating a panel of philosophies (e.g., for Pareto-frontier
  comparison in a methodology paper). This wrapper makes the safe
  workflow the default for that use case, rather than requiring users
  to assemble it from low-level pieces.

### Backward compatibility

All v0.4.0 changes are additive and opt-in:

- `multi_seed_verify` defaults to `FALSE`. Existing pipelines see no
  behavior change.
- `warmstart_from` defaults to `NULL`. Existing pipelines see no
  behavior change.
- `bo_calibrate_philosophies()` is new; existing calls to
  `bo_calibrate()` are unaffected.
- The existing `initial_history` parameter (introduced for v0.4.0
  scaffolding in v0.3.0) continues to work; `warmstart_from` is a
  higher-level convenience layer that converts user-friendly design
  specifications to `initial_history` internally. Specifying both
  parameters in the same call is an error.

### Mental model for the four "stages"

The three existing stages (Stage 1, 2, 3) are **fidelity tiers** within
a single BO run: low-fidelity exploration, medium-fidelity refinement,
high-fidelity verification of the BO surrogate's belief about the
calibrated design at one seed. The new Stage 4 (multi-seed verification)
addresses an **orthogonal axis** — seed-dependent variability of the
calibrated design's true operating characteristics. Stage 1/2/3 ask
"how confident is BO that this design is good?"; Stage 4 asks "given
that BO settled on this design, how confident are we that the design's
true OCs match what BO estimated?".

Both questions matter; they are answered by different mechanisms.
Stage 1/2/3 fidelity tier hardening alone does not catch
seed-dependent OC variability at the calibrated point. Stage 4
multi-seed verification does.

### Recommended workflow for production calibrations

```r
fit <- BATON::bo_calibrate(
  sim_fun = your_simulator,
  bounds = your_bounds,
  objective = "weighted_EN_N",
  constraints = list(power = c("ge", 0.80), type1 = c("le", 0.10)),
  budget = 300,
  fidelity_levels = c(low = 3000, med = 5000, high = 10000),
  multi_seed_verify = TRUE,
  multi_seed_n = 5,
  multi_seed_strict = TRUE
)

# NOTE: `objective` and every constraint name must be the name of a metric
# `your_simulator` returns. There is no built-in weighted-loss objective, so
# compute the weighted loss (here `weighted_EN_N`) inside the simulator and
# return it as a named metric (see the README "Design Philosophies" section).

if (fit$verdict == "MULTI_SEED_PASS") {
  adopt_design(fit$best_theta)
} else {
  # Recalibrate, possibly with warmstart_from a sibling philosophy's design
  message("Multi-seed verification failed; investigate before adoption.")
  print(fit$multi_seed_summary)
}
```

For multi-philosophy calibrations:

```r
batch <- BATON::bo_calibrate_philosophies(
  scenario_id = "cohortA_strong_effect",
  sim_fun = your_simulator,
  bounds = your_bounds,
  philosophies = list(
    Fleming   = list(objective = "weighted_EN_N",
                     constraints = list(...)),
    Minimax   = list(objective = "total_n",
                     constraints = list(...)),
    Admissible = list(objective = "weighted_EN_N", constraints = list(...)),
    Optimal   = list(objective = "EN_null",
                     constraints = list(...)),
    `Alt-Optimal` = list(objective = "EN_alt",
                         constraints = list(...))
  ),
  budget = 300
)

batch$manifest  # per-philosophy verdicts + warmstart provenance
batch$fits$Optimal$best_theta  # individual fits
```

---

## BATON 0.3.0 (prior release)

- Heteroskedastic GP surrogates via hetGP fallback
- Adaptive multi-fidelity simulation (Stage 1/2/3 fidelity tiers)
- Constraint-aware acquisition (ECI)
- Sensitivity diagnostics
- Default batch size `q = 2` (changed from `q = 8`)
- `initial_history` parameter (low-level warm-start API; superseded for
  user-facing use by v0.4.0's `warmstart_from`)
