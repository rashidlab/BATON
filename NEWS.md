# BATON Changelog

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
