# BATON Review Fixes and Performance Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix the confirmed correctness bugs from the 2026-07-22 review, lock current
numerical behavior with regression tests, then land the performance work (parallelism,
warm-start, incremental aggregation) that makes BATON fit for cloud hosting.

**Architecture:** Three sequenced phases. Phase A fixes correctness bugs test-first (each
bug gets a failing test that encodes correct behavior, then the fix). Phase B is two
public-contract changes that need a PI decision before coding. Phase C is performance:
the behavior-preserving speed work (parallelism, matrix scoring, aggregation, append) runs
FIRST against a characterization-test harness that pins the current calibration output, so a
"pure speedup" cannot silently change results. The Phase A fixes that intentionally change
the search trajectory (warm-start, penalization, expected violation) come after and
re-baseline that harness. See the C0 execution-order note.

**Tech Stack:** R 4.6.0, testthat 3e (3.3.2), DiceKriging 1.6.1, hetGP 1.1.9,
`parallel` (base) for multicore. Test command:
`Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/<file>.R")'`
Full check: `Rscript -e 'devtools::test()'`.

**Branch:** `fix/review-2026-07` (off dev main). Commit after every task.

**Reference:** `CODE_REVIEW_2026-07-22.md` (item numbers below map to its Section 1/2).

---

## Cross-cutting rules

- **TDD**: every correctness task writes the test first, watches it fail, then fixes.
- **Regression safety for Phase C**: before any performance refactor, Task C0 builds a
  golden-master test that runs a small fixed-seed calibration and snapshots the history
  and best_theta. Every Phase C task must leave that snapshot byte-identical (a speedup
  must not change answers). Parallelism tasks additionally assert seed-determinism.
- **Commits**: one per task, message `fix:`/`test:`/`perf:`/`refactor:` as appropriate,
  ending with the Co-Authored-By trailer.
- **No public-repo touch**: all work stays in the dev repo; releases are snapshotted later
  via `scripts/release-to-public.sh`.

---

# PHASE A — Correctness fixes (test-first)

### Task A1: Warm-start passes `parinit` as a top-level km argument

**Files:**
- Modify: `R/surrogates.R:374-378` and `R/surrogates.R:386-390`
- Test: `tests/testthat/test-warmstart.R` (add case)

**Step 1 — Write the failing test.** Fit a km on a small design, extract its hyperparams,
then assert that passing them as warm-start actually seeds `model@parinit` (proving the
argument is honored, not ignored):

```r
test_that("km warm-start actually seeds parinit (not swallowed by control)", {
  set.seed(1)
  X <- matrix(runif(30), ncol = 3); y <- rowSums(X) + rnorm(10, sd = 0.01)
  m1 <- DiceKriging::km(design = as.data.frame(X), response = y,
                        covtype = "matern5_2", nugget.estim = TRUE,
                        control = list(trace = FALSE))
  hp <- as.numeric(m1@covariance@range.val)
  # BATON's fitter must forward hp so the fitted model's parinit equals hp
  aggr <- list(value = y)
  m2 <- fit_dicekriging_surrogate(X_unique = as.data.frame(X), aggr = aggr,
                                  noise_vec = NULL, covtype = "matern5_2",
                                  prev_model = m1)
  expect_equal(as.numeric(m2@parinit), hp, tolerance = 1e-8)
})
```

(Confirm the real signature of `fit_dicekriging_surrogate` first; adapt the call. If it is
not separately callable, test through `fit_surrogates` with a `prev_surrogates` arg.)

**Step 2 — Run, expect FAIL** ("parinit is ... numeric(0)" or a random value != hp).

**Step 3 — Fix.** In both km calls move `parinit` out of `control`:

```r
DiceKriging::km(
  design = X_unique,
  response = aggr$value,
  covtype = covtype,
  nugget = nugget,          # or noise.var = noise_vec in the second call
  nugget.estim = FALSE,
  parinit = parinit,
  control = list(trace = FALSE)
)
```

Guard: only pass `parinit` when `!is.null(parinit) && length(parinit) == ncol(X_unique)`;
otherwise omit it (km cold-starts). Build the arg list conditionally with `do.call`.

**Step 4 — Run, expect PASS.**

**Step 5 — Commit** `fix: forward parinit as top-level km argument (warm-start was a no-op)`.

---

### Task A2: hetGP warm-start uses previous theta

**Files:**
- Modify: `R/surrogates.R:126` (`fit_hetgp_surrogate` body, the `mleHetGP`/`mleHomGP` calls)
- Test: `tests/testthat/test-warmstart.R`

**Step 1 — Failing test:** assert that when `prev_model` is supplied, the hetGP fit is
seeded (e.g., fewer optimizer iterations, or `init$theta` reflected). Simplest observable:
patch the call to pass `init = list(theta = prev_model$theta)` and assert the fit converges
to the same lengthscale from the warm seed on identical data (`expect_equal(fit$theta, ...)`).
If no robust observable exists, assert `prev_model` is referenced by fitting twice and
checking the second fit's `msg`/iteration count is <= cold. Keep the test tolerance-based.

**Step 2 — Run, expect FAIL.**

**Step 3 — Fix.** In the `mleHetGP`/`mleHomGP` calls add
`init = if (!is.null(prev_model)) list(theta = prev_model$theta) else NULL`. Wrap in
`tryCatch` so a stale/incompatible seed falls back to a cold fit rather than erroring.

**Step 4 — PASS. Step 5 — Commit** `fix: warm-start hetGP fits from previous theta`.

---

### Task A3: Correct local-penalization sign (batch diversity)

**Files:**
- Modify: `R/acquisition.R:209-215`
- Test: `tests/testthat/test-phase1-improvements.R`

**Step 1 — Failing test** encoding the intended behavior (near duplicates suppressed, far
points kept):

```r
test_that("local penalization spreads the batch (suppresses near-duplicates)", {
  # scores decreasing left to right; first pick is candidate 1 (x=0.00)
  cand <- matrix(c(0.00, 0.01, 0.99), ncol = 1)
  scores <- c(1.0, 0.95, 0.80)
  sel <- select_batch_local_penalization(cand, scores, q = 2, lipschitz = 10)
  # must pick the FAR point (index 3), not the near-duplicate (index 2)
  expect_equal(sort(sel), c(1, 3))
})
```

**Step 2 — Run, expect FAIL** (current code returns c(1, 2)).

**Step 3 — Fix.** Replace the inverted penalty with the González hammer that suppresses a
ball *around* the selected point (penalty largest at r=0, decaying with distance):

```r
# Suppress candidates NEAR the selected point: penalty is large at r=0,
# fades to 0 beyond radius acq_best / L.
radius <- penalized_scores[best_idx] / lipschitz
penalty <- pmax(0, penalized_scores[best_idx] - lipschitz * distances)
penalized_scores <- penalized_scores - penalty
penalized_scores[best_idx] <- -Inf
```

Update the now-wrong comment. Remove the 1.5x lipschitz inflation at the call site
(`bo_calibrate.R:714`) since it no longer means "stronger diversity" under the corrected
sign — verify direction with a quick reasoning note in the commit body.

**Step 4 — PASS. Step 5 — Commit** `fix: correct inverted local-penalization sign (batches were clustering)`.

---

### Task A4: `compute_expected_violation` computes the real expected violation

**Files:**
- Modify: `R/acquisition.R:144-157`
- Test: `tests/testthat/test-phase1-improvements.R`

**Step 1 — Failing test** against the closed form E[max(0, t - Y)] = sd*(z*Φ(z)+φ(z)):

```r
test_that("expected violation matches closed form and grows when deeply infeasible", {
  mu <- 0.2; sd <- 0.05; thr <- 0.8  # ge-constraint, far infeasible
  z <- (thr - mu) / sd
  expect_true(compute_expected_violation(<pred>, <constraints>) >
              compute_expected_violation(<pred at mu=0.79>, <constraints>))
  # boundary point must have SMALL expected violation, deep point LARGE
})
```

(Build the small `pred`/`constraints` fixtures the function expects; assert monotonic
increase in violation as mu moves further below the threshold.)

**Step 2 — Run, expect FAIL** (current symmetric form peaks at the boundary).

**Step 3 — Fix.** Both branches:

```r
# ge: violation when metric < threshold; z = (threshold - mu)/sd
z <- (threshold - mu_vec) / sd_vec
violations <- violations + sd_vec * (z * stats::pnorm(z) + stats::dnorm(z))
```

```r
# le: violation when metric > threshold; z = (mu - threshold)/sd
z <- (mu_vec - threshold) / sd_vec
violations <- violations + sd_vec * (z * stats::pnorm(z) + stats::dnorm(z))
```

**Step 4 — PASS. Step 5 — Commit** `fix: compute true expected constraint violation (was sd*dnorm)`.

---

### Task A5: Welford name-matching and n=1 pooling

**Files:**
- Modify: `R/welford.R` (`welford_mean_var` sample loop ~94-103; `pool_welford_results` ~194-200)
- Test: `tests/testthat/test-welford.R` (NEW — currently zero coverage)

**Step 1 — Failing tests:**

```r
test_that("welford_mean_var reorders samples by name", {
  i <- 0
  sf <- function() { i <<- i + 1; if (i %% 2) c(a=1, b=10) else c(b=10, a=1) }
  res <- welford_mean_var(sf, n_samples = 100)
  expect_equal(res$mean[["a"]], 1); expect_equal(res$mean[["b"]], 10)
})
test_that("pool_welford_results handles n=1 chunks without NA", {
  big <- welford_mean_var(function() rnorm(1, 5, 1) |> setNames("m"), 100)
  one <- list(mean = c(m = 4.0), variance = c(m = NA), n = 1)
  pooled <- pool_welford_results(list(big, one))
  expect_false(is.na(pooled$variance[["m"]]))
})
test_that("pooled mean/variance match direct computation", {
  set.seed(3); x <- rnorm(200)
  chunks <- lapply(split(x, rep(1:4, each=50)), function(v)
    welford_mean_var(local({j<-0; function() {j<<-j+1; c(m=v[j])}}), length(v)))
  pooled <- pool_welford_results(chunks)
  expect_equal(pooled$mean[["m"]], mean(x), tolerance = 1e-10)
  expect_equal(pooled$variance[["m"]], var(x), tolerance = 1e-8)
})
```

**Step 2 — Run, expect FAIL** (first two fail on current code).

**Step 3 — Fix.** (a) In `welford_mean_var`, after the first sample stores `names0`,
reorder every later sample: `x <- x[names0]` with an `if (any(is.na(match(names0, names(x)))))`
guard that `stop()`s on a name mismatch. (b) In `pool_welford_results`, compute each chunk's
`M2 <- ifelse(n > 1, variance * (n - 1), 0)` (treat n=1 as M2=0) instead of `variance*n*(n-1)`.

**Step 4 — PASS (all three). Step 5 — Commit** `fix: welford name-matching + n=1 pooling; add tests`.

---

### Task A6: Stable `theta_to_id` formatting

**Files:**
- Modify: `R/utils.R:61-63`
- Test: `tests/testthat/test-basic.R` (add) or NEW `test-utils.R`

**Step 1 — Failing test:**

```r
test_that("theta_to_id is stable at bounds edges and scale-invariant", {
  a <- theta_to_id(c(0.0, 0.5, 1.0))
  b <- theta_to_id(c(1e-7, 0.5, 1.0))     # tiny coord must not reformat 0.5
  expect_true(grepl("0.500000", a, fixed = TRUE))
  expect_true(grepl("0.500000", b, fixed = TRUE))
  old <- options(scipen = -5); on.exit(options(old))
  expect_identical(theta_to_id(c(0.0, 0.5, 1.0)), a)  # option-independent
})
```

**Step 2 — Run, expect FAIL** (scientific-notation flip).

**Step 3 — Fix.** Replace `format(round(x,6), nsmall=6)` with per-element fixed formatting:
`paste(formatC(round(unit_theta, 6), format = "f", digits = 6), collapse = "_")`.

**Step 4 — PASS. Step 5 — Commit** `fix: stable per-coordinate theta_to_id formatting`.

---

### Task A7: Multi-seed verification survives the documented contract

**Files:**
- Modify: `R/bo_v04_helpers.R:216, 233` (normalize `res`; move access inside tryCatch)
- Test: `tests/testthat/test-warmstart.R` or NEW `test-multiseed.R` (zero coverage today)

**Step 1 — Failing test:** run `run_multi_seed_verification` (or `bo_calibrate(..., multi_seed_verify = TRUE)`)
with a simulator returning a **named numeric vector** (the documented contract) and assert
it returns a verdict rather than erroring:

```r
test_that("multi-seed verification handles vector-returning simulators", {
  sim <- function(theta, fidelity = "high", seed = NULL, ...) {
    v <- c(power = 0.9, type1 = 0.02, EN = 40); attr(v, "variance") <- c(power=1e-4, type1=1e-5, EN=1)
    v
  }
  expect_no_error(
    run_multi_seed_verification(theta = list(x = 0.5), sim_fun = sim,
                                constraints = list(power = c("ge", 0.8)),
                                objective = "EN", n = 3, rng_seed = 1)
  )
})
```

**Step 2 — Run, expect FAIL** (`$ operator is invalid for atomic vectors`).

**Step 3 — Fix.** Normalize the simulator return the same way `invoke_simulator` does:

```r
res <- tryCatch(sim_fun(theta = theta, fidelity = "high", seed = seed_i, ...),
                error = function(e) e)
if (inherits(res, "error")) { ... record NA seed ...; next }
metrics_i <- if (is.list(res)) res$metrics else res   # vector OR list contract
```

Also drop the undocumented `n_rep` argument (or gate it behind a `formals` check — see B1).
Ensure every seed access is inside the tryCatch.

**Step 4 — PASS. Step 5 — Commit** `fix: multi-seed verification accepts vector-returning simulators; wrap in tryCatch`.

---

### Task A8: Incumbent falls back to any-fidelity feasible best

**Files:**
- Modify: `R/bo_calibrate.R:1367-1374` (`best_feasible_objective`, high_fidelity_only branch)
- Test: `tests/testthat/test-BATON-core.R`

**Step 1 — Failing test:** construct a history with feasible points only at low/med and one
infeasible high-fidelity row; assert `best_feasible_objective(..., high_fidelity_only=TRUE)`
returns the low/med feasible best, not `Inf`.

**Step 2 — Run, expect FAIL** (returns Inf).

**Step 3 — Fix.** When the feasible-and-high intersection is empty, fall back to the best
feasible at any fidelity (emit a one-line note under `progress`), instead of returning Inf.

**Step 4 — PASS. Step 5 — Commit** `fix: incumbent falls back to any-fidelity feasible best`.

---

### Task A9: Early-stopping robustness (transition + relative threshold)

**Files:**
- Modify: `R/bo_calibrate.R:866-931`
- Test: `tests/testthat/test-phase3-performance.R`

**Step 1 — Failing tests:** (a) a run that becomes feasible late does not trip early-stop
on the infeasible→feasible transition; (b) the acquisition-flatline stop uses a threshold
relative to |incumbent|, so a [0,1]-scale objective is not stopped at iteration 6.

**Step 2 — Run, expect FAIL.**

**Step 3 — Fix.** Track feasible and infeasible phases separately (only compare
improvements within the feasible phase; reset the counter at first feasibility). Replace
`selected_acq_max < 0.1` with `selected_acq_max < rel_tol * max(abs(incumbent), 1e-8)`
(rel_tol default 1e-3, configurable via `early_stop$acq_rel`). Append the incumbent AFTER
this iteration's evaluations are folded in (fixes the one-iteration lag).

**Step 4 — PASS. Step 5 — Commit** `fix: early-stopping transition handling and scale-relative acq threshold`.

---

### Task A10: Per-philosophy error isolation

**Files:**
- Modify: `R/bo_calibrate_philosophies.R:137-146, 231-241`
- Test: `tests/testthat/` (NEW small test or extend existing philosophies test)

**Step 1 — Failing test:** run two philosophies where the first errors (inject a
constraint-name typo) and assert the second still completes and the result carries a
per-philosophy status rather than the whole call aborting.

**Step 2 — Run, expect FAIL.**

**Step 3 — Fix.** Wrap each `bo_calibrate()` in `tryCatch`; store
`list(status = "errored", error = conditionMessage(e))` for failures and continue. Fix the
donor-name/spec misalignment (A12 territory) if touched.

**Step 4 — PASS. Step 5 — Commit** `fix: isolate per-philosophy failures; keep completed fits`.

---

### Task A11: Donor warmstart pairs best_theta with its own metrics

**Files:**
- Modify: `R/bo_calibrate_philosophies.R:177-186`
- Test: extend the philosophies test

**Step 1 — Failing test:** donor fit whose `best_theta` differs from the last high-fidelity
feasible row; assert the built spec's metrics come from best_theta's own row (match by
`theta_id`).

**Step 3 — Fix.** Look up the row via `theta_to_id(scale_to_unit(best_theta, bounds))`
against `history$theta_id`; use that row's metrics/feasible. Fall back to current behavior
only if no match.

**Step 5 — Commit** `fix: donor warmstart uses best_theta's own evaluation`.

---

### Task A12: Cheap correctness cleanups (batch commit)

**Files/edits:**
- `R/benchmark.R:254,289,322` — `unname()` the fidelity source before `c(high = ...)`.
  Test: benchmark arm smoke test with the documented named-vector input, expect no
  "subscript out of bounds".
- `R/bo_parameter_importance.R:272` — remove or rewire `compute_global_sensitivity` to
  `sa_sobol`; if removed, drop its `@export` and NAMESPACE entry. Test: `expect_error`
  gone / function returns a data frame.
- `R/case_study.R:55` — use `fit$bounds` instead of `infer_bounds_from_history`; add a
  zero-width guard in `scale_to_unit` callers. Test: constant-parameter history no longer
  yields NaN.
- `R/surrogates.R:346` — impute NA noise with `max(noise_vec, na.rm=TRUE)` (conservative).
- `R/bo_calibrate.R:304`, `R/surrogates.R:95` — gate one-time messages on `progress`;
  leave the option flags but note they are process-global (service uses subprocess isolation).
- Add input asserts: `parse_constraints` errors on non-numeric/non-finite threshold and on
  unnamed constraint lists (`R/constraints.R:14-19`); `bo_calibrate` errors when
  `objective %in% names(bounds)` (name collision, item 1 table).

One test file `tests/testthat/test-bug-fixes.R` gains a case per edit. Commit as
`fix: batch of low-risk correctness cleanups (benchmark, sensitivity, case_study, validation)`.

---

# PHASE B — Public-contract changes (REQUIRES PI DECISION before coding)

These change behavior users and the paper's simulator depend on. Do not start until the
two questions are answered.

### Task B1: Pass `n_rep` to the simulator (review item 1.4) — DECISION: add optional arg

**PI decision 2026-07-22: extend the contract with an optional `n_rep` argument.**
In `invoke_simulator` detect `"n_rep" %in% names(formals(sim_fun)) || "..." %in% names(formals(sim_fun))`
and pass `n_rep` only when accepted (backward-compatible; simulators without it are
unaffected). Update the roxygen contract to document the optional argument, update the
`inst/examples/*` simulators to accept and honor it, and align default fidelity levels with
the examples. Test: a simulator that records its received `n_rep` sees the escalated value
under `hybrid_staged`; a legacy simulator without the argument still runs unchanged.

### Task B2: `warmstart_from` augments the LHS design (item 1.6) — DECISION: augment

**PI decision 2026-07-22: donor seeds augment a full LHS initial design.**
When `nrow(initial_history) < n_init`, run LHS for the remaining `n_init - nrow(initial_history)`
points, dedup against donor `theta_id`s, and combine (donor rows carry iter 0). Test: a
2-row donor history with `n_init = 60` yields a 60-row initial design that includes the 2
donor rows and 58 fresh LHS points with no duplicate `theta_id`s. Update the docs to state
the augment behavior explicitly.

---

# PHASE C — Performance (behavior-preserving, regression-gated, runs first)

### Task C0: Golden-master calibration snapshot (regression harness)

**Files:** NEW `tests/testthat/test-regression-golden.R` with inlined reference values
(no `_snaps/` directory; the references live in the test for transparency and version control).

**Step 1 — Build it.** A deterministic toy calibration (fixed seed, small budget, the
existing `toy_sim_fun`) that records `history$objective`, `history$feasible`, `best_theta`,
and `best_objective` via `expect_equal` against inlined reference values.

**Execution-order note (2026-07-22):** we run behavior-preserving performance work FIRST,
so the reference is captured on the **current baseline** (`main @ dc25ca3`), not
post-Phase-A. Every Phase C task must keep it byte-identical. The Phase A tasks that
intentionally change the search trajectory (A1/A2 warm-start, A3 penalization, A4 expected
violation) re-baseline this snapshot in the SAME commit that makes the change, with the
delta explained in the commit message. Fixes that do not affect the toy run (e.g. welford,
theta_id edge formatting, multi-seed) leave it untouched.

**Step 2 — Confirm determinism:** run twice, assert identical. Commit the test.

Every Phase C task below re-runs this and must keep it identical (except tasks that
intentionally change results — none here; performance work is behavior-preserving).

**Commit** `test: golden-master calibration snapshot to guard performance refactors`.

---

### Task C1: Parallelize the per-metric GP fits — DONE (commit f74cc00)

**Implemented 2026-07-22.** `run_seeded()` wraps each km/hetGP fit under a deterministic
seed threaded from `bo_calibrate` (`fit_seed = rng_seed + 10000*iter`, per-metric offset).
`fit_surrogates` uses an opt-in `mclapply` mapper gated on `options(BATON.cores)`. Verified:
serial == cores=4 (new core-independence regression test); golden master re-baselined; full
suite 62 pass. Original blocked-analysis retained below for context.

---
_Original analysis (now resolved):_

**Files:** `R/surrogates.R:58-118` (`fit_surrogates` metric loop).

**Finding (2026-07-22):** naive `mclapply` here is NOT result-neutral. km/hetGP use random
hyperparameter starts; forked workers get independent RNG streams, so `options(BATON.cores=4)`
produced a different (though per-config deterministic) calibration than serial
(best_obj 0.02051 vs 0.02026 on the golden run). Core-count-dependent results are
unacceptable for a reproducible cloud backbone.

**Resolution:** parallelize only after making each fit reproducible by construction.
Threading (the missing API path codex flagged):
1. `fit_surrogates(...)` gains a `fit_seed = NULL` argument.
2. `bo_calibrate` passes `fit_seed = rng_seed + 10000L * iter_counter` at the in-loop call
   sites and `rng_seed + 99999L` for the single post-loop final fit (the `10000*iter` offset
   keeps per-iteration fit seeds disjoint from the per-eval simulator seeds
   `rng_seed + eval_counter`).
3. Inside the metric loop, derive a per-fit seed `metric_seed <- fit_seed + match(metric, metrics_needed)`
   and pass it to `fit_dicekriging_surrogate`/`fit_hetgp_surrogate`, which each add a
   `fit_seed = NULL` parameter and wrap the km/mleHetGP call in
   `run_seeded(fit_seed, ...)`, a base-R save/restore-`.Random.seed` wrapper (no withr dep).
   When `fit_seed` is NULL the helpers behave exactly as today (no seeding).
This makes serial and parallel identical regardless of core count. As implemented, the
`mclapply` mapper is opt-in via `options(BATON.cores)` and the golden master was re-baselined
in the same commit. Test (shipped): identical results for `BATON.cores` in {1, 4}.

---

### Task C2: Parallelize the simulator batch and initial design

**Files:** `R/bo_calibrate.R:535-594` (init loop), `:730-860` (batch loop). Split
`record_evaluation` into `invoke_evaluation` (pure, returns a row) and `append_rows`.

Evaluate the `q` batch points and the `n_init` initial points with `mclapply` (seeds are
already `rng_seed + eval_counter`, so results are seed-deterministic and order-independent),
then append all rows in one bind. Test: golden-master identical serial vs `BATON.cores=2`;
assert the produced history rows match the sequential order after sorting by `eval_id`.

**Commit** `perf: evaluate simulator batch and initial design in parallel`.

---

### Task C3: Warm-start actually reduces fit cost (validate A1/A2 payoff)

**Files:** none new; a benchmark test.

Add `tests/testthat/test-warmstart.R` timing assertion (skip_on_cran) that a warm-started
km fit does fewer BFGS iterations / is not slower than cold on a mid-size design. This
confirms A1/A2 delivered the intended speedup rather than just correctness.

**Commit** `test: confirm warm-start reduces surrogate fit cost`.

---

### Task C4: Objective-only incumbent prediction

**Files:** `R/bo_calibrate.R:1378-1383` (`best_feasible_objective`).

Predict only the objective surrogate (`surrogates[objective]`) over feasible points instead
of all `m`. Cache the feasible design matrix and extend it by the new rows each iteration.
Test: golden-master identical; micro-benchmark shows fewer predict calls.

**Commit** `perf: predict only the objective surrogate for the incumbent`.

---

### Task C5: Matrix end-to-end candidate scoring

**Files:** `R/surrogates.R:476-515` (`predict_surrogates`), `R/bo_calibrate.R:1340-1346`
(`lhs_candidate_pool`), call sites.

Let `predict_surrogates` accept a numeric matrix (columns in `param_names` order) directly,
skipping the `asplit`→named-vector→rebuild round trip and the per-candidate tryCatch. Keep
the list-input path for back-compat. Convert the invalid-row silent drop into an error
(item 40-score). Test: golden-master identical; assert matrix and list inputs give equal
predictions.

**Commit** `perf: vectorized matrix path through predict_surrogates`.

---

### Task C6: Incremental / matrix history aggregation

**Files:** `R/surrogates.R:49-50, 278-305`.

Replace the per-theta one-row-tibble + `bind_rows` aggregation with a single matrix pass:
build `values` (n×m) and `variance` (n×m) matrices from the list-columns once, then compute
grouped means. **Must preserve the current `mean(..., na.rm = TRUE)` semantics** (codex
review flagged that plain `rowsum`/`tabulate` propagates NA and would poison grouped
surrogate inputs). Use an NA-aware group mean: e.g. replace NA with 0 in a copy for the
numerator, count non-NA per group with `rowsum(!is.na(values), theta_ids)` for the
denominator, and divide (guarding zero denominators → NA, matching `mean(NA, na.rm=TRUE)` =
NaN for an all-NA group). Optionally cache the aggregated design across iterations and
update only new `theta_id`s. Test: golden-master identical; unit test that the matrix
aggregation equals the old group-by on a fixture **that includes NA metric and NA noise
entries**.

**Commit** `perf: matrix-based history aggregation (was ~1200 tibble allocs/iter)`.

---

### Task C7: Row-list history append + hoisted loop invariants

**Files:** `R/bo_calibrate.R:1206-1222` (append), `:780` (budget sum), `:784-793`
(distance_to_best).

Accumulate new rows in a preallocated list and `bind_rows` once per iteration instead of
per evaluation (kills the O(n²) copy). Read `cumulative_budget_used` instead of
`sum(history$n_rep)` in the inner loop. Compute `distance_to_best` by row index once per
iteration (also fixes the float-equality dead code, item 65-score). Test: golden-master
identical.

**Commit** `perf: batch history append and hoist per-iteration invariants`.

---

### Task C8: Batched sensitivity/diagnostics predicts (optional, for case studies)

**Files:** `R/sensitivity.R:265-279` (`cov_effects`/`sa_gradients`), `:257-263` (`sa_sobol`).

Batch all `2 * n_mc * d` gradient query points into one predict per metric; drop the
per-row list conversion in `sa_sobol`. Test: results equal the current per-call version
within tolerance; micro-benchmark shows the speedup.

**Commit** `perf: batch predicts in sensitivity/diagnostics`.

---

# PHASE D — Service readiness (after C; separate plan/PR)

Deferred to a follow-up plan: partial-fit return with `status`, `checkpoint_fun`,
`callback`+cancellation, `max_walltime_s`, RNG save/restore (or callr-per-job at the service
layer), dependency demotion of ggplot2/purrr/tidyr/rlang/randtoolbox to Suggests, and a
`slim = TRUE` return. These are cloud-hosting features, not correctness or single-run speed,
and are best done once the calibration core is settled.

---

## Verification gate (before merging the branch)

1. `Rscript -e 'devtools::test()'` — all green, including new welford/theta_id/multiseed
   tests and the golden-master snapshot.
2. `Rscript -e 'devtools::check(document = TRUE)'` — no new WARNINGs/ERRORs.
3. Re-run the golden-master after the full Phase C stack: byte-identical to the
   post-Phase-A reference (proves the speedups did not change answers).
4. Update `NEWS.md`, bump DESCRIPTION Version, and note the corrected warm-start /
   penalization behavior (the old CLAUDE.md/NEWS efficiency claims were describing
   behavior that now actually exists).
