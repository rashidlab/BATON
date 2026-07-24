# Phase D: Service Readiness Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make `bo_calibrate()` safe to host as a long-running cloud service: no
work lost on error, cooperative cancellation, walltime caps, periodic
checkpoints, RNG hygiene, a lightweight return option, and a slimmer hard
dependency set.

**Architecture:** One structural refactor unlocks everything: the `BATON_fit`
construction (currently inlined at the end of `bo_calibrate()`,
R/bo_calibrate.R:1183) is extracted into a `build_baton_fit()` helper that can
be called (a) at normal completion, (b) on error/cancel/walltime with whatever
history exists, and (c) by the checkpoint hook. A `status` field records why
the run ended. All new behavior is opt-in with back-compatible defaults; the
golden-master snapshot must stay byte-identical throughout (every task here is
either additive or off-by-default). Dependency demotion is last because it
touches only analysis modules.

**Tech Stack:** R 4.6.0, testthat 3e. Test command:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet = TRUE)); testthat::test_file("tests/testthat/<file>.R")'`
Full suite: same with `testthat::test_dir("tests/testthat")`; NOT_CRAN runs
prefix `NOT_CRAN=true`. Check:
`R CMD build . && NOT_CRAN=true _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual <tarball>`
(devtools is NOT installed on this machine; use pkgload/roxygen2 directly.)

**Branch:** `feat/phase-d-service` off `main` (v0.6.0). Commit after every task.

**Reference:** `CODE_REVIEW_2026-07-22.md` Section on cloud blockers;
`docs/plans/2026-07-22-review-fixes-and-performance.md` PHASE D stub.

---

## Cross-cutting rules

- **TDD:** every task writes its failing test first (except pure-extraction
  Task D1, which is behavior-preserving and gated by the existing suite).
- **Golden master:** `tests/testthat/test-regression-golden.R` must pass
  UNCHANGED after every task. No re-baselining in Phase D — nothing here may
  alter numerical results.
- **Back-compat:** every new `bo_calibrate()` argument defaults to current
  behavior (`on_error = "stop"`, `max_walltime_s = NULL`, `callback = NULL`,
  `checkpoint_fun = NULL`, `slim = FALSE`).
- **Docs:** roxygen for each new argument in the same commit; regenerate with
  `Rscript -e 'suppressMessages(roxygen2::roxygenise())'`.

---

### Task D1: Extract `build_baton_fit()` (behavior-preserving)

**Files:**
- Modify: `R/bo_calibrate.R:1183-1230` (the `structure(list(...), class = "BATON_fit")` block)
- Test: existing suite is the gate (pure extraction, no new test)

**Step 1:** Read the current return block (starts `structure(` at
R/bo_calibrate.R:1183; runs through the `class = "BATON_fit"` line). Note every
local it references (`history`, `best_theta`, `best_objective`,
`final_surrogates`, `acquisition`, `q`, `budget`, `rng_seed`,
`fidelity_levels`, `fidelity_method`, `fidelity_costs`, `candidate_pool`,
`objective`, `early_stop_config`, `multi_seed_*`, `diagnostics`, `verdict`,
`bounds`, `constraints`, ...).

**Step 2:** Create internal helper directly above `initialise_history()`:

```r
#' Assemble a BATON_fit object (Task D1)
#'
#' Single construction point so normal completion, partial returns
#' (error/cancel/walltime), and checkpoints produce structurally identical
#' objects. `status` records why the run ended.
#' @keywords internal
build_baton_fit <- function(history, best_theta, best_objective, surrogates,
                            policies, diagnostics, multi_seed_summary,
                            multi_seed_runs, verdict, bounds, constraints,
                            constraint_tbl, status, error_message = NULL) {
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
    class = c("BATON_fit", "list")   # keep BOTH classes: code and tests
                                     # dispatch on the existing class vector
  )
}
```

CAUTION: the field list, field ORDER, and the two-element class vector above
must match the current inline block at R/bo_calibrate.R:1183 exactly (as of
v0.6.0 it ends `bounds, constraints, constraint_tbl` with
`class = c("BATON_fit", "list")`); re-read that block before writing the
helper in case it has drifted, and append only `status`/`error_message` after
the existing fields. Build the `policies` list at the call site and pass it
in.

**Step 3:** Replace the inline `structure(...)` with a `build_baton_fit(...)`
call passing `status = "completed"` (placeholder until D2 refines it).

**Step 4:** Run the full suite:
`Rscript -e 'suppressMessages(pkgload::load_all(".", quiet = TRUE)); testthat::test_dir("tests/testthat")'`
Expected: same pass count as before (golden master included).

**Step 5:** Commit `refactor: extract build_baton_fit construction helper (D1)`.

---

### Task D2: `status` field records why the run ended

**Files:**
- Modify: `R/bo_calibrate.R` (main loop exit points; `build_baton_fit` call)
- Test: NEW `tests/testthat/test-status.R`

**Step 1 — Failing tests.** A `BATON_fit` must carry `status` with one of:
`"budget_exhausted"` (while-loop ran out), `"early_stopped"` (improvement
plateau), `"acq_flatline"` (acquisition threshold stop):

```r
test_that("status is budget_exhausted when the budget is consumed", {
  fit <- <small run, early_stop disabled, budget reached>
  expect_equal(fit$status, "budget_exhausted")
})
test_that("status is early_stopped when the improvement plateau trips", {
  fit <- <run with early_stop = list(enabled TRUE, patience 1, threshold 10, consecutive 1)
          on a flat deterministic objective so it stops immediately>
  expect_equal(fit$status, "early_stopped")
})
```

(For `acq_flatline`, trigger via `early_stop$acq_rel` with a huge value; if
too fiddly to trigger reliably, assert only that the recorded status is one of
the documented values and add a unit test on the small helper that maps the
loop-exit reason to a status string.)

**Step 2:** Run, expect FAIL (`fit$status` is `"completed"`/NULL).

**Step 3 — Implement.** Introduce `run_status <- "budget_exhausted"` before the
while loop; each `break` site sets it first (`"early_stopped"`,
`"acq_flatline"`). Pass `run_status` to `build_baton_fit`. Document the field
in the roxygen `@return`.

**Step 4:** Run tests + golden. **Step 5:** Commit
`feat: BATON_fit$status records termination reason (D2)`.

---

### Task D3: RNG hygiene — do not clobber the caller's stream

**Files:**
- Modify: `R/bo_calibrate.R` (function ENTRY — see the seed = NULL trap below)
- Test: NEW `tests/testthat/test-rng-hygiene.R`

**Trap:** when `seed = NULL`, the code draws
`rng_seed <- sample.int(1e6, 1)` (R/bo_calibrate.R:352) BEFORE `set.seed` —
that draw itself advances the caller's stream. The restore hook must be
installed as the FIRST statement of `bo_calibrate`, before any RNG use, not
merely before `set.seed(rng_seed)`.

**Step 1 — Failing tests** (cover BOTH seed paths):

```r
test_that("bo_calibrate leaves the caller's RNG stream untouched (explicit seed)", {
  set.seed(123); expected <- runif(3)
  set.seed(123)
  invisible(suppressMessages(bo_calibrate(<tiny deterministic run>, seed = 42)))
  expect_equal(runif(3), expected)
})
test_that("caller's stream is untouched even with seed = NULL", {
  set.seed(123); expected <- runif(3)
  set.seed(123)
  invisible(suppressMessages(bo_calibrate(<tiny deterministic run>, seed = NULL)))
  expect_equal(runif(3), expected)
})
```

**Step 2:** Run, expect FAIL (stream advanced by internal draws; the
seed = NULL case fails even with a hook placed after the sample.int draw).

**Step 3 — Fix.** As the first statement of `bo_calibrate`:

```r
if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
  old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
} else {
  on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)), add = TRUE)
}
```

(Same pattern as `run_seeded()`, R/surrogates.R:304. Internal behavior is
unchanged — `set.seed(rng_seed)` still runs — only the caller's stream is
restored on exit, so the golden master cannot move.)

**Step 4:** PASS + golden. **Step 5:** Commit
`fix: bo_calibrate restores the caller's RNG stream (D3)`.

---

### Task D4: `max_walltime_s` graceful stop

**Files:**
- Modify: `R/bo_calibrate.R` (signature + top-of-iteration check)
- Test: NEW `tests/testthat/test-service-controls.R`

**Step 1 — Failing test.** Use a simulator with a built-in `Sys.sleep(0.05)`
and `max_walltime_s = 0.1`, generous budget; assert the run returns (not
errors) with `status = "walltime"` and a non-empty partial history:

```r
test_that("max_walltime_s stops gracefully with a partial fit", {
  slow_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    Sys.sleep(0.05)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01); v
  }
  fit <- suppressMessages(bo_calibrate(
    slow_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 50, seed = 1, max_walltime_s = 0.4,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "walltime")
  expect_gte(nrow(fit$history), 3)
  expect_lt(nrow(fit$history), 50)
})
```

**Step 2:** FAIL (unused argument).

**Step 3 — Implement.** `max_walltime_s = NULL` in the signature; record
`start_time <- Sys.time()` at entry; at the TOP of each while iteration:

```r
if (!is.null(max_walltime_s) &&
    as.numeric(difftime(Sys.time(), start_time, units = "secs")) > max_walltime_s) {
  run_status <- "walltime"
  if (progress) message(sprintf("Walltime cap (%.0fs) reached at iteration %d.",
                                max_walltime_s, iter_counter))
  break
}
```

The check runs between iterations only (never mid-batch), so a partial fit is
always internally consistent. Best-theta/final-fit code after the loop already
handles whatever history exists. Document that the cap is advisory
(iteration-granular), not preemptive.

**Step 4:** PASS + golden + full suite. **Step 5:** Commit
`feat: max_walltime_s graceful stop (D4)`.

---

### Task D5: `callback` with cooperative cancellation

**Files:**
- Modify: `R/bo_calibrate.R` (signature + end-of-iteration hook)
- Test: `tests/testthat/test-service-controls.R`

**Step 1 — Failing tests:**

```r
test_that("callback sees per-iteration progress and can cancel", {
  seen <- list()
  cb <- function(info) { seen[[length(seen) + 1]] <<- info; TRUE }
  fit <- <small run with callback = cb>
  expect_gte(length(seen), 1)
  expect_true(all(vapply(seen, function(s)
    all(c("iter", "eval_count", "best_objective", "elapsed_s") %in% names(s)),
    logical(1))))

  cancel_after_2 <- function(info) info$iter < 2
  fit2 <- <same run with callback = cancel_after_2>
  expect_equal(fit2$status, "cancelled")
  expect_lt(nrow(fit2$history), <budget>)
})
```

**Step 2:** FAIL.

**Step 3 — Implement.** `callback = NULL` in signature. At the END of each
iteration (after `append_rows` and budget update):

```r
if (!is.null(callback)) {
  keep_going <- tryCatch(
    isTRUE(callback(list(
      iter = iter_counter,
      eval_count = eval_counter,
      best_objective = if (is.finite(best_feasible_value)) best_feasible_value else NA_real_,
      elapsed_s = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    ))),
    error = function(e) {
      warning(sprintf("callback error (continuing): %s", conditionMessage(e)))
      TRUE
    })
  if (!keep_going) { run_status <- "cancelled"; break }
}
```

Contract: callback returns `TRUE` to continue; anything else cancels after the
current iteration (no mid-batch interruption; no work discarded). A callback
error warns and continues — a monitoring hook must not kill a paid run.

**Step 4:** PASS + golden. **Step 5:** Commit
`feat: per-iteration callback with cooperative cancellation (D5)`.

---

### Task D6: `checkpoint_fun` periodic snapshots

**Files:**
- Modify: `R/bo_calibrate.R` (signature + hook next to the callback)
- Test: `tests/testthat/test-service-controls.R`

**Step 1 — Failing test:**

```r
test_that("checkpoint_fun receives resumable partial fits", {
  checkpoints <- list()
  ckpt <- function(fit) checkpoints[[length(checkpoints) + 1]] <<- fit
  fit <- <run with checkpoint_fun = ckpt, checkpoint_every = 1, budget 8, n_init 4>
  expect_gte(length(checkpoints), 2)
  last <- checkpoints[[length(checkpoints)]]
  expect_s3_class(last, "BATON_fit")
  expect_equal(last$status, "checkpoint")
  # the checkpoint's history must be a usable initial_history for a resume
  resumed <- suppressMessages(bo_calibrate(
    <same problem>, initial_history = last$history,
    n_init = nrow(last$history), budget = nrow(last$history) + 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_gte(nrow(resumed$history), nrow(last$history) + 1)
})
```

**Step 2:** FAIL.

**Step 3 — Implement.** `checkpoint_fun = NULL`, `checkpoint_every = 5L` in
signature. After the callback hook, every `checkpoint_every` iterations, call
`checkpoint_fun(build_baton_fit(..., surrogates = NULL, diagnostics = NULL,
status = "checkpoint"))` — checkpoints are slim by construction (no GP
objects: they are large, environment-laden, and not needed to resume, since a
resume refits from history via `initial_history`). Wrap in `tryCatch` →
warn-and-continue, same rationale as D5. Document the resume recipe
(`initial_history = ckpt$history`) in roxygen; note the B2 augment rule (pass
`n_init <= nrow(history)` to avoid LHS top-up on resume).

**Step 4:** PASS + golden. **Step 5:** Commit
`feat: checkpoint_fun periodic resumable snapshots (D6)`.

---

### Task D7: Partial fit on simulator error (`on_error = "return_partial"`)

**Files:**
- Modify: `R/bo_calibrate.R` (signature; wrap the init/BO evaluation calls)
- Test: `tests/testthat/test-service-controls.R`

**Step 1 — Failing tests:**

```r
test_that("on_error='return_partial' preserves completed work", {
  n_calls <- 0
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    if (n_calls > 6) stop("simulator infrastructure died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01); v
  }
  fit <- suppressMessages(bo_calibrate(
    flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 12, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "infrastructure died")
  expect_gte(nrow(fit$history), 4)   # completed evaluations survive
})
test_that("default on_error='stop' behavior is unchanged", {
  <same flaky sim>
  expect_error(<same call without on_error>, "infrastructure died")
})
test_that("return_partial also survives a surrogate-fit failure", {
  # fit_dicekriging_surrogate stop()s on ill-conditioned data
  # (R/surrogates.R:476) -- an error AFTER evaluations exist must not lose
  # them. Force it with testthat::local_mocked_bindings on fit_surrogates
  # erroring from the second call on (the first call happens in-loop only
  # after the init design is complete).
  calls <- 0
  testthat::local_mocked_bindings(
    fit_surrogates = function(...) {
      calls <<- calls + 1
      if (calls >= 2) stop("Failed to fit surrogate for metric 'EN': fake")
      <real fit_surrogates>(...)
    },
    .package = "BATON"
  )
  fit <- <run with on_error = "return_partial", budget comfortably above n_init>
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "Failed to fit surrogate")
  expect_gte(nrow(fit$history), <n_init>)  # completed evaluations survive
})
```

**Step 2:** FAIL.

**Step 3 — Implement.** `on_error = c("stop", "return_partial")` +
`match.arg`. Coverage must be the WHOLE per-iteration body, not just the
simulator calls: surrogate fitting can `stop()` mid-loop on ill-conditioned
data (R/surrogates.R:476) after plenty of history exists. Wrap (a) the
init-design chunk loop and (b) the entire while-iteration body (fit ->
acquisition -> fidelity -> evaluate -> append) in `tryCatch`; on error with
`on_error == "return_partial"`, set `run_status <- "errored"`,
`run_error <- conditionMessage(e)`, and break (an `interrupted <- TRUE` flag
both loops check). NOTE: a failed parallel batch loses that batch
(evaluate_points re-throws after mclapply) — completed PRIOR rows survive
because appends are per-round; document exactly that. After the loop, skip
multi-seed verification when `run_status == "errored"`, attempt the final
surrogate fit inside its own `tryCatch` (surrogates may be NULL in the
partial fit), and return via `build_baton_fit(..., status = run_status,
error_message = run_error)`.

**Step 4:** PASS + golden + full suite. **Step 5:** Commit
`feat: on_error='return_partial' preserves completed evaluations (D7)`.

---

### Task D8: `slim = TRUE` lightweight return

**Files:**
- Modify: `R/bo_calibrate.R` (signature; final `build_baton_fit` call)
- Test: `tests/testthat/test-service-controls.R`

**Step 1 — Failing test:**

```r
test_that("slim = TRUE drops heavy components but keeps decisions", {
  full <- <small run>
  slim <- <same run with slim = TRUE>
  expect_null(slim$surrogates)
  expect_null(slim$diagnostics)
  expect_equal(slim$best_theta, full$best_theta)
  expect_equal(slim$history$objective, full$history$objective)
  expect_lt(as.numeric(object.size(slim)), as.numeric(object.size(full)) / 2)
})
```

**Step 2:** FAIL.

**Step 3 — Implement.** `slim = FALSE` in signature. When TRUE, the final
`build_baton_fit` gets `surrogates = NULL, diagnostics = NULL,
multi_seed_runs = NULL` (keep `multi_seed_summary` and `verdict` — they are
the adoption gate). Roxygen: state that `sa_sobol`/`sensitivity_diagnostics`/
`extract_case_study` need a non-slim fit, and that km/hetGP objects carry
closures/environments that inflate `saveRDS` payloads — the reason a service
wants `slim`.

**Step 4:** PASS + golden. **Step 5:** Commit
`feat: slim = TRUE lightweight BATON_fit (D8)`.

---

### Task D9: Demote analysis-only dependencies to Suggests

Current usage (verified 2026-07-23): `purrr` — 21 calls in ablation.R (2),
benchmark.R (10), case_study.R (6), reliability.R (3); `tidyr` — 3 calls
(benchmark.R:106-107, case_study.R:131); `ggplot2` — plot functions in 6
analysis files; `randtoolbox` — init_stopping.R:271 (already
requireNamespace-guarded); `rlang` — only `@importFrom rlang .data` (`.data`
is re-exported by dplyr). The calibration core (bo_calibrate, surrogates,
acquisition, constraints, utils, welford) is already purrr/tidyr/ggplot2-free.

**Files:**
- Modify: `DESCRIPTION`, `R/ablation.R`, `R/benchmark.R`, `R/case_study.R`,
  `R/reliability.R`, `R/sensitivity.R`, `R/init_stopping.R`
- Test: existing analysis tests are the gate + NEW guard tests

**Step 1 — Failing test** (NEW `tests/testthat/test-suggests-guards.R`):

```r
test_that("plot functions error informatively, not cryptically, without ggplot2", {
  # simulate absence via mocked requireNamespace is brittle; instead assert
  # every exported plot_* function starts with a require_suggests() guard
  for (fn_name in grep("^plot_", getNamespaceExports("BATON"), value = TRUE)) {
    body_txt <- paste(deparse(body(get(fn_name, asNamespace("BATON")))), collapse = " ")
    expect_match(body_txt, "require_suggests", info = fn_name)
  }
})
```

**Step 2:** FAIL (no `require_suggests`).

**Step 3 — Implement**, one sub-commit per bullet:
1. Add internal `require_suggests(pkg)` in R/utils.R:
   `if (!requireNamespace(pkg, quietly = TRUE)) stop(sprintf("Package '%s' is required for this function. Install it with install.packages(\"%s\").", pkg, pkg), call. = FALSE)`
   and call it at the top of every `ggplot2`-using plot function (6 files).
2. Replace the 21 `purrr::` calls with base R (`lapply` + `dplyr::bind_rows`
   for `map_dfr`/`imap_dfr`; `vapply` for `map_dbl`/`map_lgl`; `Map` for
   `map2`). The same mechanical translation was already done in the core
   files — copy that idiom (see R/utils.R comments).
3. Replace the 3 `tidyr::` calls: `unnest_wider(metrics)` becomes the
   vapply-to-columns idiom used in `record_evaluation`'s `metrics_df` block
   (R/bo_calibrate.R:1281); `replace_na` becomes `ifelse(is.na(x), NA_real_, x)`
   or plain indexing.
4. Swap `@importFrom rlang .data` to `@importFrom dplyr .data` in the 5 files
   using it; drop rlang from Imports.
5. DESCRIPTION: move `ggplot2`, `purrr`, `tidyr`, `randtoolbox` from Imports to
   Suggests; drop `rlang`. Regenerate NAMESPACE (roxygenise).

**Step 4:** Full NOT_CRAN suite + golden. Then the real gate:
`R CMD build . && NOT_CRAN=true _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual <tarball>`
Expected: Status OK (imports actually unused would surface here).

**Step 5:** Commit
`refactor: demote ggplot2/purrr/tidyr/randtoolbox to Suggests; drop rlang (D9)`.

---

### Task D10: Verification gate + release prep

1. Full NOT_CRAN suite green; golden master UNCHANGED from v0.6.0 values.
2. `R CMD check` Status OK (as in D9 Step 4).
3. NEWS.md: new `## BATON 0.7.0` section — service controls
   (`status`/`on_error`/`max_walltime_s`/`callback`/`checkpoint_fun`/`slim`),
   RNG hygiene, dependency slimming (installation gets lighter; plotting and
   Sobol init need Suggests packages).
4. DESCRIPTION Version 0.7.0; README header/highlight/citation to 0.7.0
   (three spots: lines ~5, ~7, ~718 — the v0.6.0 release missed these once).
5. Commit `release: v0.7.0 - service readiness (Phase D)`.

**Explicitly out of scope** (per the original plan): callr-per-job process
isolation — that belongs to the service layer, not the package; the RNG
save/restore (D3) plus `options(BATON.cores)` isolation make in-process
hosting safe enough for the paper's cloud demo.
