# Tests for service-readiness run controls (Task D4): max_walltime_s;
# (Task D5): per-iteration callback with cooperative cancellation;
# (Task D6): checkpoint_fun periodic resumable snapshots.

test_that("effective_cores centralizes the parallel-forking guard", {
  old <- options(BATON.cores = NULL)
  on.exit(options(old), add = TRUE)
  # Option unset: serial.
  expect_identical(BATON:::effective_cores(), 1L)
  # Explicit serial opt-out.
  options(BATON.cores = 1)
  expect_identical(BATON:::effective_cores(), 1L)
  # Opted in: the worker count only when the platform can actually fork
  # (Unix + parallel); otherwise 1L. The Windows case exercises the same
  # else-branch as parallel-unavailable, so no platform shadowing is needed.
  options(BATON.cores = 4)
  if (.Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    expect_identical(BATON:::effective_cores(), 4L)
  } else {
    expect_identical(BATON:::effective_cores(), 1L)
  }
})

test_that("max_walltime_s stops gracefully with a partial fit", {
  skip_if_not_installed("DiceKriging")
  slow_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    Sys.sleep(0.05)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    slow_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 50, seed = 1, max_walltime_s = 0.4,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "walltime")
  expect_gte(nrow(fit$history), 3)
  expect_lt(nrow(fit$history), 50)
})

test_that("max_walltime_s bounds the initialization phase", {
  skip_if_not_installed("DiceKriging")
  slow_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    Sys.sleep(0.05)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    slow_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 40, budget = 50, seed = 1, max_walltime_s = 0.3,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$status, "walltime")
  expect_lt(nrow(fit$history), 40)
})

test_that("multi-seed verification is skipped when the walltime cap trips", {
  skip_if_not_installed("DiceKriging")
  calls <- new.env()
  calls$n <- 0L
  counting_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    calls$n <- calls$n + 1L
    Sys.sleep(0.05)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    counting_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 50, seed = 1, max_walltime_s = 0.4,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_reps = 10000))
  expect_equal(fit$status, "walltime")
  # No simulator calls beyond the recorded history rows: Stage 4 verification
  # (which would add multi_seed_n extra high-fidelity calls) was skipped.
  expect_equal(calls$n, nrow(fit$history))
  expect_null(fit$multi_seed_summary)
  expect_true(is.na(fit$verdict))
})

test_that("verification is skipped when the cap is exceeded on a budget-filling final chunk", {
  skip_if_not_installed("DiceKriging")
  # Deterministic construction of the "condition-exit" hole: the donor top-up
  # fills the budget exactly (n_init = budget = 4, 1 donor + 3 top-up sim
  # calls), and only the FINAL top-up call crosses the cap. The BO while
  # condition is then false before its walltime check can run, so run_status
  # stays "budget_exhausted" (status records why the run ended); Stage 4 must
  # still be suppressed because elapsed time exceeds the cap.
  calls <- new.env()
  calls$n <- 0L
  final_sleep_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    calls$n <- calls$n + 1L
    if (calls$n == 3L) Sys.sleep(0.4)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    final_sleep_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    warmstart_from = list(list(theta = list(x = 0.5), power = 0.95, EN = 45)),
    n_init = 4, budget = 4, seed = 1, max_walltime_s = 0.25,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_reps = 10000))
  expect_equal(fit$status, "budget_exhausted")
  # 4 history rows = 1 donor (no sim call) + 3 top-up calls; verification
  # would have added multi_seed_n further calls.
  expect_equal(calls$n, nrow(fit$history) - 1L)
  expect_null(fit$multi_seed_summary)
  expect_true(is.na(fit$verdict))
})

test_that("cap crossed on the final budget-filling LHS init chunk keeps status budget_exhausted", {
  skip_if_not_installed("DiceKriging")
  # Standard LHS-path twin of the top-up test above: n_init = budget = 4, so
  # completing the init design ends the run for budget reasons. Only the
  # FINAL init call crosses the cap; the run must report "budget_exhausted"
  # (status = why the run ended, consistent across both init paths and the
  # docs) while Stage 4 verification is still suppressed.
  calls <- new.env()
  calls$n <- 0L
  final_sleep_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    calls$n <- calls$n + 1L
    if (calls$n == 4L) Sys.sleep(0.4)
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    final_sleep_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 4, seed = 1, max_walltime_s = 0.25,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_reps = 10000))
  expect_equal(fit$status, "budget_exhausted")
  expect_equal(calls$n, nrow(fit$history))
  expect_null(fit$multi_seed_summary)
  expect_true(is.na(fit$verdict))
})

test_that("max_walltime_s = NULL (default) imposes no cap", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 5, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_false(identical(fit$status, "walltime"))
})

test_that("max_walltime_s is validated", {
  dummy_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    c(power = 0.95, EN = 50)
  }
  args <- list(
    dummy_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 5, seed = 1, progress = FALSE
  )
  for (bad in list(-1, 0, Inf, NaN, NA_real_, "10", c(1, 2))) {
    expect_error(
      do.call(bo_calibrate, c(args, list(max_walltime_s = bad))),
      "max_walltime_s",
      info = paste("value:", paste(bad, collapse = ","))
    )
  }
})

test_that("callback sees per-iteration progress and can cancel", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  seen <- list()
  cb <- function(info) { seen[[length(seen) + 1]] <<- info; TRUE }
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 9, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    callback = cb))
  expect_gte(length(seen), 1)
  expect_true(all(vapply(seen, function(s)
    all(c("iter", "eval_count", "best_objective", "best_observed",
          "elapsed_s") %in% names(s)),
    logical(1))))
  # Everything is feasible with this simulator, so the observed best
  # (including the just-completed batch) is finite from the first invocation.
  expect_true(all(vapply(seen, function(s) is.finite(s$best_observed),
                         logical(1))))
  # Callback fires once per completed BO iteration: iter strictly increasing,
  # eval_count non-decreasing across invocations.
  iters <- vapply(seen, function(s) as.numeric(s$iter), numeric(1))
  evals <- vapply(seen, function(s) as.numeric(s$eval_count), numeric(1))
  expect_true(all(diff(iters) > 0))
  expect_true(all(diff(evals) >= 0))

  cancel_after_2 <- function(info) info$iter < 2
  fit2 <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 9, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    callback = cancel_after_2))
  expect_equal(fit2$status, "cancelled")
  # Cancellation is cooperative and end-of-iteration: the completed
  # iteration-2 batch is kept (3 init + 2 batches of q = 2 = 7 rows), and
  # nothing beyond iteration 2 ran.
  expect_equal(max(fit2$history$iter), 2L)
  expect_equal(nrow(fit2$history), 7L)
})

test_that("callback best_objective and best_observed are NA before feasibility", {
  skip_if_not_installed("DiceKriging")
  # power >= 2 is unsatisfiable (power is a probability), so no evaluation is
  # ever feasible and both incumbent fields must stay NA_real_.
  infeasible_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  infos <- list()
  fit <- suppressWarnings(suppressMessages(bo_calibrate(
    infeasible_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 2)),
    n_init = 3, budget = 7, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    callback = function(info) { infos[[length(infos) + 1]] <<- info; TRUE })))
  expect_gte(length(infos), 1)
  expect_true(all(vapply(infos, function(s) is.na(s$best_objective),
                         logical(1))))
  expect_true(all(vapply(infos, function(s) is.na(s$best_observed),
                         logical(1))))
})

test_that("a callback error warns and the run continues", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # n_init = 3, budget = 5, q = 2: exactly one BO iteration, so exactly one
  # callback invocation (and one warning). A monitoring-hook error must not
  # kill the run: full budget consumed, natural status.
  expect_warning(
    fit <- suppressMessages(bo_calibrate(
      fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 3, budget = 5, q = 2, seed = 1,
      progress = FALSE, early_stop = list(enabled = FALSE),
      callback = function(info) stop("monitoring backend down"))),
    "callback error")
  expect_equal(fit$status, "budget_exhausted")
  expect_equal(nrow(fit$history), 5)
})

test_that("a cancelled run skips Stage 4 multi-seed verification", {
  skip_if_not_installed("DiceKriging")
  calls <- new.env()
  calls$n <- 0L
  counting_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    calls$n <- calls$n + 1L
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    counting_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 9, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_reps = 10000,
    callback = function(info) FALSE))
  expect_equal(fit$status, "cancelled")
  # No simulator calls beyond the recorded history rows: Stage 4 verification
  # (which would add multi_seed_n extra high-fidelity calls) was skipped.
  expect_equal(calls$n, nrow(fit$history))
  expect_null(fit$multi_seed_summary)
  expect_true(is.na(fit$verdict))
})

test_that("callback must be a function or NULL", {
  dummy_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    c(power = 0.95, EN = 50)
  }
  expect_error(
    bo_calibrate(
      dummy_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 3, budget = 5, seed = 1, progress = FALSE,
      callback = "not a function"),
    "callback")
})

test_that("checkpoint_fun receives resumable partial fits", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  problem <- list(list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)))
  checkpoints <- list()
  ckpt <- function(fit) checkpoints[[length(checkpoints) + 1]] <<- fit
  fit <- suppressMessages(bo_calibrate(
    fast_sim, problem[[1]], problem[[2]], problem[[3]],
    n_init = 4, budget = 8, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    checkpoint_fun = ckpt, checkpoint_every = 1))
  # n_init = 4, budget = 8, q = 2: two BO iterations, one checkpoint each.
  expect_gte(length(checkpoints), 2)
  last <- checkpoints[[length(checkpoints)]]
  expect_s3_class(last, "BATON_fit")
  expect_equal(last$status, "checkpoint")
  # Checkpoints are slim by construction: no GP objects, no diagnostics.
  expect_null(last$surrogates)
  expect_null(last$diagnostics)
  # Incumbent info is the observed feasible best by row.
  expect_true(is.finite(last$best_objective))
  expect_true(is.list(last$best_theta))
  # The checkpoint's history must be a usable initial_history for a resume:
  # n_init = nrow(history) avoids the B2 LHS top-up, so the resume refits
  # from the recorded history and continues the run. Counting the resumed
  # run's simulator calls proves the donor rows are REUSED, not
  # discarded-and-resimulated.
  resume_calls <- 0L
  counting_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    resume_calls <<- resume_calls + 1L
    fast_sim(theta, fidelity = fidelity, seed = seed, ...)
  }
  resumed <- suppressMessages(bo_calibrate(
    counting_sim, problem[[1]], problem[[2]], problem[[3]],
    initial_history = last$history,
    n_init = nrow(last$history), budget = nrow(last$history) + 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_gte(nrow(resumed$history), nrow(last$history) + 1)
  # Only the post-resume evaluations hit the simulator.
  expect_equal(resume_calls, nrow(resumed$history) - nrow(last$history))
  # The donor rows survive verbatim, in order, at the head of the history.
  expect_identical(resumed$history$theta_id[seq_len(nrow(last$history))],
                   last$history$theta_id)
})

test_that("checkpoint_every controls snapshot frequency in BO iterations", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  checkpoints <- list()
  ckpt <- function(fit) checkpoints[[length(checkpoints) + 1]] <<- fit
  # n_init = 3, budget = 11, q = 2: four BO iterations; checkpoint_every = 2
  # fires at iterations 2 and 4 only (iteration-granular, not per-evaluation).
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 11, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    checkpoint_fun = ckpt, checkpoint_every = 2))
  expect_equal(length(checkpoints), 2L)
})

test_that("a checkpoint_fun error warns and the run continues", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  expect_warning(
    fit <- suppressMessages(bo_calibrate(
      fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 3, budget = 5, q = 2, seed = 1,
      progress = FALSE, early_stop = list(enabled = FALSE),
      checkpoint_fun = function(fit) stop("disk full"),
      checkpoint_every = 1)),
    "checkpoint_fun error")
  expect_equal(fit$status, "budget_exhausted")
  expect_equal(nrow(fit$history), 5L)
})

test_that("checkpoint_fun and checkpoint_every are validated", {
  dummy_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    c(power = 0.95, EN = 50)
  }
  args <- list(
    dummy_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 5, seed = 1, progress = FALSE
  )
  expect_error(
    do.call(bo_calibrate, c(args, list(checkpoint_fun = "not a function"))),
    "checkpoint_fun")
  for (bad in list(0, -1, 1.5, NA, c(1, 2), "3")) {
    expect_error(
      do.call(bo_calibrate, c(args, list(
        checkpoint_fun = function(fit) NULL, checkpoint_every = bad))),
      "checkpoint_every",
      info = paste("value:", paste(bad, collapse = ",")))
  }
})

# Task D7: on_error = "return_partial" preserves completed work.

test_that("on_error='return_partial' preserves completed work", {
  skip_if_not_installed("DiceKriging")
  n_calls <- 0
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    if (n_calls > 6) stop("simulator infrastructure died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 12, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "infrastructure died")
  expect_gte(nrow(fit$history), 4)   # completed evaluations survive
  # The partial fit still reports the observed feasible best.
  expect_true(is.finite(fit$best_objective))
  expect_true(is.list(fit$best_theta))
})

test_that("default on_error='stop' propagates the simulator error unchanged", {
  skip_if_not_installed("DiceKriging")
  n_calls <- 0
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    if (n_calls > 6) stop("simulator infrastructure died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # Capture the condition message so we can assert it is BYTE-IDENTICAL to
  # the simulator's own error (no tryCatch wrapper altered it in stop mode).
  err <- tryCatch(
    suppressMessages(bo_calibrate(
      flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 4, budget = 12, seed = 1,
      progress = FALSE, early_stop = list(enabled = FALSE))),
    error = conditionMessage)
  expect_identical(err, "simulator infrastructure died")
})

test_that("return_partial also survives a surrogate-fit failure", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # An error arising AFTER evaluations exist, in the surrogate-fit stage
  # rather than the simulator, must not lose the completed history. The mock
  # errors from the second fit_surrogates call on: call 1 is the iteration-1
  # in-loop fit (real), call 2 the iteration-2 fit (errors), call 3 the
  # in-loop no-warm-start retry (errors -> in-loop stop), call 4 the final
  # refit (errors -> surrogates NULL in the partial fit).
  calls <- 0
  real_fit <- BATON:::fit_surrogates
  testthat::local_mocked_bindings(
    fit_surrogates = function(...) {
      calls <<- calls + 1
      if (calls >= 2) stop("Failed to fit surrogate for metric 'EN': fake",
                           call. = FALSE)
      real_fit(...)
    },
    .package = "BATON"
  )
  fit <- suppressWarnings(suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 10, q = 2, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE))))
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "Failed to fit surrogate")
  expect_gte(nrow(fit$history), 4)   # completed evaluations survive
  # The final refit also fails under the mock, so surrogates are NULL; the
  # ORIGINAL in-loop error message is preserved (not the final-refit one).
  expect_null(fit$surrogates)
  expect_null(fit$diagnostics)
})

test_that("return_partial with an init-phase failure returns an empty-history fit", {
  n_calls <- 0
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    if (n_calls > 2) stop("simulator infrastructure died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 6, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "infrastructure died")
  # Without a walltime cap the whole LHS init design is ONE evaluation round,
  # and a failed round loses that round's rows (documented: appends are
  # per-round), so the partial fit has an empty history but is still a
  # structurally complete BATON_fit.
  expect_equal(nrow(fit$history), 0L)
  expect_null(fit$best_theta)
  expect_true(is.na(fit$best_objective))
  expect_null(fit$surrogates)
})

test_that("an errored run skips Stage 4 multi-seed verification", {
  skip_if_not_installed("DiceKriging")
  calls <- new.env()
  calls$n <- 0L
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    calls$n <- calls$n + 1L
    if (calls$n > 6L) stop("simulator infrastructure died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 12, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 5, multi_seed_reps = 10000))
  expect_equal(fit$status, "errored")
  # 6 successful calls recorded + 1 failing call not recorded; Stage 4
  # (which would add multi_seed_n further calls) was skipped.
  expect_equal(calls$n, nrow(fit$history) + 1L)
  expect_null(fit$multi_seed_summary)
  expect_true(is.na(fit$verdict))
})

test_that("on_error is validated", {
  dummy_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    c(power = 0.95, EN = 50)
  }
  expect_error(
    bo_calibrate(
      dummy_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 3, budget = 5, seed = 1, progress = FALSE,
      on_error = "retry"),
    "'arg' should be one of")
})

test_that("checkpoint_every accepts huge integral values without overflow", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # 2^31 overflows as.integer(); the integrality check must not produce NA.
  n_ckpt <- 0L
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 5, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE),
    checkpoint_fun = function(fit) n_ckpt <<- n_ckpt + 1L,
    checkpoint_every = 2^31))
  expect_equal(fit$status, "budget_exhausted")
  expect_equal(n_ckpt, 0L)
})

test_that("return_partial preserves donor rows when the top-up round fails", {
  skip_if_not_installed("DiceKriging")
  flaky <- function(theta, fidelity = "low", seed = NULL, ...) {
    stop("simulator infrastructure died")
  }
  fit <- suppressMessages(bo_calibrate(
    flaky, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    warmstart_from = list(list(theta = list(x = 0.5), power = 0.95, EN = 45)),
    n_init = 4, budget = 6, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "infrastructure died")
  # The donor row (no simulator call) survives; the failed top-up round is
  # lost as one unit (appends are per-round).
  expect_equal(nrow(fit$history), 1L)
  expect_true(is.list(fit$best_theta))
})

test_that("return_partial reports a final-refit failure after a clean run", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # n_init = budget = 4: the BO loop never runs, so the ONLY fit_surrogates
  # call is the post-loop final refit. Its failure must not lose the
  # completed history: status flips from budget_exhausted to errored with
  # the refit's message, surrogates NULL.
  testthat::local_mocked_bindings(
    fit_surrogates = function(...) {
      stop("Failed to fit surrogate for metric 'EN': fake", call. = FALSE)
    },
    .package = "BATON"
  )
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 4, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE)))
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "Failed to fit surrogate")
  expect_equal(nrow(fit$history), 4L)
  expect_null(fit$surrogates)
  expect_null(fit$diagnostics)
  expect_true(is.finite(fit$best_objective))
})

# Task D8: slim = TRUE lightweight return + D7 review carry-forwards.

test_that("slim = TRUE drops heavy components but keeps decisions", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  args <- list(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 6, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE))
  full <- suppressMessages(do.call(bo_calibrate, args))
  slim <- suppressMessages(do.call(bo_calibrate, c(args, list(slim = TRUE))))
  expect_null(slim$surrogates)
  expect_null(slim$diagnostics)
  # The decisions are untouched: same history, same selected design.
  expect_equal(slim$best_theta, full$best_theta)
  expect_equal(slim$best_objective, full$best_objective)
  expect_equal(slim$history$objective, full$history$objective)
  expect_equal(slim$status, full$status)
  # Persisted fits are unambiguous about why surrogates are NULL.
  expect_true(slim$policies$slim)
  expect_false(full$policies$slim)
  # The point of slim: a materially smaller serialization payload. Compare
  # serialized sizes, not object.size(): object.size() does not follow
  # environments, which is exactly the weight (GP model closures) slim drops.
  expect_lt(length(serialize(slim, NULL)), length(serialize(full, NULL)) / 2)
})

test_that("slim keeps the multi-seed adoption gate but drops per-seed runs", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 6, q = 2, seed = 1, slim = TRUE,
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 3, multi_seed_reps = 200))
  expect_null(fit$surrogates)
  expect_null(fit$diagnostics)
  # multi_seed_runs (per-seed evaluations) are dropped, but the adoption gate
  # (summary + verdict) survives: a service gates on these fields.
  expect_null(fit$multi_seed_runs)
  expect_false(is.null(fit$multi_seed_summary))
  expect_equal(fit$verdict, "MULTI_SEED_PASS")
})

test_that("slim is validated", {
  dummy_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    c(power = 0.95, EN = 50)
  }
  args <- list(
    dummy_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 3, budget = 5, seed = 1, progress = FALSE
  )
  for (bad in list(NA, "TRUE", c(TRUE, TRUE), 1)) {
    expect_error(
      do.call(bo_calibrate, c(args, list(slim = bad))),
      "slim",
      info = paste("value:", paste(bad, collapse = ",")))
  }
})

test_that("Stage 4 per-seed simulator failures never lose the fit", {
  skip_if_not_installed("DiceKriging")
  # Verification-phase calls are identifiable: they come AFTER budget
  # exhaustion. budget = 6, so call 7 is the first Stage 4 per-seed
  # re-evaluation; the whole optimization ran clean and only verification
  # dies. run_multi_seed_verification catches per-seed simulator errors
  # itself (warning + NA row), so a dying simulator yields an intact fit
  # flagged MULTI_SEED_FAIL (exercised here under return_partial; the
  # per-seed catch does not depend on on_error). The run-level guard
  # (next test) only handles failures outside that per-seed catch.
  n_calls <- 0
  dies_on_verify <- function(theta, fidelity = "low", seed = NULL, ...) {
    n_calls <<- n_calls + 1
    if (n_calls > 6) stop("verification backend died")
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  expect_warning(
    fit <- suppressMessages(bo_calibrate(
      dies_on_verify, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 4, budget = 6, q = 2, seed = 1, on_error = "return_partial",
      progress = FALSE, early_stop = list(enabled = FALSE),
      multi_seed_verify = TRUE, multi_seed_n = 3, multi_seed_reps = 200)),
    "Multi-seed eval")
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$status, "budget_exhausted")
  expect_null(fit$error_message)
  # The optimization itself is intact: all budgeted rows survive.
  expect_equal(nrow(fit$history), 6L)
  expect_true(is.list(fit$best_theta))
  expect_true(is.finite(fit$best_objective))
  # All-NA seeds fail the strict gate: the design is flagged, not adopted.
  expect_equal(fit$verdict, "MULTI_SEED_FAIL")
  expect_equal(fit$multi_seed_summary$strict_feas, 0L)
})

test_that("return_partial guards a Stage 4 verification crash", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # A failure that escapes run_multi_seed_verification's per-seed tryCatch
  # (e.g. a malformed simulator return crashing the summary construction)
  # arrives AFTER the most expensive work in the run. Under return_partial
  # it must not lose the fit: best_theta/history are kept, status flips to
  # "errored", and the verification gate fields are cleanly absent.
  testthat::local_mocked_bindings(
    run_multi_seed_verification = function(...) stop("verification summary crashed"),
    .package = "BATON"
  )
  fit <- suppressMessages(bo_calibrate(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 6, q = 2, seed = 1, on_error = "return_partial",
    progress = FALSE, early_stop = list(enabled = FALSE),
    multi_seed_verify = TRUE, multi_seed_n = 3, multi_seed_reps = 200))
  expect_s3_class(fit, "BATON_fit")
  expect_equal(fit$status, "errored")
  expect_match(fit$error_message, "verification summary crashed")
  expect_equal(nrow(fit$history), 6L)
  expect_true(is.list(fit$best_theta))
  expect_true(is.finite(fit$best_objective))
  expect_null(fit$multi_seed_summary)
  expect_null(fit$multi_seed_runs)
  expect_true(is.na(fit$verdict))

  # Stop mode is unchanged: the same failure propagates as a plain error.
  err <- tryCatch(
    suppressMessages(bo_calibrate(
      fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 4, budget = 6, q = 2, seed = 1,
      progress = FALSE, early_stop = list(enabled = FALSE),
      multi_seed_verify = TRUE, multi_seed_n = 3, multi_seed_reps = 200)),
    error = conditionMessage)
  expect_identical(err, "verification summary crashed")
})

test_that("return_partial drops diagnostics with a warning when generation fails", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  testthat::local_mocked_bindings(
    generate_diagnostics = function(...) stop("posterior draw failure"),
    .package = "BATON"
  )
  # Under return_partial a diagnostics failure AFTER a successful final refit
  # must not lose the fit: diagnostics are dropped with a warning, everything
  # else (including the fitted surrogates) is returned, status untouched.
  expect_warning(
    fit <- suppressMessages(bo_calibrate(
      fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 4, budget = 6, q = 2, seed = 1, on_error = "return_partial",
      progress = FALSE, early_stop = list(enabled = FALSE))),
    "diagnostics")
  expect_equal(fit$status, "budget_exhausted")
  expect_null(fit$error_message)
  expect_null(fit$diagnostics)
  expect_false(is.null(fit$surrogates))
  expect_equal(nrow(fit$history), 6L)

  # Stop mode is unchanged: the same failure propagates as a plain error.
  err <- tryCatch(
    suppressMessages(bo_calibrate(
      fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
      n_init = 4, budget = 6, q = 2, seed = 1,
      progress = FALSE, early_stop = list(enabled = FALSE))),
    error = conditionMessage)
  expect_identical(err, "posterior draw failure")
})

test_that("return_partial on a clean run is equivalent to stop mode", {
  skip_if_not_installed("DiceKriging")
  fast_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 50 - 10 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  args <- list(
    fast_sim, list(x = c(0, 1)), "EN", list(power = c("ge", 0.8)),
    n_init = 4, budget = 8, q = 2, seed = 1,
    progress = FALSE, early_stop = list(enabled = FALSE))
  fit_stop <- suppressMessages(do.call(bo_calibrate, args))
  fit_rp <- suppressMessages(do.call(
    bo_calibrate, c(args, list(on_error = "return_partial"))))
  expect_equal(fit_rp$status, "budget_exhausted")
  expect_null(fit_rp$error_message)
  expect_equal(fit_rp$history, fit_stop$history)
  expect_equal(fit_rp$best_theta, fit_stop$best_theta)
  expect_equal(fit_rp$best_objective, fit_stop$best_objective)
})
