# Task B2: warmstart_from / initial_history donor seeds AUGMENT the LHS
# initial design instead of replacing it (PI decision 2026-07-22).

ws_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  x <- theta$x
  v <- c(power = 0.85 + 0.1 * x + rnorm(1, sd = 0.005),
         EN = 50 - 10 * x + rnorm(1, sd = 0.05))
  attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
  v
}

test_that("warmstart_from donors are augmented with fresh LHS up to n_init", {
  donors <- list(
    list(theta = list(x = 0.123), power = 0.87, EN = 48.8),
    list(theta = list(x = 0.456), power = 0.90, EN = 45.4)
  )
  fit <- suppressMessages(bo_calibrate(
    sim_fun = ws_sim,
    bounds = list(x = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    warmstart_from = donors,
    n_init = 8,
    budget = 10,
    seed = 7,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  ))
  init_rows <- fit$history[fit$history$iter == 0, ]
  # 2 donor rows + 6 fresh LHS points = n_init
  expect_equal(nrow(init_rows), 8)
  # donor thetas present among the initial rows
  donor_ids <- vapply(donors, function(d)
    BATON:::theta_to_id(BATON:::scale_to_unit(d$theta, list(x = c(0, 1)))),
    character(1))
  expect_true(all(donor_ids %in% init_rows$theta_id))
  # no duplicate theta_ids in the initial design
  expect_false(any(duplicated(init_rows$theta_id)))
  # budget still respected overall
  expect_lte(nrow(fit$history), 10)
})

test_that("augment with integer_params neither duplicates nor silently underfills", {
  # bounds c(1, 3) with integer coercion admits only 3 distinct designs.
  # Donor covers x=2; the top-up can add at most x=1 and x=3. It must not
  # evaluate coerced duplicates at iter 0, and must stop cleanly when the
  # discrete space is exhausted (3 rows, not n_init = 6).
  int_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 40 + theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  # suppressWarnings: km fits on a 3-point discrete space legitimately fall
  # back (L-BFGS-B non-finite fn); the fallback path is what handles it.
  fit <- suppressWarnings(suppressMessages(bo_calibrate(
    sim_fun = int_sim,
    bounds = list(x = c(1, 3)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    warmstart_from = list(list(theta = list(x = 2), power = 0.95, EN = 42)),
    integer_params = "x",
    n_init = 6,
    budget = 8,
    seed = 7,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  )))
  init_rows <- fit$history[!is.na(fit$history$iter) & fit$history$iter == 0, ]
  # no duplicate designs in the initial design
  expect_false(any(duplicated(init_rows$theta_id)))
  # all three distinct integer designs present, and nothing beyond them
  expect_equal(nrow(init_rows), 3)
})

test_that("top-up fills n_init on a large integer grid despite collisions", {
  # 100 distinct designs, 1 donor, n_init = 60: rounded LHS draws collide with
  # already-seen designs, so the top-up must keep replenishing until n_init is
  # reached, not stop after a fixed number of rounds.
  int_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 40 + 0.1 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressWarnings(suppressMessages(bo_calibrate(
    sim_fun = int_sim,
    bounds = list(x = c(1, 100)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    warmstart_from = list(list(theta = list(x = 50), power = 0.95, EN = 45)),
    integer_params = "x",
    n_init = 60,
    budget = 60,
    seed = 1,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  )))
  init_rows <- fit$history[!is.na(fit$history$iter) & fit$history$iter == 0, ]
  expect_equal(nrow(init_rows), 60)
  expect_false(any(duplicated(init_rows$theta_id)))
})

test_that("exhaustion count matches round()'s reachable set on non-integral bounds", {
  # round() over [0.4, 3.6] reaches {0, 1, 2, 3, 4}: values in [0.4, 0.5)
  # round to 0 and values in (3.5, 3.6] round to 4. An exhaustion count of
  # ceiling(lower):floor(upper) (= 3 designs) would stop the top-up early.
  int_sim <- function(theta, fidelity = "low", seed = NULL, ...) {
    v <- c(power = 0.95, EN = 40 + 0.1 * theta$x)
    attr(v, "variance") <- c(power = 1e-4, EN = 0.01)
    v
  }
  fit <- suppressWarnings(suppressMessages(bo_calibrate(
    sim_fun = int_sim,
    bounds = list(x = c(0.4, 3.6)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    warmstart_from = list(list(theta = list(x = 2), power = 0.95, EN = 40.2)),
    integer_params = "x",
    n_init = 5,
    budget = 5,
    seed = 3,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  )))
  init_rows <- fit$history[!is.na(fit$history$iter) & fit$history$iter == 0, ]
  expect_equal(nrow(init_rows), 5)
  expect_false(any(duplicated(init_rows$theta_id)))
})

test_that("initial_history with >= n_init rows is used as-is (no augmentation)", {
  set.seed(11)
  xs <- seq(0.05, 0.95, length.out = 5)
  ih <- data.frame(
    x = xs,
    power = 0.85 + 0.1 * xs,
    EN = 50 - 10 * xs,
    objective = 50 - 10 * xs,
    fidelity = "high",
    feasible = TRUE
  )
  fit <- suppressMessages(bo_calibrate(
    sim_fun = ws_sim,
    bounds = list(x = c(0, 1)),
    objective = "EN",
    constraints = list(power = c("ge", 0.8)),
    initial_history = ih,
    n_init = 3,
    budget = 7,
    seed = 7,
    progress = FALSE,
    early_stop = list(enabled = FALSE)
  ))
  # all 5 provided rows kept, none added at iter 0
  init_rows <- fit$history[fit$history$iter %in% c(0L, NA), ]
  expect_equal(nrow(init_rows), 5)
})
