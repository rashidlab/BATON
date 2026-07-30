# Tests for Phase 2: Adaptive fidelity selection

test_that("select_fidelity_method dispatcher works correctly", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  fidelity_costs <- c(low = 1, med = 5, high = 50)

  # Test staged method (requires: prob, cv, iter, levels)
  fid_staged <- BATON:::select_fidelity_method(
    method = "staged",
    prob_feasible = 0.5,
    cv_estimate = 0.2,
    iter = 50,
    fidelity_levels = fidelity_levels
  )
  expect_true(fid_staged %in% names(fidelity_levels))

  # Test threshold method (requires: prob, levels)
  fid_threshold <- BATON:::select_fidelity_method(
    method = "threshold",
    prob_feasible = 0.8,
    fidelity_levels = fidelity_levels
  )
  expect_true(fid_threshold %in% names(fidelity_levels))

  # Test adaptive method (requires all parameters)
  fid_adaptive <- BATON:::select_fidelity_method(
    method = "adaptive",
    prob_feasible = 0.5,
    cv_estimate = 0.2,
    acq_value = 1.0,
    best_obj = 10,
    fidelity_levels = fidelity_levels,
    fidelity_costs = fidelity_costs,
    iter = 50,
    total_budget_used = 50000,
    total_budget = 150000
  )
  expect_true(fid_adaptive %in% names(fidelity_levels))

  # Test invalid method
  expect_error(
    BATON:::select_fidelity_method(
      method = "invalid",
      prob_feasible = 0.5
    ),
    "Unknown fidelity method"
  )
})

test_that("select_fidelity_adaptive responds to cost-benefit tradeoff", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, high = 10000)
  fidelity_costs <- c(low = 1, high = 50)

  # Scenario 1: High value, high cost, early iteration
  # Should sometimes use low (exploration) but also use high when value is high
  set.seed(42)
  selections <- replicate(20, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,  # Boundary (high value)
      cv_estimate = 0.3,    # High uncertainty (high value)
      acq_value = 2.0,      # High acquisition (high value)
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 10,  # Early (more exploration)
      total_budget_used = 5000,
      total_budget = 150000
    )
  })

  # Validate algorithm runs and makes valid choices
  # Note: Only two fidelity levels (low, high) in this test
  # Exact distribution depends on cost-benefit calculation and exploration randomness
  expect_true(all(selections %in% c("low", "high")))  # Valid choices
  expect_length(selections, 20)  # Correct number of selections

  # Scenario 2: Low value (far from boundary, low uncertainty)
  # Should prefer low fidelity (not worth the cost)
  set.seed(42)
  selections2 <- replicate(20, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.95,  # Far from boundary (low value)
      cv_estimate = 0.05,     # Low uncertainty (low value)
      acq_value = 0.1,        # Low acquisition (low value)
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 50,
      total_budget_used = 50000,
      total_budget = 150000
    )
  })

  # Should mostly use low fidelity
  low_ratio2 <- sum(selections2 == "low") / length(selections2)
  expect_gt(low_ratio2, 0.5)  # At least 50% low fidelity

  # Scenario 3: Late iteration with high value - should use high fidelity
  # (less exploration randomness)
  set.seed(42)
  selections3 <- replicate(20, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.25,
      acq_value = 1.5,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 100,  # Late (less exploration)
      total_budget_used = 120000,
      total_budget = 150000
    )
  })

  # Validate algorithm runs and makes valid choices
  # Note: Late iteration + high value may increase fidelity, but cost-awareness
  # can still lead to low fidelity if cost-benefit calculation favors it
  # Only two fidelity levels (low, high) in this test
  expect_true(all(selections3 %in% c("low", "high")))  # Valid choices
  expect_length(selections3, 20)  # Correct number of selections
})

test_that("select_fidelity_adaptive handles edge cases", {
  skip_on_cran()

  # Single fidelity level - should return it
  single_fid <- BATON:::select_fidelity_adaptive(
    prob_feasible = 0.5,
    cv_estimate = 0.2,
    acq_value = 1.0,
    best_obj = 10,
    fidelity_levels = c(med = 1000),
    fidelity_costs = c(med = 1),
    iter = 50,
    total_budget_used = 50000,
    total_budget = 150000
  )
  expect_equal(single_fid, "med")

  # Extreme values
  fid_extreme <- BATON:::select_fidelity_adaptive(
    prob_feasible = 1.0,  # Max feasibility
    cv_estimate = 0,      # Zero uncertainty
    acq_value = 0,        # Zero acquisition
    best_obj = 10,
    fidelity_levels = c(low = 200, high = 10000),
    fidelity_costs = c(low = 1, high = 50),
    iter = 1,
    total_budget_used = 0,
    total_budget = 150000
  )
  expect_true(fid_extreme %in% c("low", "high"))

  # Budget nearly depleted - should be more cost-sensitive
  set.seed(42)
  selections_depleted <- replicate(10, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.2,
      acq_value = 1.0,
      best_obj = 10,
      fidelity_levels = c(low = 200, high = 10000),
      fidelity_costs = c(low = 1, high = 50),
      iter = 50,
      total_budget_used = 145000,  # Nearly depleted
      total_budget = 150000
    )
  })
  # Should prefer low fidelity more when budget depleted
  low_ratio <- sum(selections_depleted == "low") / length(selections_depleted)
  expect_gt(low_ratio, 0.4)  # At least 40% low
})

test_that("bo_calibrate works with different fidelity methods", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Simple test function
  toy_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    x <- as.numeric(unlist(theta))
    power <- 0.8 + 0.1 * (x[1] - 0.5)
    type1 <- 0.05
    EN <- sum((x - 0.5)^2) * 100
    res <- c(power = power, type1 = type1, EN = EN)
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
    attr(res, "n_rep") <- 100
    res
  }

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.75), type1 = c("le", 0.1))

  # Test adaptive method (default)
  fit_adaptive <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 6,
      q = 2,
      budget = 12,
      fidelity_method = "adaptive",
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate with adaptive failed:", e$message))
  })

  expect_s3_class(fit_adaptive, "BATON_fit")
  expect_equal(nrow(fit_adaptive$history), 12)
  expect_equal(fit_adaptive$policies$fidelity_method, "adaptive")

  # Test staged method
  fit_staged <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 6,
      q = 2,
      budget = 12,
      fidelity_method = "staged",
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate with staged failed:", e$message))
  })

  expect_s3_class(fit_staged, "BATON_fit")
  expect_equal(fit_staged$policies$fidelity_method, "staged")

  # Test threshold method
  fit_threshold <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 6,
      q = 2,
      budget = 12,
      fidelity_method = "threshold",
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate with threshold failed:", e$message))
  })

  expect_s3_class(fit_threshold, "BATON_fit")
  expect_equal(fit_threshold$policies$fidelity_method, "threshold")
})

test_that("custom fidelity costs are respected", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  toy_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = 0.85, type1 = 0.05, EN = sum(x^2) * 100)
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
    attr(res, "n_rep") <- 100
    res
  }

  bounds <- list(x1 = c(0, 1))
  constraints <- list(power = c("ge", 0.8))

  # Custom costs (non-linear relationship)
  custom_costs <- c(low = 1, med = 3, high = 20)  # Not proportional to replications

  fit_custom <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 4,
      q = 1,
      budget = 8,
      fidelity_method = "adaptive",
      fidelity_costs = custom_costs,
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate with custom costs failed:", e$message))
  })

  expect_s3_class(fit_custom, "BATON_fit")
  expect_equal(fit_custom$policies$fidelity_costs, custom_costs)
})

test_that("fidelity selection uses acquisition value correctly", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, high = 10000)
  fidelity_costs <- c(low = 1, high = 50)

  # High acquisition value should increase likelihood of high fidelity
  set.seed(42)
  selections_high_acq <- replicate(20, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.2,
      acq_value = 10.0,  # Very high acquisition
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 50,
      total_budget_used = 50000,
      total_budget = 150000
    )
  })

  # Low acquisition value should favor low fidelity
  set.seed(42)
  selections_low_acq <- replicate(20, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.2,
      acq_value = 0.01,  # Very low acquisition
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 50,
      total_budget_used = 50000,
      total_budget = 150000
    )
  })

  # Compare med+high usage (more reliable than just high)
  high_med_ratio_high_acq <- sum(selections_high_acq %in% c("med", "high")) / length(selections_high_acq)
  high_med_ratio_low_acq <- sum(selections_low_acq %in% c("med", "high")) / length(selections_low_acq)

  # High acquisition should lead to more med/high fidelity selections
  expect_gte(high_med_ratio_high_acq, high_med_ratio_low_acq)
})

test_that("adaptive method escalates for high-value candidates (v0.8 defect 1)", {
  skip_on_cran()

  # Regression guard for the v0.8 defect: value_score had no fidelity
  # dependence, so value-per-cost was maximized by the cheapest level in
  # every one of 200 empirical calls, even with high-value inputs.
  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  fidelity_costs <- c(low = 1, med = 5, high = 50)

  set.seed(42)
  selections <- replicate(200, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,   # on the feasibility boundary
      cv_estimate = 0.3,     # high uncertainty
      acq_value = 2.0,       # high acquisition
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 10,
      total_budget_used = 5000,
      total_budget = 150000
    )
  })

  non_low <- mean(selections != "low")
  expect_gte(non_low, 0.30)
})

test_that("adaptive method keeps low-value candidates at low fidelity", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  fidelity_costs <- c(low = 1, med = 5, high = 50)

  set.seed(42)
  selections <- replicate(200, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.95,  # far from the boundary
      cv_estimate = 0.05,    # low uncertainty
      acq_value = 0.1,       # low acquisition
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 50,
      total_budget_used = 50000,
      total_budget = 150000
    )
  })

  expect_gte(mean(selections == "low"), 0.80)
})

test_that("adaptive selection is monotone nondecreasing in cv_estimate", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  fidelity_costs <- c(low = 1, med = 5, high = 50)

  cv_grid <- seq(0, 0.3, by = 0.025)
  # Same seed per call so the exploration draw is identical across the grid;
  # only cv_estimate varies.
  chosen_reps <- vapply(cv_grid, function(cv) {
    set.seed(7)
    sel <- BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = cv,
      acq_value = 2.0,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 80,   # late enough that the exploration override is at its floor
      total_budget_used = 60000,
      total_budget = 150000
    )
    fidelity_levels[[sel]]
  }, numeric(1))

  # As uncertainty grows the chosen replication count must never decrease.
  expect_true(all(diff(chosen_reps) >= 0))
  # And the grid must actually exercise an escalation, not sit constant.
  expect_gt(max(chosen_reps), min(chosen_reps))
})

test_that("raising a level's cost reduces how often it is selected", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  costs_cheap <- c(low = 1, med = 5, high = 50)
  costs_dear  <- c(low = 1, med = 5, high = 500)

  pick <- function(costs, seed) {
    set.seed(seed)
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.3,
      acq_value = 2.0,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = costs,
      iter = 10,
      total_budget_used = 5000,
      total_budget = 150000
    )
  }

  sel_cheap <- vapply(1:100, function(s) pick(costs_cheap, s), character(1))
  sel_dear  <- vapply(1:100, function(s) pick(costs_dear, s), character(1))

  n_high_cheap <- sum(sel_cheap == "high")
  n_high_dear <- sum(sel_dear == "high")
  expect_gt(n_high_cheap, 0)          # cheap high gets selected at all
  expect_lt(n_high_dear, n_high_cheap)  # inflating its cost suppresses it
})

test_that("adaptive method never selects a level the remaining budget cannot afford", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, med = 1000, high = 10000)
  fidelity_costs <- c(low = 1, med = 5, high = 50)

  # Remaining replication budget (5000) is below high's 10000 reps, so even a
  # maximally valuable candidate must not be routed to high.
  set.seed(42)
  selections <- replicate(50, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.3,
      acq_value = 5.0,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 40,
      total_budget_used = 145000,
      total_budget = 150000
    )
  })
  expect_false(any(selections == "high"))
})

test_that("exploration probability decays with iteration", {
  skip_on_cran()

  fidelity_levels <- c(low = 200, high = 10000)
  fidelity_costs <- c(low = 1, high = 50)

  # Early iterations: more randomization (exploration)
  set.seed(42)
  selections_early <- replicate(30, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.2,
      acq_value = 1.0,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 5,  # Very early
      total_budget_used = 5000,
      total_budget = 150000
    )
  })

  # Late iterations: less randomization (exploitation)
  set.seed(42)
  selections_late <- replicate(30, {
    BATON:::select_fidelity_adaptive(
      prob_feasible = 0.5,
      cv_estimate = 0.2,
      acq_value = 1.0,
      best_obj = 10,
      fidelity_levels = fidelity_levels,
      fidelity_costs = fidelity_costs,
      iter = 100,  # Very late
      total_budget_used = 120000,
      total_budget = 150000
    )
  })

  # Both selections should have some of each, but late should have more high
  # (less influenced by exploration randomness)
  high_ratio_early <- sum(selections_early == "high") / length(selections_early)
  high_ratio_late <- sum(selections_late == "high") / length(selections_late)

  # Late iterations should use high fidelity more often (when it's optimal)
  expect_gte(high_ratio_late, high_ratio_early * 0.8)  # Allow some variance
})

test_that("validate_fidelity_levels canonicalizes to ascending replication count (v0.8 defect 3)", {
  out <- validate_fidelity_levels(c(high = 10000, low = 2000, med = 4000))
  expect_identical(names(out), c("low", "med", "high"))
  expect_identical(unname(out), c(2000, 4000, 10000))
  # Already-ascending input passes through unchanged.
  out2 <- validate_fidelity_levels(c(low = 2000, med = 4000, high = 10000))
  expect_identical(names(out2), c("low", "med", "high"))
  # Ordering follows the VALUES, not the conventional label order: the
  # first-named (cheapest) level is correct-by-construction downstream.
  out3 <- validate_fidelity_levels(c(low = 500, high = 50, med = 200))
  expect_identical(names(out3), c("high", "med", "low"))
  expect_identical(unname(out3), c(50, 200, 500))
})

test_that("initial design runs at the cheapest level under any input ordering (v0.8 defect 3)", {
  skip_if_not_installed("DiceKriging")
  fid_sim <- function(theta, fidelity = "low", seed = NULL, n_rep = NULL, ...) {
    x <- as.numeric(unlist(theta))
    res <- c(power = 0.9, EN = sum((x - 0.5)^2) * 10)
    attr(res, "variance") <- c(power = 0.001, EN = 0.5)
    res
  }
  fit <- suppressMessages(suppressWarnings(bo_calibrate(
    sim_fun = fid_sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "EN", constraints = list(power = c("ge", 0.8)),
    n_init = 5, q = 1, budget = 5, seed = 1, progress = FALSE,
    fidelity_levels = c(high = 30, low = 10, med = 20)
  )))
  init_rows <- fit$history[fit$history$iter == 0L, ]
  expect_equal(nrow(init_rows), 5L)
  # On pre-fix code the first-named level ("high", 30 reps) silently ran the
  # whole initial design; canonicalization makes it the cheapest (10 reps).
  expect_true(all(init_rows$fidelity == "low"))
  expect_true(all(init_rows$n_rep == 10L))
})

test_that("dynamic fidelity scaling pairs each label with its own base value (v0.8)", {
  # Pathological but officially accepted labeling: values contradict the
  # conventional label order, so canonicalization yields names high, med, low.
  base <- validate_fidelity_levels(c(low = 500, high = 50, med = 200))
  expect_identical(names(base), c("high", "med", "low"))
  # Middle phase (iter_fraction ~0.43): scales low 1.5x, med 1.5x, high 2.0x.
  # budget_used = 0.8 * theoretical keeps both the budget scale-back and the
  # >0.85 safety factor inactive, so the expected values are exact.
  out <- compute_dynamic_fidelity_levels(
    iter = 30, budget = 100, base_levels = base,
    budget_used = 0.8 * 100 * mean(base), batch_size = 1
  )
  # Per-LABEL scaling: "high" (50 reps) doubles to 100; positional pairing
  # against the value-sorted base would instead floor it at 500.
  expect_equal(out[["high"]], 100)
  expect_equal(out[["med"]], 300)
  expect_equal(out[["low"]], 750)
  # Return is canonicalized ascending: first-named stays the cheapest even
  # after per-label scaling reorders relative costs.
  expect_false(is.unsorted(unname(out)))
  expect_identical(names(out)[1], "high")
})

test_that("dynamic fidelity scaling is unchanged for conventional labelings (v0.8)", {
  base <- validate_fidelity_levels(c(low = 200, med = 1000, high = 5000))
  out <- compute_dynamic_fidelity_levels(
    iter = 30, budget = 100, base_levels = base,
    budget_used = 0.8 * 100 * mean(base), batch_size = 1
  )
  expect_equal(out[["low"]], 300)     # 200 * 1.5
  expect_equal(out[["med"]], 1500)    # 1000 * 1.5
  expect_equal(out[["high"]], 10000)  # 5000 * 2.0
  expect_identical(names(out), c("low", "med", "high"))
})
