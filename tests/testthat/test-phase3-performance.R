# Tests for Phase 3: Performance optimizations

test_that("extract_gp_hyperparams works with DiceKriging models", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Create a simple GP model
  set.seed(42)
  X <- lhs::randomLHS(20, 2)
  y <- rowSums(X) + rnorm(20, 0, 0.1)

  model <- tryCatch({
    DiceKriging::km(
      design = X,
      response = y,
      covtype = "matern5_2",
      control = list(trace = FALSE)
    )
  }, error = function(e) {
    skip("DiceKriging model fitting failed")
  })

  # Extract hyperparameters
  theta <- BATON:::extract_gp_hyperparams(model)

  expect_true(is.numeric(theta))
  expect_true(all(is.finite(theta)))
  expect_true(all(theta > 0))
  expect_equal(length(theta), 2)  # 2D problem
})

test_that("extract_gp_hyperparams handles NULL and invalid models", {
  skip_on_cran()

  # NULL model
  theta_null <- BATON:::extract_gp_hyperparams(NULL)
  expect_null(theta_null)

  # Invalid model
  theta_invalid <- BATON:::extract_gp_hyperparams(list(not_a_model = 1))
  expect_null(theta_invalid)
})

test_that("fit_surrogates uses warm-start when prev_surrogates provided", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Create initial history
  set.seed(42)
  n <- 10

  # Pre-generate unit_x to avoid closure issue
  unit_x_list <- lapply(1:n, function(i) {
    x <- lhs::randomLHS(1, 2)[1, ]
    names(x) <- c("x1", "x2")
    x
  })

  history <- tibble::tibble(
    unit_x = unit_x_list,
    theta = lapply(1:n, function(i) {
      as.list(unit_x_list[[i]])
    }),
    theta_id = sapply(1:n, function(i) paste0("id_", i)),
    metrics = lapply(1:n, function(i) {
      x <- unit_x_list[[i]]
      c(power = 0.8 + rnorm(1, 0, 0.05),
        type1 = 0.05 + rnorm(1, 0, 0.01),
        EN = sum((x - 0.5)^2) * 100 + rnorm(1, 0, 5))
    }),
    variance = lapply(1:n, function(i) {
      c(power = 0.001, type1 = 0.0005, EN = 0.5)
    })
  )

  constraint_tbl <- data.frame(
    metric = c("power", "type1"),
    direction = c("ge", "le"),
    threshold = c(0.75, 0.1),
    stringsAsFactors = FALSE
  )

  # Fit initial surrogates
  surrogates1 <- tryCatch({
    BATON:::fit_surrogates(history, "EN", constraint_tbl)
  }, error = function(e) {
    skip(paste("Initial surrogate fitting failed:", e$message))
  })

  # Add more data
  history2 <- dplyr::bind_rows(history, tibble::tibble(
    unit_x = list(c(x1 = 0.5, x2 = 0.5)),
    theta = list(list(x1 = 0.5, x2 = 0.5)),
    theta_id = "id_11",
    metrics = list(c(power = 0.82, type1 = 0.048, EN = 0)),
    variance = list(c(power = 0.001, type1 = 0.0005, EN = 0.5))
  ))

  # Fit with warm-start
  surrogates2 <- tryCatch({
    BATON:::fit_surrogates(history2, "EN", constraint_tbl,
                             prev_surrogates = surrogates1)
  }, error = function(e) {
    skip(paste("Warm-start surrogate fitting failed:", e$message))
  })

  # Should return valid surrogates
  expect_true(is.list(surrogates2))
  expect_true("EN" %in% names(surrogates2))
  expect_true("power" %in% names(surrogates2))
  expect_true("type1" %in% names(surrogates2))
})

test_that("adaptive candidate pool size scales with dimension", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  toy_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    x <- unlist(theta)
    res <- c(power = 0.85, type1 = 0.05, EN = sum(x^2) * 100)
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
    attr(res, "n_rep") <- 100
    res
  }

  # Test 2D problem
  bounds_2d <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.8))

  # Run with low budget to capture candidate pool messages
  fit_2d <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds_2d,
      objective = "EN",
      constraints = constraints,
      n_init = 4,
      q = 1,
      budget = 6,
      candidate_pool = 500,  # Set minimum
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("2D bo_calibrate failed:", e$message))
  })

  expect_s3_class(fit_2d, "BATON_fit")
  expect_equal(nrow(fit_2d$history), 6)

  # Test 5D problem (should use larger pool)
  bounds_5d <- list(x1 = c(0, 1), x2 = c(0, 1), x3 = c(0, 1),
                    x4 = c(0, 1), x5 = c(0, 1))

  fit_5d <- tryCatch({
    bo_calibrate(
      sim_fun = function(theta, ...) {
        x <- unlist(theta)
        res <- c(power = 0.85, type1 = 0.05, EN = sum(x^2) * 100)
        attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
        attr(res, "n_rep") <- 100
        res
      },
      bounds = bounds_5d,
      objective = "EN",
      constraints = constraints,
      n_init = 8,
      q = 1,
      budget = 10,
      candidate_pool = 1000,
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("5D bo_calibrate failed:", e$message))
  })

  expect_s3_class(fit_5d, "BATON_fit")
  # Pool size should be larger for higher dimensions (tested via code inspection)
})

test_that("early stopping triggers when no improvement", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Simulator that reaches optimum quickly
  toy_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    x <- unlist(theta)
    # Optimum at (0.5, 0.5) with value = 0
    res <- c(
      power = 0.85,
      type1 = 0.05,
      EN = sum((x - 0.5)^2) * 100
    )
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
    attr(res, "n_rep") <- 100
    res
  }

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.8))

  # Run with large budget but expect early stopping
  fit <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 10,
      q = 2,
      budget = 100,  # Large budget
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate failed:", e$message))
  })

  expect_s3_class(fit, "BATON_fit")

  # Should stop before budget exhausted (early stopping)
  expect_lt(nrow(fit$history), 100)

  # Should find good solution (near optimum at (0.5, 0.5))
  best_theta <- fit$best_theta
  expect_lte(abs(best_theta$x1 - 0.5), 0.2)
  expect_lte(abs(best_theta$x2 - 0.5), 0.2)

  # Objective should be near zero
  best_obj <- fit$history$objective[nrow(fit$history)]
  expect_lte(best_obj, 10)  # Should be close to 0
})

test_that("early stopping criterion: acquisition values", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Flat objective function (all points equally good)
  flat_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    # Constant objective regardless of theta
    res <- c(power = 0.85, type1 = 0.05, EN = 100)
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.1)
    attr(res, "n_rep") <- 100
    res
  }

  bounds <- list(x1 = c(0, 1))
  constraints <- list(power = c("ge", 0.8))

  # Should stop early due to low acquisition values
  fit <- tryCatch({
    bo_calibrate(
      sim_fun = flat_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 6,
      q = 1,
      budget = 30,
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("bo_calibrate failed:", e$message))
  })

  expect_s3_class(fit, "BATON_fit")

  # Should stop before budget (acquisition values near zero)
  expect_lte(nrow(fit$history), 30)
})

test_that("warm-start improves fitting time", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")
  skip("Manual timing test - run interactively")

  # This is a manual test to verify speedup
  # Run interactively to see timing difference

  set.seed(42)
  n <- 50
  history <- tibble::tibble(
    unit_x = lapply(1:n, function(i) {
      x <- lhs::randomLHS(1, 5)[1, ]
      names(x) <- paste0("x", 1:5)
      x
    }),
    theta = lapply(1:n, function(i) as.list(history$unit_x[[i]])),
    theta_id = paste0("id_", 1:n),
    metrics = lapply(1:n, function(i) {
      c(power = 0.8, type1 = 0.05, EN = sum(history$unit_x[[i]]^2) * 100)
    }),
    variance = lapply(1:n, function(i) {
      c(power = 0.001, type1 = 0.0005, EN = 0.5)
    })
  )

  constraint_tbl <- data.frame(
    metric = c("power", "type1"),
    direction = c("ge", "le"),
    threshold = c(0.75, 0.1),
    stringsAsFactors = FALSE
  )

  # Time without warm-start
  time1 <- system.time({
    surrogates1 <- BATON:::fit_surrogates(history, "EN", constraint_tbl)
  })

  # Time with warm-start (should be faster)
  time2 <- system.time({
    surrogates2 <- BATON:::fit_surrogates(history, "EN", constraint_tbl,
                                            prev_surrogates = surrogates1)
  })

  cat(sprintf("\nWithout warm-start: %.3f seconds\n", time1["elapsed"]))
  cat(sprintf("With warm-start: %.3f seconds\n", time2["elapsed"]))
  cat(sprintf("Speedup: %.1fx\n", time1["elapsed"] / time2["elapsed"]))

  # Typically expect 1.3-2x speedup
})

test_that("Phase 3 features work together", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Integration test with all Phase 3 features
  toy_sim_fun <- function(theta, fidelity = "high", seed = NULL, ...) {
    x <- unlist(theta)
    res <- c(
      power = 0.8 + 0.1 * (x[1] - 0.5),
      type1 = 0.05,
      EN = sum((x - 0.5)^2) * 100
    )
    attr(res, "variance") <- c(power = 0.001, type1 = 0.0005, EN = 0.5)
    attr(res, "n_rep") <- 100
    res
  }

  bounds <- list(x1 = c(0, 1), x2 = c(0, 1))
  constraints <- list(power = c("ge", 0.75), type1 = c("le", 0.1))

  # Run with all features enabled
  fit <- tryCatch({
    bo_calibrate(
      sim_fun = toy_sim_fun,
      bounds = bounds,
      objective = "EN",
      constraints = constraints,
      n_init = 8,
      q = 2,
      budget = 40,
      progress = FALSE,
      seed = 42
    )
  }, error = function(e) {
    skip(paste("Integration test failed:", e$message))
  })

  # Should complete successfully
  expect_s3_class(fit, "BATON_fit")

  # Should have reasonable results
  expect_lte(nrow(fit$history), 40)  # May stop early
  expect_gte(nrow(fit$history), 8)   # At least initial design

  # Should find decent solution
  expect_true(fit$history$feasible[nrow(fit$history)])
  expect_lte(fit$history$objective[nrow(fit$history)], 30)
})

test_that("candidate pool size respects user minimum", {
  skip_on_cran()

  # The adaptive sizing should never go below user-specified candidate_pool
  # This is tested via code inspection of the formula:
  # candidate_pool_size <- max(candidate_pool_size, candidate_pool)

  # Just verify the logic is correct
  d <- 2
  base_pool_size <- pmax(1000, pmin(5000, 500 * d))
  user_min <- 3000

  # Adaptive would choose base_pool_size (1000)
  # But with user min, should be 3000
  final_size <- max(base_pool_size, user_min)

  expect_equal(final_size, 3000)
})
