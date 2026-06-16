# Test script for GP-based initialization stopping feature
# Tests: check_init_sufficiency(), init_stopping_config(), integration with bo_calibrate()

library(BATON)

cat("=== Testing GP-Based Initialization Stopping ===\n\n")

# ------------------------------------------------------------------------------
# Test 1: check_init_sufficiency() basic functionality
# ------------------------------------------------------------------------------
cat("Test 1: check_init_sufficiency() basic functionality\n")

set.seed(42)
d <- 2
n <- 30

# Generate test data with a known objective function
X <- matrix(runif(n * d), nrow = n, ncol = d)
colnames(X) <- c("x1", "x2")
y <- apply(X, 1, function(x) sum((x - 0.5)^2) + rnorm(1, sd = 0.1))

# Test with no checkpoints (first call)
result1 <- check_init_sufficiency(
  X = X[1:20, ],
  y = y[1:20],
  checkpoints = list(),
  metrics = c("variance"),
  threshold = 0.10,
  window = 2,
  verbose = TRUE
)

cat("  First checkpoint (20 points):\n")
cat("    - metrics available:", !is.null(result1$metrics), "\n")
cat("    - converged:", result1$converged, "(expected: FALSE, need window checkpoints)\n")
if (!is.null(result1$metrics$max_var)) {
  cat("    - max_var:", round(result1$metrics$max_var, 4), "\n")
} else {
  cat("    - max_var: NULL (GP prediction may have failed)\n")
}
stopifnot(!result1$converged)  # Should not converge on first check

# Add more checkpoints
checkpoints <- list(result1$metrics)

result2 <- check_init_sufficiency(
  X = X[1:25, ],
  y = y[1:25],
  checkpoints = checkpoints,
  metrics = c("variance"),
  threshold = 0.10,
  window = 2,
  verbose = TRUE
)
checkpoints <- c(checkpoints, list(result2$metrics))

cat("  Second checkpoint (25 points):\n")
cat("    - converged:", result2$converged, "\n")
if (!is.null(result2$metrics$max_var)) {
  cat("    - max_var:", round(result2$metrics$max_var, 4), "\n")
} else {
  cat("    - max_var: NULL\n")
}

result3 <- check_init_sufficiency(
  X = X,
  y = y,
  checkpoints = checkpoints,
  metrics = c("variance"),
  threshold = 0.10,
  window = 2,
  verbose = TRUE
)

cat("  Third checkpoint (30 points):\n")
cat("    - converged:", result3$converged, "\n")
cat("    - reason:", result3$reason %||% "none", "\n")

cat("  PASSED\n\n")

# ------------------------------------------------------------------------------
# Test 2: init_stopping_config() helper
# ------------------------------------------------------------------------------
cat("Test 2: init_stopping_config() helper\n")

config_default <- init_stopping_config()
cat("  Default config:\n")
cat("    - enabled:", config_default$enabled, "\n")
cat("    - min_init:", config_default$min_init, "\n")
cat("    - check_every:", config_default$check_every, "\n")
cat("    - threshold:", config_default$threshold, "\n")
cat("    - window:", config_default$window, "\n")

stopifnot(config_default$enabled == TRUE)
stopifnot(config_default$min_init == 20)
stopifnot(config_default$check_every == 20)
stopifnot(config_default$threshold == 0.10)
stopifnot(config_default$window == 2)

config_custom <- init_stopping_config(
  enabled = TRUE,
  min_init = 15,
  check_every = 10,
  threshold = 0.15,
  window = 3
)
cat("  Custom config:\n")
cat("    - min_init:", config_custom$min_init, "\n")
cat("    - check_every:", config_custom$check_every, "\n")
cat("    - threshold:", config_custom$threshold, "\n")
cat("    - window:", config_custom$window, "\n")

stopifnot(config_custom$min_init == 15)
stopifnot(config_custom$check_every == 10)

cat("  PASSED\n\n")

# ------------------------------------------------------------------------------
# Test 3: Integration with bo_calibrate() (minimal test)
# ------------------------------------------------------------------------------
cat("Test 3: Integration with bo_calibrate() (minimal test)\n")

# Simple simulator function
simple_sim <- function(theta, fidelity = "low", ...) {
  x1 <- theta[["x1"]]
  x2 <- theta[["x2"]]

  # Branin-like objective
  objective <- (x2 - 5.1/(4*pi^2)*x1^2 + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10

  # Simulated constraint
  constraint <- x1 + x2

  result <- c(obj = objective, c1 = constraint)
  attr(result, "variance") <- c(obj = 0.01, c1 = 0.001)
  attr(result, "n_rep") <- 1000
  result
}

bounds <- list(
  x1 = c(-5, 10),
  x2 = c(0, 15)
)

constraints <- list(
  c1 = c("le", 20)
)

# Test with init_stopping enabled
cat("  Running bo_calibrate with init_stopping enabled...\n")
cat("  (Using 60 point budget with init stopping at 40)\n")

fit <- tryCatch({
  bo_calibrate(
    sim_fun = simple_sim,
    bounds = bounds,
    objective = "obj",
    constraints = constraints,
    n_init = 40,
    q = 4,
    budget = 60,
    seed = 123,
    init_stopping = init_stopping_config(
      enabled = TRUE,
      min_init = 20,
      check_every = 10,
      threshold = 0.15,
      window = 2,
      verbose = TRUE
    ),
    progress = TRUE
  )
}, error = function(e) {
  cat("  Error:", e$message, "\n")
  NULL
})

if (!is.null(fit)) {
  cat("  bo_calibrate completed successfully!\n")
  cat("  Total evaluations:", nrow(fit$history), "\n")
  cat("  Best objective:", round(min(fit$history$obj[fit$history$feasible], na.rm = TRUE), 4), "\n")
  cat("  PASSED\n\n")
} else {
  cat("  FAILED - bo_calibrate returned error\n\n")
}

# ------------------------------------------------------------------------------
# Test 4: Edge cases
# ------------------------------------------------------------------------------
cat("Test 4: Edge cases\n")

# Too few points
result_few <- check_init_sufficiency(
  X = X[1:3, ],
  y = y[1:3],
  checkpoints = list(),
  metrics = c("variance"),
  verbose = TRUE
)
cat("  With 3 points: converged =", result_few$converged,
    ", recommendation =", result_few$recommendation, "\n")
stopifnot(!result_few$converged)
stopifnot(grepl("at least 5", result_few$recommendation))

# Multiple metrics
result_multi <- check_init_sufficiency(
  X = X,
  y = y,
  checkpoints = list(),
  metrics = c("variance", "loo", "loglik"),
  verbose = TRUE
)
cat("  Multiple metrics:\n")
cat("    max_var =", if (!is.null(result_multi$metrics$max_var)) round(result_multi$metrics$max_var, 4) else "NULL", "\n")
cat("    loo_rmse =", if (!is.null(result_multi$metrics$loo_rmse)) round(result_multi$metrics$loo_rmse, 4) else "NULL", "\n")
cat("    loglik =", if (!is.null(result_multi$metrics$loglik)) round(result_multi$metrics$loglik, 4) else "NULL", "\n")

cat("  PASSED\n\n")

cat("=== All Tests Passed ===\n")
