#' Compute mean and variance incrementally using Welford's algorithm
#'
#' This function provides a memory-efficient way to compute sample means and
#' variances when you have a function that generates individual samples. Instead
#' of storing all samples in memory, it updates running statistics incrementally.
#'
#' @param sample_fn A function with signature `function(i, ...)` that returns a
#'   named numeric vector for the i-th sample. The function receives the
#'   iteration index and any additional arguments passed via `...`.
#' @param n_samples Number of samples to generate.
#' @param ... Additional arguments passed to `sample_fn`.
#'
#' @return A list with two components:
#'   \describe{
#'     \item{mean}{Named numeric vector of sample means}
#'     \item{variance}{Named numeric vector of sample variances (variance of the mean)}
#'   }
#'
#' @details
#' Welford's algorithm computes mean and variance in a single pass without
#' storing all samples. This is numerically stable and memory-efficient.
#'
#' The variance returned is the variance of the sample mean (i.e., divided by n),
#' not the sample variance. This is appropriate for Monte Carlo simulation where
#' we want to quantify uncertainty in the mean estimate.
#'
#' Memory usage: O(m) where m is the number of metrics, regardless of n_samples.
#'
#' @references
#' Welford, B. P. (1962). "Note on a method for calculating corrected sums of
#' squares and products". Technometrics. 4 (3): 419–420.
#'
#' @examples
#' \dontrun{
#' # Example: Simulate a clinical trial outcome
#' simulate_trial <- function(i, theta) {
#'   # Simulate one trial with given parameters
#'   success <- rbinom(1, 1, theta$prob)
#'   sample_size <- rpois(1, theta$expected_n)
#'   c(power = success, EN = sample_size)
#' }
#'
#' # Compute mean and variance over 1000 simulations
#' result <- welford_mean_var(
#'   sample_fn = simulate_trial,
#'   n_samples = 1000,
#'   theta = list(prob = 0.8, expected_n = 100)
#' )
#'
#' # Use in your simulator:
#' sim_fun <- function(theta, fidelity = "high", ...) {
#'   n_rep <- switch(fidelity, low = 200, med = 1000, high = 10000)
#'
#'   result <- welford_mean_var(
#'     sample_fn = function(i, theta) {
#'       trial <- run_one_trial(theta)
#'       c(power = trial$power, EN = trial$sample_size)
#'     },
#'     n_samples = n_rep,
#'     theta = theta
#'   )
#'
#'   metrics <- result$mean
#'   attr(metrics, "variance") <- result$variance
#'   attr(metrics, "n_rep") <- n_rep
#'   return(metrics)
#' }
#' }
#'
#' @export
welford_mean_var <- function(sample_fn, n_samples, ...) {
  if (n_samples <= 0) {
    stop("`n_samples` must be positive.", call. = FALSE)
  }

  # Generate first sample to determine dimension
  first_sample <- sample_fn(1, ...)
  if (!is.numeric(first_sample)) {
    stop("`sample_fn` must return a numeric vector.", call. = FALSE)
  }

  n_metrics <- length(first_sample)
  metric_names <- names(first_sample)

  # Initialize accumulators
  mean_vec <- as.numeric(first_sample)
  M2_vec <- numeric(n_metrics)

  # Process remaining samples with Welford's algorithm
  if (n_samples > 1) {
    for (i in 2:n_samples) {
      x <- sample_fn(i, ...)

      if (length(x) != n_metrics) {
        stop(sprintf("Sample %d has different length (%d) than first sample (%d).",
                     i, length(x), n_metrics), call. = FALSE)
      }

      # Welford's online update
      delta <- as.numeric(x) - mean_vec
      mean_vec <- mean_vec + delta / i
      delta2 <- as.numeric(x) - mean_vec
      M2_vec <- M2_vec + delta * delta2
    }
  }

  # Compute variance of the mean
  if (n_samples > 1) {
    variance <- M2_vec / (n_samples * (n_samples - 1))
  } else {
    variance <- rep(NA_real_, n_metrics)
  }

  # Restore names
  names(mean_vec) <- metric_names
  names(variance) <- metric_names

  list(
    mean = mean_vec,
    variance = variance,
    n = n_samples
  )
}


#' Compute mean and variance from parallel simulation chunks
#'
#' When running simulations in parallel, each worker computes summary statistics
#' independently. This function pools results from multiple workers, correctly
#' accounting for both within-chunk and between-chunk variance.
#'
#' @param chunk_results A list where each element is a list with components
#'   `mean` (numeric vector), `M2` (sum of squared deviations), and `n`
#'   (number of samples). Typically produced by calling [welford_mean_var()]
#'   on each parallel worker.
#'
#' @return A list with components:
#'   \describe{
#'     \item{mean}{Pooled mean across all chunks}
#'     \item{variance}{Pooled variance of the mean}
#'     \item{n}{Total number of samples}
#'   }
#'
#' @details
#' Uses Chan's parallel variance algorithm to combine statistics from
#' independent chunks. This correctly accounts for variance both within and
#' between chunks.
#'
#' @references
#' Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1979), "Updating
#' Formulae and a Pairwise Algorithm for Computing Sample Variances"
#'
#' @examples
#' \dontrun{
#' # Parallel simulation example
#' library(parallel)
#' n_cores <- 4
#' n_rep <- 10000
#'
#' # Run simulations in parallel chunks
#' chunk_results <- mclapply(1:n_cores, function(core_id) {
#'   chunk_size <- n_rep / n_cores
#'   welford_mean_var(
#'     sample_fn = function(i, theta) run_trial(theta),
#'     n_samples = chunk_size,
#'     theta = theta
#'   )
#' }, mc.cores = n_cores)
#'
#' # Pool results
#' pooled <- pool_welford_results(chunk_results)
#'
#' # Use in simulator
#' metrics <- pooled$mean
#' attr(metrics, "variance") <- pooled$variance
#' attr(metrics, "n_rep") <- pooled$n
#' }
#'
#' @export
pool_welford_results <- function(chunk_results) {
  if (length(chunk_results) == 0) {
    stop("No chunk results to pool.", call. = FALSE)
  }

  if (length(chunk_results) == 1) {
    return(chunk_results[[1]])
  }

  # Start with first chunk
  pooled_mean <- chunk_results[[1]]$mean
  pooled_n <- chunk_results[[1]]$n
  pooled_M2 <- chunk_results[[1]]$M2 %||%
    (chunk_results[[1]]$variance * pooled_n * (pooled_n - 1))

  # Sequentially combine remaining chunks using Chan's algorithm
  for (i in 2:length(chunk_results)) {
    chunk <- chunk_results[[i]]
    chunk_M2 <- chunk$M2 %||% (chunk$variance * chunk$n * (chunk$n - 1))

    # Chan's parallel variance formula
    delta <- chunk$mean - pooled_mean
    total_n <- pooled_n + chunk$n

    # Update pooled statistics
    pooled_mean <- (pooled_n * pooled_mean + chunk$n * chunk$mean) / total_n
    pooled_M2 <- pooled_M2 + chunk_M2 +
      (delta^2 * pooled_n * chunk$n / total_n)
    pooled_n <- total_n
  }

  # Compute variance of the mean
  variance <- pooled_M2 / (pooled_n * (pooled_n - 1))

  list(
    mean = pooled_mean,
    variance = variance,
    n = pooled_n,
    M2 = pooled_M2
  )
}
