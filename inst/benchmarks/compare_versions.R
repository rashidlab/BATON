# Benchmark: Compare v0.2.0 vs v0.3.0 Behavior
#
# This script compares the performance of BATON v0.2.0 (baseline) vs v0.3.0
# (with all improvements) across multiple test problems.
#
# v0.2.0 equivalent configuration:
#   - fidelity_method = "staged"
#   - q = 1 (no batch diversity)
#   - No early stopping (runs full budget)
#   - No warm-start (cold start each iteration)
#
# v0.3.0 configuration:
#   - fidelity_method = "adaptive"
#   - q = 4 (batch diversity via local penalization)
#   - Early stopping enabled (automatic)
#   - Warm-start enabled (automatic)
#   - Adaptive candidate pool (automatic)
#
# Expected improvements (v0.3.0 vs v0.2.0):
#   - 10-20% fewer evaluations to convergence (batch diversity + early stop)
#   - 30-60% faster wall-clock time (warm-start + adaptive pool + early stop)
#   - Similar or better final objective values (no regression)

library(BATON)
library(tictoc)
library(dplyr)
library(ggplot2)

# Load test problems
source("inst/benchmarks/test_problems.R")

#' Run Single Benchmark
#'
#' Compares v0.2.0 vs v0.3.0 on a single problem and seed.
#'
#' @param problem problem configuration
#' @param seed RNG seed
#' @param n_init initial design size
#' @param budget total evaluation budget
#' @return list with comparison results
run_single_benchmark <- function(problem, seed = 42, n_init = 10, budget = 50) {
  cat("\n=== Benchmarking:", problem$name, "(seed =", seed, ") ===\n")

  # v0.2.0 equivalent: staged, no batch diversity, full budget
  cat("Running v0.2.0 equivalent...\n")
  tic()
  fit_v02 <- bo_calibrate(
    sim_fun = problem$sim_fun,
    bounds = problem$bounds,
    objective = problem$objective,
    constraints = problem$constraints,
    n_init = n_init,
    q = 1,  # No batch diversity
    budget = budget,
    fidelity_method = "staged",
    seed = seed,
    progress = FALSE
  )
  time_v02 <- toc(quiet = TRUE)

  # v0.3.0: adaptive, batch diversity, early stopping
  cat("Running v0.3.0...\n")
  tic()
  fit_v03 <- bo_calibrate(
    sim_fun = problem$sim_fun,
    bounds = problem$bounds,
    objective = problem$objective,
    constraints = problem$constraints,
    n_init = n_init,
    q = 4,  # Batch diversity
    budget = budget,
    fidelity_method = "adaptive",
    seed = seed,
    progress = FALSE
  )
  time_v03 <- toc(quiet = TRUE)

  # Extract metrics
  feasible_v02 <- fit_v02$history %>% filter(feasible)
  feasible_v03 <- fit_v03$history %>% filter(feasible)

  results <- list(
    problem = problem$name,
    seed = seed,

    # v0.2.0 metrics
    v02_time = time_v02$toc - time_v02$tic,
    v02_evals = nrow(fit_v02$history),
    v02_best = if (nrow(feasible_v02) > 0) min(feasible_v02$objective) else NA,
    v02_feasible_rate = mean(fit_v02$history$feasible),
    v02_fidelity_dist = table(fit_v02$history$fidelity),

    # v0.3.0 metrics
    v03_time = time_v03$toc - time_v03$tic,
    v03_evals = nrow(fit_v03$history),
    v03_best = if (nrow(feasible_v03) > 0) min(feasible_v03$objective) else NA,
    v03_feasible_rate = mean(fit_v03$history$feasible),
    v03_fidelity_dist = table(fit_v03$history$fidelity),
    v03_early_stop = nrow(fit_v03$history) < budget,

    # Improvements
    time_speedup = (time_v02$toc - time_v02$tic) / (time_v03$toc - time_v03$tic),
    eval_reduction = (nrow(fit_v02$history) - nrow(fit_v03$history)) / nrow(fit_v02$history),
    objective_improvement = if (!is.na(fit_v02$best_obj) && !is.na(fit_v03$best_obj)) {
      (fit_v02$best_obj - fit_v03$best_obj) / abs(fit_v02$best_obj)
    } else {
      NA
    },

    # Full fit objects (for detailed analysis)
    fit_v02 = fit_v02,
    fit_v03 = fit_v03
  )

  # Print summary
  cat(sprintf("\nResults for %s (seed %d):\n", problem$name, seed))
  cat(sprintf("  v0.2.0: %d evals, %.1f sec, best = %.2f\n",
              results$v02_evals, results$v02_time, results$v02_best))
  cat(sprintf("  v0.3.0: %d evals, %.1f sec, best = %.2f\n",
              results$v03_evals, results$v03_time, results$v03_best))
  cat(sprintf("  Speedup: %.2fx faster\n", results$time_speedup))
  cat(sprintf("  Eval reduction: %.1f%%\n", results$eval_reduction * 100))
  cat(sprintf("  Early stop: %s\n", ifelse(results$v03_early_stop, "YES", "NO")))

  return(results)
}


#' Run Complete Benchmark Suite
#'
#' Compares v0.2.0 vs v0.3.0 across all test problems and multiple seeds.
#'
#' @param problems list of problem configurations
#' @param seeds vector of RNG seeds
#' @param n_init initial design size
#' @param budget total evaluation budget
#' @param save_results whether to save results to RDS file
#' @return data frame of results
run_benchmark_suite <- function(problems = get_test_problems(),
                                 seeds = 1:10,
                                 n_init = 10,
                                 budget = 50,
                                 save_results = TRUE) {

  cat("\n======================================\n")
  cat("  BATON v0.2.0 vs v0.3.0 Benchmark\n")
  cat("======================================\n")
  cat(sprintf("Problems: %d\n", length(problems)))
  cat(sprintf("Seeds: %d\n", length(seeds)))
  cat(sprintf("Budget: %d evaluations per run\n", budget))
  cat(sprintf("Total runs: %d\n", length(problems) * length(seeds) * 2))
  cat("======================================\n\n")

  all_results <- list()
  idx <- 1

  for (prob in problems) {
    for (seed in seeds) {
      result <- run_single_benchmark(prob, seed = seed, n_init = n_init, budget = budget)
      all_results[[idx]] <- result
      idx <- idx + 1
    }
  }

  # Convert to data frame (excluding fit objects)
  results_df <- do.call(rbind, lapply(all_results, function(r) {
    data.frame(
      problem = r$problem,
      seed = r$seed,
      v02_time = r$v02_time,
      v02_evals = r$v02_evals,
      v02_best = r$v02_best,
      v02_feasible_rate = r$v02_feasible_rate,
      v03_time = r$v03_time,
      v03_evals = r$v03_evals,
      v03_best = r$v03_best,
      v03_feasible_rate = r$v03_feasible_rate,
      v03_early_stop = r$v03_early_stop,
      time_speedup = r$time_speedup,
      eval_reduction = r$eval_reduction,
      objective_improvement = r$objective_improvement,
      stringsAsFactors = FALSE
    )
  }))

  # Save results
  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    saveRDS(list(
      results_df = results_df,
      full_results = all_results,
      metadata = list(
        n_init = n_init,
        budget = budget,
        seeds = seeds,
        problems = names(problems),
        timestamp = timestamp
      )
    ), file = sprintf("inst/benchmarks/version_comparison_%s.rds", timestamp))
    cat("\n\nResults saved to: inst/benchmarks/version_comparison_", timestamp, ".rds\n", sep = "")
  }

  # Print summary statistics
  print_benchmark_summary(results_df)

  return(results_df)
}


#' Print Benchmark Summary Statistics
#'
#' @param results_df data frame from run_benchmark_suite()
print_benchmark_summary <- function(results_df) {
  cat("\n\n======================================\n")
  cat("  BENCHMARK SUMMARY STATISTICS\n")
  cat("======================================\n\n")

  # Overall statistics
  cat("Overall Performance (across all problems and seeds):\n")
  cat(sprintf("  Mean time speedup: %.2fx (SD = %.2f)\n",
              mean(results_df$time_speedup, na.rm = TRUE),
              sd(results_df$time_speedup, na.rm = TRUE)))
  cat(sprintf("  Mean evaluation reduction: %.1f%% (SD = %.1f%%)\n",
              mean(results_df$eval_reduction, na.rm = TRUE) * 100,
              sd(results_df$eval_reduction, na.rm = TRUE) * 100))
  cat(sprintf("  Mean objective improvement: %.1f%% (SD = %.1f%%)\n",
              mean(results_df$objective_improvement, na.rm = TRUE) * 100,
              sd(results_df$objective_improvement, na.rm = TRUE) * 100))
  cat(sprintf("  Early stopping rate: %.1f%%\n",
              mean(results_df$v03_early_stop, na.rm = TRUE) * 100))

  # Per-problem statistics
  cat("\n\nPer-Problem Performance:\n")
  for (prob_name in unique(results_df$problem)) {
    prob_data <- results_df[results_df$problem == prob_name, ]
    cat(sprintf("\n%s:\n", prob_name))
    cat(sprintf("  Time speedup: %.2fx (%.2f - %.2f)\n",
                mean(prob_data$time_speedup, na.rm = TRUE),
                min(prob_data$time_speedup, na.rm = TRUE),
                max(prob_data$time_speedup, na.rm = TRUE)))
    cat(sprintf("  Eval reduction: %.1f%% (%.1f%% - %.1f%%)\n",
                mean(prob_data$eval_reduction, na.rm = TRUE) * 100,
                min(prob_data$eval_reduction, na.rm = TRUE) * 100,
                max(prob_data$eval_reduction, na.rm = TRUE) * 100))
    cat(sprintf("  Obj improvement: %.1f%% (SD = %.1f%%)\n",
                mean(prob_data$objective_improvement, na.rm = TRUE) * 100,
                sd(prob_data$objective_improvement, na.rm = TRUE) * 100))
  }

  # Statistical tests
  cat("\n\nStatistical Significance Tests:\n")

  # Paired t-test for time
  time_test <- t.test(results_df$v02_time, results_df$v03_time, paired = TRUE)
  cat(sprintf("  Time difference (paired t-test): p = %.4f %s\n",
              time_test$p.value,
              ifelse(time_test$p.value < 0.05, "***", "")))

  # Paired t-test for evaluations
  eval_test <- t.test(results_df$v02_evals, results_df$v03_evals, paired = TRUE)
  cat(sprintf("  Evaluation difference (paired t-test): p = %.4f %s\n",
              eval_test$p.value,
              ifelse(eval_test$p.value < 0.05, "***", "")))

  # Paired t-test for objective (where both found feasible solution)
  valid_obj <- !is.na(results_df$v02_best) & !is.na(results_df$v03_best)
  if (sum(valid_obj) > 0) {
    obj_test <- t.test(results_df$v02_best[valid_obj],
                       results_df$v03_best[valid_obj],
                       paired = TRUE)
    cat(sprintf("  Objective difference (paired t-test): p = %.4f %s\n",
                obj_test$p.value,
                ifelse(obj_test$p.value < 0.05, "***", "")))
  }

  cat("\n(*** indicates p < 0.05)\n")
  cat("\n======================================\n\n")
}


#' Plot Benchmark Results
#'
#' Creates visualization of v0.2.0 vs v0.3.0 comparison.
#'
#' @param results_df data frame from run_benchmark_suite()
#' @param output_file optional PDF file path for saving
plot_benchmark_results <- function(results_df, output_file = NULL) {
  require(ggplot2)
  require(gridExtra)

  # Plot 1: Time comparison
  p1 <- ggplot(results_df, aes(x = v02_time, y = v03_time, color = problem)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, alpha = 0.3) +
    labs(
      title = "Wall-Clock Time Comparison",
      x = "v0.2.0 Time (seconds)",
      y = "v0.3.0 Time (seconds)",
      subtitle = "Points below diagonal = v0.3.0 faster"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot 2: Evaluation count comparison
  p2 <- ggplot(results_df, aes(x = v02_evals, y = v03_evals, color = problem)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, alpha = 0.3) +
    labs(
      title = "Evaluation Count Comparison",
      x = "v0.2.0 Evaluations",
      y = "v0.3.0 Evaluations",
      subtitle = "Points below diagonal = v0.3.0 more efficient"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot 3: Speedup distribution
  p3 <- ggplot(results_df, aes(x = time_speedup, fill = problem)) +
    geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    geom_vline(aes(xintercept = mean(time_speedup, na.rm = TRUE)),
               linetype = "solid", color = "blue", size = 1) +
    labs(
      title = "Time Speedup Distribution",
      x = "Speedup Factor (v0.2.0 / v0.3.0)",
      y = "Count",
      subtitle = sprintf("Mean = %.2fx faster", mean(results_df$time_speedup, na.rm = TRUE))
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot 4: Objective comparison (where both feasible)
  valid_obj <- !is.na(results_df$v02_best) & !is.na(results_df$v03_best)
  p4 <- ggplot(results_df[valid_obj, ], aes(x = v02_best, y = v03_best, color = problem)) +
    geom_point(alpha = 0.6, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_smooth(method = "lm", se = FALSE, alpha = 0.3) +
    labs(
      title = "Objective Quality Comparison",
      x = "v0.2.0 Best Objective",
      y = "v0.3.0 Best Objective",
      subtitle = "Points below diagonal = v0.3.0 better"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Combine plots
  combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)

  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 12, height = 10)
    cat("Plots saved to:", output_file, "\n")
  }

  return(combined)
}


# ============================================================================
# Example Usage
# ============================================================================

if (FALSE) {
  # Quick test (1 problem, 3 seeds, small budget)
  problems <- get_test_problems()
  results_quick <- run_benchmark_suite(
    problems = problems["toy_2d"],
    seeds = 1:3,
    n_init = 8,
    budget = 30,
    save_results = TRUE
  )
  plot_benchmark_results(results_quick)

  # Full benchmark (all problems, 10 seeds)
  results_full <- run_benchmark_suite(
    problems = get_test_problems(),
    seeds = 1:10,
    n_init = 10,
    budget = 50,
    save_results = TRUE
  )
  plot_benchmark_results(results_full, output_file = "inst/benchmarks/version_comparison.pdf")

  # Load and analyze saved results
  saved <- readRDS("inst/benchmarks/version_comparison_TIMESTAMP.rds")
  print_benchmark_summary(saved$results_df)
  plot_benchmark_results(saved$results_df)
}
