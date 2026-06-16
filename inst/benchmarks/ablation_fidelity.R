# Ablation Study: Multi-Fidelity Selection Strategies
#
# This script performs an ablation study comparing different fidelity selection
# methods implemented in BATON v0.3.0:
#
# 1. adaptive: Cost-aware value-per-cost optimization (NEW, recommended)
# 2. staged: Fixed iteration thresholds (EXISTING)
# 3. threshold: Simple feasibility probability thresholds (EXISTING)
# 4. fixed_low: Always use low fidelity (baseline)
# 5. fixed_high: Always use high fidelity (baseline)
#
# Research Questions:
# - Does adaptive method outperform staged and threshold?
# - How much benefit from multi-fidelity vs fixed fidelity?
# - Is the cost-awareness worth the complexity?
# - How does performance vary across problem characteristics?

library(BATON)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load test problems
source("inst/benchmarks/test_problems.R")

#' Run Single Ablation Experiment
#'
#' Compares all fidelity methods on a single problem and seed.
#'
#' @param problem problem configuration
#' @param seed RNG seed
#' @param n_init initial design size
#' @param q batch size
#' @param budget total evaluation budget
#' @return list with results for all methods
run_single_ablation <- function(problem, seed = 42, n_init = 10, q = 4, budget = 50) {
  cat(sprintf("\n=== Ablation: %s (seed = %d) ===\n", problem$name, seed))

  methods <- list(
    adaptive = list(fidelity_method = "adaptive", levels = c(low = 200, med = 1000, high = 10000)),
    staged = list(fidelity_method = "staged", levels = c(low = 200, med = 1000, high = 10000)),
    threshold = list(fidelity_method = "threshold", levels = c(low = 200, med = 1000, high = 10000)),
    fixed_low = list(fidelity_method = "adaptive", levels = c(low = 200)),
    fixed_high = list(fidelity_method = "adaptive", levels = c(high = 10000))
  )

  results <- list()

  for (method_name in names(methods)) {
    cat(sprintf("  Running %s...\n", method_name))

    method_config <- methods[[method_name]]

    tryCatch({
      start_time <- Sys.time()

      fit <- bo_calibrate(
        sim_fun = problem$sim_fun,
        bounds = problem$bounds,
        objective = problem$objective,
        constraints = problem$constraints,
        n_init = n_init,
        q = q,
        budget = budget,
        fidelity_method = method_config$fidelity_method,
        fidelity_levels = method_config$levels,
        seed = seed,
        progress = FALSE
      )

      end_time <- Sys.time()

      # Compute metrics
      feasible <- fit$history %>% filter(feasible)
      best_obj <- if (nrow(feasible) > 0) min(feasible$objective) else NA

      # Compute total simulation cost (weighted by fidelity)
      costs <- sapply(fit$history$fidelity, function(f) {
        method_config$levels[f]
      })
      total_cost <- sum(costs)

      # Fidelity distribution
      fidelity_dist <- table(fit$history$fidelity)
      fidelity_pct <- prop.table(fidelity_dist) * 100

      results[[method_name]] <- list(
        method = method_name,
        problem = problem$name,
        seed = seed,
        time = as.numeric(difftime(end_time, start_time, units = "secs")),
        evals = nrow(fit$history),
        best_obj = best_obj,
        feasible_rate = mean(fit$history$feasible),
        total_cost = total_cost,
        fidelity_dist = fidelity_dist,
        fidelity_pct = fidelity_pct,
        fit = fit
      )

      cat(sprintf("    %d evals, cost = %d, best = %.2f\n",
                  results[[method_name]]$evals,
                  results[[method_name]]$total_cost,
                  results[[method_name]]$best_obj))
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
      results[[method_name]] <- NULL
    })
  }

  return(results)
}


#' Run Complete Ablation Study
#'
#' Compares all fidelity methods across problems and seeds.
#'
#' @param problems list of problem configurations
#' @param seeds vector of RNG seeds
#' @param n_init initial design size
#' @param q batch size
#' @param budget total evaluation budget
#' @param save_results whether to save results to RDS file
#' @return data frame of results
run_ablation_study <- function(problems = get_test_problems(),
                                seeds = 1:20,
                                n_init = 10,
                                q = 4,
                                budget = 50,
                                save_results = TRUE) {

  cat("\n================================================\n")
  cat("  Multi-Fidelity Ablation Study\n")
  cat("================================================\n")
  cat(sprintf("Problems: %d\n", length(problems)))
  cat(sprintf("Seeds: %d\n", length(seeds)))
  cat(sprintf("Budget: %d evaluations per run\n", budget))
  cat(sprintf("Total runs: %d\n", length(problems) * length(seeds) * 5))
  cat("================================================\n")

  all_results <- list()
  idx <- 1

  for (prob in problems) {
    for (seed in seeds) {
      result <- run_single_ablation(prob, seed = seed, n_init = n_init, q = q, budget = budget)
      all_results[[idx]] <- result
      idx <- idx + 1
    }
  }

  # Convert to data frame
  results_df <- do.call(rbind, lapply(all_results, function(run) {
    do.call(rbind, lapply(run, function(r) {
      if (is.null(r)) return(NULL)
      data.frame(
        method = r$method,
        problem = r$problem,
        seed = r$seed,
        time = r$time,
        evals = r$evals,
        best_obj = r$best_obj,
        feasible_rate = r$feasible_rate,
        total_cost = r$total_cost,
        stringsAsFactors = FALSE
      )
    }))
  }))

  # Save results
  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    saveRDS(list(
      results_df = results_df,
      full_results = all_results,
      metadata = list(
        n_init = n_init,
        q = q,
        budget = budget,
        seeds = seeds,
        problems = names(problems),
        timestamp = timestamp
      )
    ), file = sprintf("inst/benchmarks/ablation_fidelity_%s.rds", timestamp))
    cat("\n\nResults saved to: inst/benchmarks/ablation_fidelity_", timestamp, ".rds\n", sep = "")
  }

  # Print summary
  print_ablation_summary(results_df)

  return(results_df)
}


#' Print Ablation Study Summary
#'
#' @param results_df data frame from run_ablation_study()
print_ablation_summary <- function(results_df) {
  cat("\n\n================================================\n")
  cat("  ABLATION STUDY SUMMARY\n")
  cat("================================================\n\n")

  # Overall performance by method
  cat("Overall Performance by Method:\n")
  method_summary <- results_df %>%
    group_by(method) %>%
    summarize(
      n_runs = n(),
      mean_time = mean(time, na.rm = TRUE),
      mean_evals = mean(evals, na.rm = TRUE),
      mean_best = mean(best_obj, na.rm = TRUE),
      sd_best = sd(best_obj, na.rm = TRUE),
      mean_cost = mean(total_cost, na.rm = TRUE),
      feasible_rate = mean(feasible_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_best)

  print(as.data.frame(method_summary), row.names = FALSE)

  # Per-problem analysis
  cat("\n\nPer-Problem Performance:\n")
  for (prob_name in unique(results_df$problem)) {
    prob_data <- results_df[results_df$problem == prob_name, ]

    cat(sprintf("\n%s:\n", prob_name))

    prob_summary <- prob_data %>%
      group_by(method) %>%
      summarize(
        mean_best = mean(best_obj, na.rm = TRUE),
        sd_best = sd(best_obj, na.rm = TRUE),
        mean_cost = mean(total_cost, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(mean_best)

    print(as.data.frame(prob_summary), row.names = FALSE)
  }

  # Statistical comparisons
  cat("\n\nStatistical Comparisons (Adaptive vs Others):\n")
  for (method_name in c("staged", "threshold", "fixed_low", "fixed_high")) {
    # Paired comparison (same problem + seed)
    adaptive_data <- results_df %>%
      filter(method == "adaptive") %>%
      arrange(problem, seed)

    method_data <- results_df %>%
      filter(method == method_name) %>%
      arrange(problem, seed)

    if (nrow(adaptive_data) == nrow(method_data)) {
      # Objective comparison
      valid <- !is.na(adaptive_data$best_obj) & !is.na(method_data$best_obj)
      if (sum(valid) > 0) {
        obj_test <- t.test(adaptive_data$best_obj[valid],
                           method_data$best_obj[valid],
                           paired = TRUE)

        # Cost comparison
        cost_test <- t.test(adaptive_data$total_cost[valid],
                            method_data$total_cost[valid],
                            paired = TRUE)

        cat(sprintf("\nAdaptive vs %s:\n", method_name))
        cat(sprintf("  Objective: mean diff = %.3f, p = %.4f %s\n",
                    mean(adaptive_data$best_obj[valid] - method_data$best_obj[valid]),
                    obj_test$p.value,
                    ifelse(obj_test$p.value < 0.05, "***", "")))
        cat(sprintf("  Cost: mean diff = %.0f, p = %.4f %s\n",
                    mean(adaptive_data$total_cost[valid] - method_data$total_cost[valid]),
                    cost_test$p.value,
                    ifelse(cost_test$p.value < 0.05, "***", "")))
      }
    }
  }

  cat("\n(*** indicates p < 0.05)\n")
  cat("\n================================================\n\n")
}


#' Plot Ablation Study Results
#'
#' Creates comprehensive visualization of ablation study.
#'
#' @param results_df data frame from run_ablation_study()
#' @param output_file optional PDF file path for saving
plot_ablation_results <- function(results_df, output_file = NULL) {
  require(ggplot2)
  require(gridExtra)

  # Plot 1: Objective quality by method
  p1 <- ggplot(results_df, aes(x = method, y = best_obj, fill = method)) +
    geom_boxplot() +
    facet_wrap(~problem, scales = "free_y") +
    labs(
      title = "Objective Quality by Method",
      x = "Fidelity Method",
      y = "Best Objective Value",
      subtitle = "Lower is better"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Plot 2: Total cost by method
  p2 <- ggplot(results_df, aes(x = method, y = total_cost, fill = method)) +
    geom_boxplot() +
    facet_wrap(~problem, scales = "free_y") +
    labs(
      title = "Total Simulation Cost by Method",
      x = "Fidelity Method",
      y = "Total Cost (sum of replications)",
      subtitle = "Lower is more efficient"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Plot 3: Cost-efficiency trade-off
  method_summary <- results_df %>%
    group_by(method, problem) %>%
    summarize(
      mean_cost = mean(total_cost, na.rm = TRUE),
      mean_obj = mean(best_obj, na.rm = TRUE),
      .groups = "drop"
    )

  p3 <- ggplot(method_summary, aes(x = mean_cost, y = mean_obj, color = method, shape = problem)) +
    geom_point(size = 4) +
    geom_text(aes(label = method), vjust = -0.5, size = 3) +
    labs(
      title = "Cost-Efficiency Trade-off",
      x = "Mean Total Cost",
      y = "Mean Best Objective",
      subtitle = "Lower-left is better (low cost, low objective)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Plot 4: Feasibility rate by method
  p4 <- ggplot(results_df, aes(x = method, y = feasible_rate, fill = method)) +
    geom_boxplot() +
    facet_wrap(~problem) +
    labs(
      title = "Feasibility Rate by Method",
      x = "Fidelity Method",
      y = "Proportion of Feasible Evaluations",
      subtitle = "Higher is better"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  # Combine plots
  combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)

  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, combined, width = 14, height = 12)
    cat("Plots saved to:", output_file, "\n")
  }

  return(combined)
}


#' Analyze Cost-Efficiency
#'
#' Computes cost-normalized objective for ranking methods.
#'
#' @param results_df data frame from run_ablation_study()
#' @return data frame with cost-efficiency metrics
analyze_cost_efficiency <- function(results_df) {
  # Normalize cost and objective within each problem
  efficiency <- results_df %>%
    group_by(problem) %>%
    mutate(
      cost_norm = (total_cost - min(total_cost)) / (max(total_cost) - min(total_cost) + 1e-10),
      obj_norm = (best_obj - min(best_obj, na.rm = TRUE)) / (max(best_obj, na.rm = TRUE) - min(best_obj, na.rm = TRUE) + 1e-10)
    ) %>%
    ungroup() %>%
    mutate(
      # Efficiency score: lower is better (balance cost and quality)
      efficiency_score = 0.5 * cost_norm + 0.5 * obj_norm
    )

  # Summarize by method
  method_efficiency <- efficiency %>%
    group_by(method) %>%
    summarize(
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      sd_efficiency = sd(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(mean_efficiency)

  cat("\nCost-Efficiency Ranking (lower is better):\n")
  print(as.data.frame(method_efficiency), row.names = FALSE)

  return(method_efficiency)
}


# ============================================================================
# Example Usage
# ============================================================================

if (FALSE) {
  # Quick test (1 problem, 3 seeds)
  problems <- get_test_problems()
  results_quick <- run_ablation_study(
    problems = problems["toy_2d"],
    seeds = 1:3,
    n_init = 8,
    q = 4,
    budget = 30,
    save_results = TRUE
  )
  plot_ablation_results(results_quick)
  analyze_cost_efficiency(results_quick)

  # Full ablation study (all problems, 20 seeds for statistical power)
  results_full <- run_ablation_study(
    problems = get_test_problems(),
    seeds = 1:20,
    n_init = 10,
    q = 4,
    budget = 50,
    save_results = TRUE
  )
  plot_ablation_results(results_full, output_file = "inst/benchmarks/ablation_fidelity.pdf")
  analyze_cost_efficiency(results_full)

  # Load and analyze saved results
  saved <- readRDS("inst/benchmarks/ablation_fidelity_TIMESTAMP.rds")
  print_ablation_summary(saved$results_df)
  plot_ablation_results(saved$results_df)
  analyze_cost_efficiency(saved$results_df)
}
