# BATON Benchmarking Suite

This directory contains scripts and utilities for benchmarking the performance improvements in BATON v0.3.0.

## Overview

Version 0.3.0 introduces substantial performance improvements across three implementation phases:
- **Phase 1**: Acquisition improvements & batch diversity (10-20% fewer evaluations)
- **Phase 2**: Adaptive fidelity selection (15-25% better budget use)
- **Phase 3**: Performance optimizations (30-60% faster wall-clock time)

**Expected combined improvement**: 50-70% overall efficiency gain

This benchmarking suite validates these claims through:
1. **Version comparison**: v0.2.0 vs v0.3.0
2. **Ablation study**: Comparison of fidelity selection methods
3. **Standardized test problems**: Reproducible benchmarks

## Files

| File | Purpose |
|------|---------|
| `test_problems.R` | Standard test problem definitions |
| `compare_versions.R` | Version comparison benchmark (v0.2.0 vs v0.3.0) |
| `ablation_fidelity.R` | Ablation study for fidelity selection methods |
| `README.md` | This file |

## Test Problems

Three standard test problems are provided:

### 1. toy_2d
- **Dimension**: 2
- **Characteristics**: Simple quadratic objective, soft constraints
- **Purpose**: Quick validation, algorithm development
- **Typical runtime**: ~1-2 minutes per run

### 2. high_dim_5d
- **Dimension**: 5
- **Characteristics**: Rosenbrock-like function, coupled constraints
- **Purpose**: Test scalability to higher dimensions
- **Typical runtime**: ~3-5 minutes per run

### 3. tight_constraints
- **Dimension**: 2
- **Characteristics**: Small feasible region, difficult constraints
- **Purpose**: Test constraint handling in infeasible regions
- **Typical runtime**: ~1-2 minutes per run

## Quick Start

### Version Comparison

Compare v0.2.0 behavior vs v0.3.0 improvements:

```r
source("inst/benchmarks/test_problems.R")
source("inst/benchmarks/compare_versions.R")

# Quick test (1 problem, 3 seeds)
problems <- get_test_problems()
results <- run_benchmark_suite(
  problems = problems["toy_2d"],
  seeds = 1:3,
  budget = 30
)

# View results
print_benchmark_summary(results)
plot_benchmark_results(results)
```

### Ablation Study

Compare fidelity selection methods:

```r
source("inst/benchmarks/test_problems.R")
source("inst/benchmarks/ablation_fidelity.R")

# Quick test (1 problem, 3 seeds)
problems <- get_test_problems()
results <- run_ablation_study(
  problems = problems["toy_2d"],
  seeds = 1:3,
  budget = 30
)

# View results
print_ablation_summary(results)
plot_ablation_results(results)
analyze_cost_efficiency(results)
```

## Full Benchmark Suite

For publication-quality results:

```r
# Version comparison (all problems, 10 seeds)
results_version <- run_benchmark_suite(
  problems = get_test_problems(),
  seeds = 1:10,
  n_init = 10,
  budget = 50,
  save_results = TRUE
)
plot_benchmark_results(results_version, output_file = "inst/benchmarks/version_comparison.pdf")

# Ablation study (all problems, 20 seeds for statistical power)
results_ablation <- run_ablation_study(
  problems = get_test_problems(),
  seeds = 1:20,
  n_init = 10,
  q = 2,
  budget = 50,
  save_results = TRUE
)
plot_ablation_results(results_ablation, output_file = "inst/benchmarks/ablation_fidelity.pdf")
analyze_cost_efficiency(results_ablation)
```

**Runtime**: ~6-8 hours for complete suite (can be parallelized by problem/seed)

## Expected Results

### Version Comparison (v0.2.0 vs v0.3.0)

**Conservative estimates**:
- Time speedup: 1.5-2x faster
- Evaluation reduction: 10-20%
- Early stopping rate: 20-40%
- Objective quality: Similar or better

**Optimistic estimates** (favorable problems):
- Time speedup: 2-3x faster
- Evaluation reduction: 20-30%
- Early stopping rate: 30-50%
- Objective quality: 5-10% better

### Ablation Study

**Expected ranking** (by cost-efficiency):
1. **adaptive**: Best cost-efficiency (recommended default)
2. **staged**: Good compromise, simple
3. **threshold**: Basic, less efficient
4. **fixed_high**: Best quality, high cost
5. **fixed_low**: Lowest cost, poor quality

**Statistical significance**:
- Adaptive should significantly outperform threshold and fixed methods (p < 0.05)
- Adaptive vs staged may be marginal on simple problems
- Differences should be clearer on high-dimensional and tight-constraint problems

## Analyzing Saved Results

Results are automatically saved with timestamps:

```r
# Load saved results
version_results <- readRDS("inst/benchmarks/version_comparison_YYYYMMDD_HHMMSS.rds")
ablation_results <- readRDS("inst/benchmarks/ablation_fidelity_YYYYMMDD_HHMMSS.rds")

# Analyze
print_benchmark_summary(version_results$results_df)
plot_benchmark_results(version_results$results_df)

print_ablation_summary(ablation_results$results_df)
plot_ablation_results(ablation_results$results_df)
analyze_cost_efficiency(ablation_results$results_df)
```

## Creating Custom Test Problems

To add a custom problem:

```r
# 1. Define simulator function
my_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  # ... implementation ...
  result <- welford_mean_var(...)
  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- n_rep
  return(metrics)
}

# 2. Define problem configuration
my_problem <- function() {
  list(
    name = "my_problem",
    sim_fun = my_sim,
    bounds = list(x1 = c(0, 1), x2 = c(0, 1)),
    objective = "metric_name",
    constraints = list(power = c("ge", 0.8)),
    optimum = list(x1 = 0.5, x2 = 0.5),  # Optional
    description = "My custom test problem"
  )
}

# 3. Run benchmarks
problems <- list(my_problem = my_problem())
results <- run_benchmark_suite(problems = problems, seeds = 1:5)
```

## Benchmark Metrics

### Primary Metrics

- **Time (seconds)**: Wall-clock time for full BO run
- **Evaluations**: Number of simulator calls
- **Best objective**: Final objective value (feasible region)
- **Total cost**: Sum of simulation replications (weighted by fidelity)

### Secondary Metrics

- **Feasible rate**: Proportion of evaluations in feasible region
- **Early stopping rate**: Proportion of runs that stopped early
- **Fidelity distribution**: Usage of low/med/high fidelity
- **Speedup**: Time ratio (baseline / improved)
- **Efficiency score**: Cost-normalized objective quality

## Statistical Analysis

All comparisons use paired t-tests (same problem + seed):
- **Objective difference**: Is v0.3.0 better?
- **Time difference**: Is v0.3.0 faster?
- **Cost difference**: Is adaptive more efficient?

Significance threshold: p < 0.05

## Troubleshooting

### Benchmarks take too long

- Reduce `budget` (default 50 → 30)
- Reduce number of seeds (default 10 → 3)
- Run single problem at a time
- Use parallel execution (future/furrr)

```r
# Parallel execution example
library(furrr)
plan(multisession, workers = 4)

results <- future_map_dfr(1:10, function(seed) {
  run_single_benchmark(problem, seed = seed)
})
```

### Out of memory

Test problems use Welford's algorithm (minimal memory), but storing full fit objects can be large.

```r
# Don't save full fit objects
results <- run_benchmark_suite(..., save_results = FALSE)

# Or manually cleanup
results$full_results <- NULL  # Remove fit objects
saveRDS(results$results_df, "results_light.rds")
```

### Inconsistent results

- Ensure same `seed` for reproducibility
- Check R version and package versions
- Verify simulator is deterministic given seed
- Use more seeds (20-30) for stable statistics

## References

For implementation details and expected improvements, see:
- `PHASE1_IMPLEMENTATION_SUMMARY.md`: Batch diversity & acquisition
- `PHASE2_IMPLEMENTATION_SUMMARY.md`: Adaptive fidelity
- `PHASE3_IMPLEMENTATION_SUMMARY.md`: Performance optimizations
- `ALL_PHASES_COMPLETE.md`: Complete overview

## Citation

If you use these benchmarks in your research:

```bibtex
@software{BATON_benchmarks,
  title = {BATON Benchmarking Suite},
  author = {Rashid, Naim},
  year = {2025},
  url = {https://github.com/naimurashid/BATON},
  version = {0.3.0},
  note = {Test problems and validation scripts for BATON performance improvements}
}
```

## Contributing

To contribute new test problems or benchmark scripts:

1. Follow the existing structure in `test_problems.R`
2. Ensure simulators use `welford_mean_var()` for variance estimation
3. Include known optimum if possible
4. Document expected characteristics
5. Test on small seed set before full run

## Contact

- **Issues**: [GitHub Issues](https://github.com/naimurashid/BATON/issues)
- **Author**: Naim Rashid (naim_rashid@unc.edu)
