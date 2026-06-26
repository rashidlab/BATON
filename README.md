# BATON

**Bayesian Optimization for Calibration of Adaptive Clinical Trials**

Version 0.4.0 | R >= 4.2 | UNC not-for-profit license | [GitHub](https://github.com/naimurashid/BATON) | [Issues](https://github.com/naimurashid/BATON/issues)

**v0.4.0 (May 2026):** hardening release. New `warmstart_from` parameter
(Stage 0 cross-philosophy seeds) and `multi_seed_verify` gate (Stage 4
verification at multiple seeds at high fidelity) on `bo_calibrate()`,
plus a new `bo_calibrate_philosophies()` batch wrapper. Backward-
compatible defaults; opt-in for new behavior. See `NEWS.md` for the
full rationale, migration notes, and the manuscript campaign that
motivated these additions.

Designing a clinical trial requires choosing design parameters -- for
example, the maximum sample size, efficacy and futility decision thresholds,
or the number of interim analyses -- so that the trial achieves adequate
power while controlling type I error, and simultaneously minimizes a
quantity of interest such as the expected sample size, the maximum sample
size, or a weighted combination of both. These performance metrics are
collectively called **operating characteristics**. For adaptive designs,
operating characteristics are often impossible to compute analytically and
must be estimated by simulating thousands of trial replicates.

BATON automates this calibration process. Given a simulator that maps design
parameters to operating characteristics, BATON uses constrained Bayesian
optimization to find configurations that minimize a user-specified objective
(e.g., expected sample size) subject to constraints (e.g., power >= 0.80,
type I error <= 0.10), without requiring closed-form gradients or analytic
solutions.

BATON works with any trial simulator, whether the underlying design is
Bayesian or frequentist. It is most useful when (1) analytic solutions for
operating characteristics do not exist, (2) the design space is complex
enough that manual calibration is difficult, and (3) the number of design
parameters is large enough that grid search becomes infeasible.

BATON uses heteroskedastic Gaussian process surrogates, the expected
constrained improvement acquisition function, multi-fidelity simulation
budgeting, and Welford-based variance estimation.

## Installation

### System Dependencies

`DiceKriging` and `hetGP` require C and Fortran compilers.

| Platform | Command |
|----------|---------|
| Ubuntu/Debian | `sudo apt-get install build-essential gfortran` |
| Fedora/RHEL | `sudo dnf install gcc gcc-c++ gcc-gfortran` |
| macOS | `xcode-select --install` then install `gfortran` from <https://mac.r-project.org/tools/> or `brew install gcc` |
| Windows | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version |

### R Package

```r
install.packages("remotes")
remotes::install_github("naimurashid/BATON")

# With vignettes (requires knitr, rmarkdown):
remotes::install_github("naimurashid/BATON", build_vignettes = TRUE)
```

BATON requires R >= 4.2 (see `DESCRIPTION`) and pairs with the companion
simulator package
[evolveTrial](https://github.com/naimurashid/evolveTrial), which requires
the same R version.

### Troubleshooting

| Error | Fix |
|-------|-----|
| `gfortran: command not found` (during `DiceKriging` install) | Install Fortran per the System Dependencies table above (macOS: `brew install gcc`; Linux: install `gfortran`) |
| `ld: library not found for -lgfortran` (macOS) | `brew install gcc` and ensure the gfortran path is on `PATH` (see <https://mac.r-project.org/tools/>) |
| `hetGP` fails to link BLAS/LAPACK (Linux) | `sudo apt-get install libblas-dev liblapack-dev` (Ubuntu/Debian) or `sudo dnf install blas-devel lapack-devel` (Fedora/RHEL) |
| `make: not found` (Windows) | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version; verify with `Sys.which("make")` |
| `plan(multicore)` runs sequentially on Windows | Use `future::plan(multisession, workers = N)` — Windows has no fork; see "Parallelism" below |
| Installation times out on CI or slow links | `remotes::install_github("naimurashid/BATON", build_vignettes = FALSE)` |
| `Error: package or namespace load failed for 'BATON' ... evolveTrial` | Install the companion [evolveTrial](https://github.com/naimurashid/evolveTrial) package first |

## Quick Start

Self-contained example using a synthetic simulator. Copy-paste and run after
installing BATON.

```r
library(BATON)

# Step 1: Define a simulator that evaluates under BOTH hypotheses
# The evaluator runs n_rep simulated trials under null (p = 0.2) and
# alternative (p = 0.4), returning power, type I error, and expected N.
# The fidelity argument controls R (number of Monte Carlo replications);
# BATON maps categorical levels to replication counts internally.
my_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  n_rep <- switch(fidelity, low = 500, med = 2000, high = 5000)
  thr   <- theta$threshold
  nmax  <- round(theta$nmax)
  if (!is.null(seed)) set.seed(seed)

  # Simulate under null hypothesis (H0: p = 0.2)
  null_results <- welford_mean_var(
    sample_fn = function(i, theta) {
      x <- rbinom(1, size = nmax, prob = 0.2)
      c(reject = as.numeric(x / nmax > theta$threshold), N = nmax)
    }, n_samples = n_rep, theta = theta
  )

  # Simulate under alternative hypothesis (H1: p = 0.4)
  alt_results <- welford_mean_var(
    sample_fn = function(i, theta) {
      x <- rbinom(1, size = nmax, prob = 0.4)
      c(reject = as.numeric(x / nmax > theta$threshold), N = nmax)
    }, n_samples = n_rep, theta = theta
  )

  metrics <- c(
    power = unname(alt_results$mean["reject"]),
    type1 = unname(null_results$mean["reject"]),
    EN    = unname(null_results$mean["N"])
  )
  attr(metrics, "variance") <- c(
    power = unname(alt_results$variance["reject"]),
    type1 = unname(null_results$variance["reject"]),
    EN    = unname(null_results$variance["N"])
  )
  attr(metrics, "n_rep") <- n_rep
  metrics
}

# Step 2: Bounds and constraints
bounds <- list(threshold = c(0.20, 0.45), nmax = c(30, 80))
constraints <- list(power = c("ge", 0.80), type1 = c("le", 0.10))

# Step 3: Run optimization
fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, q = 2, budget = 60, seed = 42
)

# Step 4: Extract results
fit$best_theta       # best feasible design parameters
fit$best_objective   # minimum EN among feasible designs

# Get metrics for the best design (no fit$best_metrics field exists)
feasible <- fit$history[fit$history$feasible, ]
best_row <- feasible[which.min(feasible$objective), ]
best_row[, c("power", "type1", "EN")]
```

## Choosing `n_init` and `budget`

Two key parameters control BATON's computational effort: `n_init` (number of
initial space-filling evaluations before surrogate-guided search begins) and
`budget` (total number of evaluations including `n_init`). The right settings
depend on the dimensionality of the design space `d` (the number of
parameters in theta), the cost per evaluation, and the complexity of the
feasibility region.

### Rules of Thumb

| Parameter | Guideline |
|-----------|-----------|
| `n_init` | 10 x `d` as a starting point. For `d = 2`, use 20; for `d = 4`-`6`, use 40-60. Fewer initial points risk poor surrogate fits; more waste budget on un-guided exploration. |
| `budget` | 20-50 x `d` total evaluations for well-behaved problems. For `d = 2`, 60-100 is typically sufficient; for `d = 6`+, budget 150-300+. |
| `q` (batch size) | 2-4 for parallel evaluation; 1 for sequential. Larger batches trade statistical efficiency for wall-clock speedup. |

### Adjusting for Problem Difficulty

- **Narrow feasibility regions** (tight constraints on power AND type I
  error) require more evaluations to locate the feasible boundary. Increase
  `budget` by 50-100% relative to the rules of thumb above.
- **Expensive evaluators** (e.g., survival trials with R = 5,000) benefit
  from multi-fidelity scheduling (`fidelity_method = "adaptive"`), which
  stretches the budget by using cheap low-fidelity evaluations for
  exploration and reserving high-fidelity evaluations for promising regions.
- **Multi-stage warm-starting** is recommended for `d >= 4`: run a broad
  Stage 1 with moderate budget, then narrow bounds around the best region
  and warm-start Stage 2 with `initial_history = fit1$history`.

### Example: Scaling with Dimension

```r
# d = 2 (e.g., threshold + nmax)
fit <- bo_calibrate(..., n_init = 20, budget = 60)

# d = 4 (e.g., eff_threshold + fut_threshold + nmax + interim_fraction)
fit <- bo_calibrate(..., n_init = 40, budget = 150)

# d = 6+ (e.g., hybrid seamless design with SA + BA parameters)
# Use multi-stage workflow
fit1 <- bo_calibrate(..., n_init = 60, budget = 300)  # Stage 1: broad
narrow <- refine_bounds(fit1, shrink_factor = 0.5)
fit2 <- bo_calibrate(..., bounds = narrow,
  initial_history = fit1$history, budget = 200)        # Stage 2: refine
```

## Writing Your Own Evaluator

BATON works with **any** trial simulator -- whether the trial uses Bayesian
decision rules, frequentist group sequential boundaries, or hybrid designs.
You supply an evaluator function that simulates trials and returns their
operating characteristics (power, type I error, expected sample size); BATON
handles the optimization. This section walks through building one from
scratch.

### Evaluator Requirements

Your evaluator must satisfy three requirements:

1. **Signature:** `function(theta, fidelity = c("low", "med", "high"), ...)`
   - `theta` is a named list of design parameters (e.g., `theta$alpha`, `theta$nmax`)
   - `fidelity` controls the number of Monte Carlo replications
2. **Return value:** A named numeric vector of operating characteristics
   (e.g., `c(power = 0.83, type1 = 0.07, EN = 52)`)
3. **Variance attributes (recommended):** Attach `attr(result, "variance")` and
   `attr(result, "n_rep")` to enable the heteroskedastic GP surrogate, which
   improves efficiency by 30-50%

The names in the return vector must match the names used in `objective` and
`constraints`. For example, if you specify `constraints = list(power = c("ge", 0.80))`,
your evaluator must return a value named `"power"`.

### Complete Example: Group Sequential Design

This evaluator simulates a two-stage group sequential trial with a continuous
endpoint. It is entirely self-contained -- no external packages required beyond
base R.

```r
# Evaluator for a two-stage group sequential trial
# Design parameters: alpha_spend (Stage 1 alpha), nmax (max per arm)
# Fixed: effect size = 0.3, interim at 50% enrollment
gs_evaluator <- function(theta, fidelity = c("low", "med", "high"),
                         seed = NULL, ...) {
  fidelity <- match.arg(fidelity)
  n_rep <- switch(fidelity, low = 500, med = 2000, high = 10000)
  if (!is.null(seed)) set.seed(seed)

  alpha_spend <- theta$alpha_spend  # alpha spent at interim
  nmax        <- round(theta$nmax)  # max patients per arm
  n_interim   <- round(nmax / 2)    # interim at 50%
  delta       <- 0.3                # true effect size (alternative)

  # Simulate one trial under a given hypothesis
  sim_one <- function(i, theta, null = FALSE) {
    effect <- if (null) 0 else delta

    # Stage 1
    y1_trt <- rnorm(n_interim, mean = effect, sd = 1)
    y1_ctl <- rnorm(n_interim, mean = 0, sd = 1)
    z1 <- (mean(y1_trt) - mean(y1_ctl)) / sqrt(2 / n_interim)

    # Interim decision
    crit1 <- qnorm(1 - alpha_spend)
    if (z1 > crit1) return(c(reject = 1, N = 2 * n_interim))

    # Stage 2: enroll remaining patients
    y2_trt <- rnorm(nmax - n_interim, mean = effect, sd = 1)
    y2_ctl <- rnorm(nmax - n_interim, mean = 0, sd = 1)
    y_trt <- c(y1_trt, y2_trt)
    y_ctl <- c(y1_ctl, y2_ctl)
    z_final <- (mean(y_trt) - mean(y_ctl)) / sqrt(2 / nmax)

    # Final decision (remaining alpha)
    alpha_remain <- 0.10 - alpha_spend
    crit_final <- qnorm(1 - alpha_remain)
    c(reject = as.numeric(z_final > crit_final), N = 2 * nmax)
  }

  # Run under alternative (power) and null (type I error)
  alt_results <- welford_mean_var(
    sample_fn = function(i, theta) sim_one(i, theta, null = FALSE),
    n_samples = n_rep, theta = theta
  )
  null_results <- welford_mean_var(
    sample_fn = function(i, theta) sim_one(i, theta, null = TRUE),
    n_samples = n_rep, theta = theta
  )

  # Combine into return vector
  metrics <- c(
    power    = unname(alt_results$mean["reject"]),
    type1    = unname(null_results$mean["reject"]),
    EN_null  = unname(null_results$mean["N"]),
    EN_alt   = unname(alt_results$mean["N"])
  )
  attr(metrics, "variance") <- c(
    power    = unname(alt_results$variance["reject"]),
    type1    = unname(null_results$variance["reject"]),
    EN_null  = unname(null_results$variance["N"]),
    EN_alt   = unname(alt_results$variance["N"])
  )
  attr(metrics, "n_rep") <- n_rep
  metrics
}

# Calibrate
library(BATON)
fit <- bo_calibrate(
  sim_fun     = gs_evaluator,
  bounds      = list(alpha_spend = c(0.001, 0.08), nmax = c(50, 200)),
  objective   = "EN_null",
  constraints = list(power = c("ge", 0.80), type1 = c("le", 0.10)),
  n_init = 20, budget = 80, seed = 2025
)
fit$best_theta
```

### Key Patterns

**Null vs alternative:** Run simulations under both hypotheses inside your
evaluator (as above). Return `power` from the alternative and `type1` from
the null in a single vector. BATON constrains both simultaneously.

**Fidelity:** Use the `fidelity` argument to set the replication count. Low
fidelity (fewer reps) is used for exploration; high fidelity for verification.
BATON's multi-fidelity engine manages the schedule automatically.

**Seeds:** Accept `seed` via `...` for reproducibility. BATON passes different
seeds across iterations.

**Welford variance:** Use `welford_mean_var()` rather than storing all
replicates. This computes both the mean and variance of each metric in a
single pass, with O(m) memory regardless of the number of replications.

**Existing simulators:** If you already have a simulator that returns a
data.frame or list, write a thin wrapper that extracts the metrics into a
named vector and attaches variance attributes.

See `system.file("examples/simulator_with_variance.R", package = "BATON")`
for additional examples including parallel execution.

### Variance Attributes and Surrogate Selection

BATON uses heteroskedastic Gaussian processes (hetGP) when per-evaluation
variance estimates are available, and falls back to homoskedastic GPs
(DiceKriging) when they are not. The surrogate choice is made **per metric**:

| Scenario | Surrogate Used | How to Enable |
|----------|---------------|---------------|
| Evaluator attaches `attr(result, "variance")` for metric | **hetGP** (recommended) | Use `welford_mean_var()` in your evaluator |
| No variance attribute, metric in [0,1] (e.g., power, type I) | **hetGP** via binomial approximation | Automatic: BATON estimates `p(1-p)/n_rep` |
| No variance attribute, metric outside [0,1] (e.g., EN, ET) | **DiceKriging** (homoskedastic) | Provide variance attribute to upgrade |
| Analytical evaluator (no MC noise) | **DiceKriging** with small nugget | Appropriate: no noise to model |

**Why this matters:** Monte Carlo variance of expected sample size (EN) is
input-dependent: designs with aggressive early stopping produce bimodal
N distributions (high variance), while designs that rarely stop early have
N concentrated near N_max (low variance). A homoskedastic GP assumes
constant noise, oversmoothing high-variance regions and undersmoothing
low-variance regions. In our benchmarks, providing EN variance improved the
best feasible design by approximately 15%.

**What to do if you cannot compute per-replication variance:**

If your simulator only returns aggregated means (e.g., from an external
program that reports `E[N]` but not individual trial sample sizes), BATON
still works. Constraint surrogates (power, type I error) automatically get
binomial variance estimates and use hetGP. Only the objective surrogate (EN)
falls back to DiceKriging. This is suboptimal but not incorrect. Options:

1. **Modify the simulator** to store per-replication values and compute
   `var(N_per_rep) / n_rep`. This is the recommended approach.
2. **Use `welford_mean_var()`** to wrap individual trial simulations if you
   have access to the per-trial simulation function.
3. **Accept the fallback** for rapid prototyping. Increase budget by 20-50%
   to compensate for the less efficient surrogate.

```r
# Example: wrapping an existing simulator that only returns means
my_wrapper <- function(theta, fidelity = "high", seed = NULL, ...) {
  n_rep <- switch(fidelity, low = 500, med = 2000, high = 5000)

  # Option A: If you can call the simulator per-trial
  result <- welford_mean_var(
    sample_fn = function(i, theta) run_one_trial(theta),
    n_samples = n_rep, theta = theta
  )
  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- n_rep
  return(metrics)

  # Option B: If you only have the aggregated output
  # BATON will use hetGP for proportions, DiceKriging for EN
  agg <- run_simulator_batch(theta, n_rep = n_rep)
  c(power = agg$power, type1 = agg$type1, EN = agg$EN)
}
```

### Parallel Evaluation Within Your Evaluator

Each call to `bo_calibrate()` evaluates your simulator many times
sequentially (or in small batches via `q`). The primary opportunity for
parallelism is **inside your evaluator**: splitting Monte Carlo replications
across CPU cores so each evaluation completes faster.

BATON provides `pool_welford_results()` to combine chunk-level Welford
statistics from independent workers. The pattern works with any parallel
backend. By default, the number of cores is autodetected, but you can
override it manually.

**Cross-platform parallel evaluator:**

```r
my_sim <- function(theta, fidelity = "high", seed = NULL, ncores = NULL, ...) {
  n_rep <- switch(fidelity, low = 500, med = 2000, high = 5000)
  if (!is.null(seed)) set.seed(seed)

  # Autodetect cores, or let user override
  if (is.null(ncores)) {
    ncores <- max(1, parallel::detectCores(logical = FALSE) - 1)
  }
  chunk_size <- ceiling(n_rep / ncores)

  # Define the single-trial function
  sim_one_null <- function(i, theta) {
    x <- rbinom(1, size = round(theta$nmax), prob = 0.2)
    c(reject = as.numeric(x / round(theta$nmax) > theta$threshold),
      N = round(theta$nmax))
  }
  sim_one_alt <- function(i, theta) {
    x <- rbinom(1, size = round(theta$nmax), prob = 0.4)
    c(reject = as.numeric(x / round(theta$nmax) > theta$threshold),
      N = round(theta$nmax))
  }

  # Run chunks in parallel (works on macOS/Linux; Windows uses socket cluster)
  run_chunk <- function(hypothesis_fn) {
    if (.Platform$OS.type == "unix") {
      # Fork-based parallelism (macOS, Linux)
      chunks <- parallel::mclapply(seq_len(ncores), function(k) {
        welford_mean_var(
          sample_fn = hypothesis_fn,
          n_samples = chunk_size,
          theta = theta
        )
      }, mc.cores = ncores)
    } else {
      # Socket-based parallelism (Windows)
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterExport(cl, c("theta", "chunk_size", "hypothesis_fn"),
                              envir = environment())
      chunks <- parallel::parLapply(cl, seq_len(ncores), function(k) {
        BATON::welford_mean_var(
          sample_fn = hypothesis_fn,
          n_samples = chunk_size,
          theta = theta
        )
      })
    }
    pool_welford_results(chunks)
  }

  null_results <- run_chunk(sim_one_null)
  alt_results  <- run_chunk(sim_one_alt)

  metrics <- c(
    power = unname(alt_results$mean["reject"]),
    type1 = unname(null_results$mean["reject"]),
    EN    = unname(null_results$mean["N"])
  )
  attr(metrics, "variance") <- c(
    power = unname(alt_results$variance["reject"]),
    type1 = unname(null_results$variance["reject"]),
    EN    = unname(null_results$variance["N"])
  )
  attr(metrics, "n_rep") <- null_results$n
  metrics
}

# Define bounds and constraints (must match the metrics returned by my_sim)
bounds      <- list(threshold = c(0.20, 0.45), nmax = c(30, 80))
constraints <- list(power = c("ge", 0.80), type1 = c("le", 0.10))

# Use all available cores (autodetect)
fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, budget = 60, seed = 42
)

# Or specify cores explicitly
fit <- bo_calibrate(
  sim_fun = function(theta, ...) my_sim(theta, ..., ncores = 4),
  bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, budget = 60, seed = 42
)
```

**Key points:**
- `parallel::detectCores(logical = FALSE)` returns physical cores; subtract
  1 to leave a core free for the BO loop itself.
- On **macOS/Linux**, `mclapply` uses fork-based parallelism (no data
  copying overhead).
- On **Windows**, forking is not available; use `makeCluster`/`parLapply`
  with socket-based parallelism instead.
- `pool_welford_results()` uses Chan's parallel variance algorithm to
  correctly combine mean and variance estimates across chunks, maintaining
  the Welford guarantee of numerically stable one-pass computation.
- Each chunk gets `ceiling(n_rep / ncores)` replications, so the total
  may slightly exceed `n_rep`. This is harmless.

## Design Philosophies via Weighted Loss

The central contribution of the JASA paper is a weighted loss framework:

```
L(w; theta) = w_N * N_max + (1 - w_N) * [w_0 * E_0[N] + w_1 * E_1[N]]
```

where `E_0[N]` and `E_1[N]` are expected sample sizes under null and
alternative, and `N_max` is the maximum sample size. Varying the weights
produces different design philosophies, all subject to the same constraints:

| Philosophy | w_N | w_0 | w_1 | Target |
|------------|:---:|:---:|:---:|--------|
| H0-Optimal | 0 | 1 | 0 | Minimize E[N] under null |
| Minimax | 1 | - | - | Minimize maximum sample size |
| Balanced (Fleming) | 0 | 0.5 | 0.5 | Balance E[N] under both hypotheses |
| H1-Optimal | 0 | 0 | 1 | Minimize E[N] under alternative |
| Admissible | varies | varies | varies | Pareto-optimal trade-offs |

BATON has no built-in weighted-loss constructor: `objective` (and every name in
`constraints`) must be the **name of a metric your simulator returns**. There is
no automatic "weighted sum" objective; if you pass an `objective` the simulator
does not return, `bo_calibrate()` errors with `Objective '...' not returned by
simulator`. To use the weighted loss, have the **simulator itself compute and
return** the combined metric, then point `objective` at it:

```r
# Choose the weights for the desired philosophy (see table above).
w_N <- 0; w_0 <- 0.5; w_1 <- 0.5   # e.g. Balanced (Fleming)

weighted_simulator <- function(theta, fidelity = "high", seed = NULL, ...) {
  # ... run trial replications, then compute the raw operating characteristics:
  EN_null <- ...   # expected sample size under the null
  EN_alt  <- ...   # expected sample size under the alternative
  N_max   <- ...   # maximum sample size
  power   <- ...
  type1   <- ...

  # Simulator computes the weighted loss and RETURNS it as a named metric:
  weighted_EN_N <- w_N * N_max + (1 - w_N) * (w_0 * EN_null + w_1 * EN_alt)

  c(weighted_EN_N = weighted_EN_N,
    EN_null = EN_null, EN_alt = EN_alt, N_max = N_max,
    power = power, type1 = type1)
}

fit <- bo_calibrate(
  sim_fun     = weighted_simulator,
  bounds      = bounds,
  objective   = "weighted_EN_N",                 # name of a returned metric
  constraints = list(power = c("ge", 0.80),      # constraint names must also
                     type1 = c("le", 0.10)),     # be returned-metric names
  budget = 300, seed = 2025)
```

To sweep philosophies, vary `w_N`/`w_0`/`w_1` per the table above and re-run
(this is what `bo_calibrate_philosophies()` automates over a grid of weights).

## Multi-Stage Warm-Start Workflow

```r
# Define broad search region and constraints
broad_bounds <- list(threshold = c(0.15, 0.50), nmax = c(20, 100))
constraints  <- list(power = c("ge", 0.80), type1 = c("le", 0.10))

# Stage 1: Broad exploration
fit1 <- bo_calibrate(sim_fun = my_sim, bounds = broad_bounds,
  objective = "EN", constraints = constraints, budget = 300, seed = 2025)
save_bo_state(fit1, "stage1.rds")

# Narrow bounds around the best region
narrow_bounds <- refine_bounds(fit1, shrink_factor = 0.5)

# Stage 2: Warm start with Stage 1 history
fit2 <- bo_calibrate(sim_fun = my_sim, bounds = narrow_bounds,
  objective = "EN", constraints = constraints,
  initial_history = fit1$history, budget = 200, seed = 2026)
```

Additional utilities: `fix_parameters()`, `remove_fixed_from_bounds()`,
`sequential_refinement()`, `load_bo_state()`.

## Multi-Fidelity Optimization

Four fidelity selection methods via `fidelity_method`:

| Method | Description |
|--------|-------------|
| `"adaptive"` (default) | Cost-aware, balances information gain vs. compute |
| `"staged"` | Fixed iteration-based schedule |
| `"threshold"` | Feasibility probability thresholds |
| `"hybrid_staged"` | MCEM-inspired with CV-based triggers |

```r
fit <- bo_calibrate(...,
  fidelity_levels = c(low = 500, med = 2000, high = 10000),
  fidelity_costs  = c(low = 1, med = 3, high = 15),
  fidelity_method = "adaptive")
```

## Working with Results

`bo_calibrate()` returns a `BATON_fit` object with these fields:

| Field | Description |
|-------|-------------|
| `history` | All evaluations (eval_id, iter, theta, metrics, variance, n_rep, objective, feasible, fidelity, acq_score, prob_feas, cv_estimate) |
| `best_theta` | Named list of best feasible design parameters |
| `best_objective` | Minimum objective among feasible designs |
| `surrogates` | Fitted GP models for each metric |
| `policies` | Acquisition, fidelity, and budget configuration |
| `diagnostics` | Posterior draws and sensitivity info |
| `bounds` | Parameter bounds |
| `constraints` | Constraint specification |
| `constraint_tbl` | Parsed constraint tibble |

There is no `best_metrics` field. To retrieve metrics for the best design:

```r
feasible <- fit$history[fit$history$feasible, ]
best_row <- feasible[which.min(feasible$objective), ]
```

Convergence check: `check_convergence(fit)`.

## Benchmarking and Sensitivity Analysis

```r
# Compare BO vs grid vs random vs heuristic
benchmark_methods(sim_fun = my_sim, bounds = bounds,
                  objective = "EN", constraints = constraints, budget = 100,
                  grid_args = list(resolution = 10))

# Sensitivity analysis (high-level, accepts BATON_fit directly)
analyze_parameter_importance(fit)

# Sensitivity analysis (low-level, requires components from fit)
sa_sobol(fit$surrogates, bounds = bounds, outcome = "EN")
sa_gradients(fit$surrogates, theta = fit$best_theta, bounds = bounds, outcome = "EN")
extract_lengthscales(fit$surrogates)
cov_effects(fit$surrogates, bounds = bounds, outcome = "EN")

# Case study tools
summarise_case_study(fit)
case_study_diagnostics(fit)
```

## Package Structure

```
R/
  ablation.R             # Multi-fidelity ablation studies
  acquisition.R          # ECI acquisition and batch selection
  benchmark.R            # Benchmark comparisons
  bo_calibrate.R         # Main optimization loop
  bo_parameter_importance.R
  bo_warmstart.R         # save/load state, refine_bounds, fix_parameters
  bounds.R               # Bound manipulation utilities
  case_study.R           # Summaries and diagnostics
  constraints.R          # Constraint parsing and feasibility
  convergence.R          # Convergence diagnostics
  init_stopping.R        # GP-based initialization early stopping
  reliability.R          # Constraint reliability estimation
  sensitivity.R          # Sobol, gradients, covariance analysis
  surrogates.R           # Heteroskedastic GP fitting
  utils.R                # Plotting (13 functions) and utilities
  welford.R              # Welford variance estimation and pooling
inst/examples/
  run_demo.R
  simulator_with_variance.R
vignettes/
  BATON-introduction.Rmd
  advanced-features.Rmd
  variance-estimation.Rmd
tests/testthat/          # 8 test files
```

## Vignettes and Documentation

```r
vignette("BATON-introduction")  # Overview
vignette("advanced-features")   # Multi-stage, warm-start, fidelity control
vignette("variance-estimation") # Welford's algorithm

?bo_calibrate           ?welford_mean_var       ?benchmark_methods
?save_bo_state          ?refine_bounds          ?sa_sobol
```

## Testing

```r
devtools::test()
testthat::test_file("tests/testthat/test-BATON-core.R")
```

## Citation

```bibtex
@article{rashid2026baton,
  title   = {Constrained {Bayesian} Optimization for Calibration of
             {Bayesian} Adaptive Clinical Trials},
  author  = {Rashid, Naim},
  journal = {Journal of the American Statistical Association},
  year    = {2026},
  note    = {R package version 0.4.0, \url{https://github.com/naimurashid/BATON}}
}
```

## License

Copyright (c) 2026, The University of North Carolina at Chapel Hill. For
not-for-profit research and educational use only; all other rights reserved.
See LICENSE. For commercial licensing, contact otc@unc.edu.

## Contact

**Naim Rashid** - Department of Biostatistics, UNC Chapel Hill; Lineberger Comprehensive Cancer Center
- Email: naim_rashid@unc.edu
- Issues: <https://github.com/naimurashid/BATON/issues>
