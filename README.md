# BATON

**Bayesian Optimization for Calibration of Adaptive Clinical Trials**

Version 0.7.1 | R >= 4.2 | UNC not-for-profit license | [GitHub](https://github.com/rashidlab/BATON) | [Issues](https://github.com/rashidlab/BATON/issues)

Designing a clinical trial requires choosing design parameters -- for
example, the maximum sample size, efficacy and futility decision thresholds,
or the number of interim analyses -- so that the trial achieves adequate
power while controlling type I error, and simultaneously minimizes a
quantity of interest such as the expected sample size. These performance
metrics are collectively called **operating characteristics**. For adaptive
designs, operating characteristics are often impossible to compute
analytically and must be estimated by simulating thousands of trial
replicates.

BATON automates this calibration. Given a simulator that maps design
parameters to operating characteristics, BATON uses constrained Bayesian
optimization to find configurations that minimize a user-specified objective
(e.g., expected sample size) subject to constraints (e.g., power >= 0.80,
type I error <= 0.10), without requiring closed-form gradients or analytic
solutions. It works with any trial simulator, Bayesian or frequentist, and
is most useful when analytic solutions do not exist and the design space is
too large for grid search. Under the hood, BATON fits a Gaussian process
surrogate (a fast statistical stand-in for your expensive simulator) to each
operating characteristic, weights each observation by its Monte Carlo noise,
scores candidate designs by expected constrained improvement (how much a
candidate is likely to improve the objective while still satisfying the
constraints), and spends cheap low-replication simulations on exploration
while reserving high-replication runs for verification.

## Releases

| Version | Date | Headline |
|---------|------|----------|
| 0.7.1 | Jul 2026 | Documentation release: self-contained Getting Started, new methods and case-study vignettes |
| 0.7.0 | Jul 2026 | Service controls: `status` field, `max_walltime_s`, `callback` cancellation, `checkpoint_fun`, `on_error = "return_partial"`, `slim`; slimmer dependencies |
| 0.6.0 | Jul 2026 | Optional `n_rep` simulator contract; parallel evaluation via `options(BATON.cores)`; matrix-based hot paths (58% faster candidate scoring) |
| 0.5.0 | Jul 2026 | Correctness release from a full code review: working warm-start, corrected batch penalization, robustness fixes |
| 0.4.0 | May 2026 | `warmstart_from` seeds, `multi_seed_verify` gate, `bo_calibrate_philosophies()` batch wrapper |

Full details and migration notes: [NEWS.md](NEWS.md).

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
remotes::install_github("rashidlab/BATON")

# With vignettes (requires knitr, rmarkdown):
remotes::install_github("rashidlab/BATON", build_vignettes = TRUE)
```

BATON requires R >= 4.2 (see `DESCRIPTION`). The companion simulator package
[evolveTrial](https://github.com/naimurashid/evolveTrial) provides ready-made
oncology trial simulators, but it is optional and not a dependency: BATON
works with any evaluator you write yourself (see below).

### Troubleshooting

| Error | Fix |
|-------|-----|
| `gfortran: command not found` (during `DiceKriging` install) | Install Fortran per the System Dependencies table above (macOS: `brew install gcc`; Linux: install `gfortran`) |
| `ld: library not found for -lgfortran` (macOS) | `brew install gcc` and ensure the gfortran path is on `PATH` (see <https://mac.r-project.org/tools/>) |
| `hetGP` fails to link BLAS/LAPACK (Linux) | `sudo apt-get install libblas-dev liblapack-dev` (Ubuntu/Debian) or `sudo dnf install blas-devel lapack-devel` (Fedora/RHEL) |
| `make: not found` (Windows) | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version; verify with `Sys.which("make")` |
| Runs too slow on one core | Parallelize the Monte Carlo loop inside your evaluator (see "Parallel Evaluation Within Your Evaluator" below), or set `options(BATON.cores = k)` to evaluate batches in parallel (Unix/macOS only) |
| Installation times out on CI or slow links | `remotes::install_github("rashidlab/BATON", build_vignettes = FALSE)` |

## Quick Start

Self-contained example using a synthetic simulator. Copy-paste and run after
installing BATON.

```r
library(BATON)

# Step 1: Define a simulator that evaluates under BOTH hypotheses
# The evaluator runs n_rep simulated trials under null (p = 0.2) and
# alternative (p = 0.4), returning power, type I error, and expected N.
# BATON calls it with a fidelity label; this simulator maps the label to a
# replication count itself. (BATON can also pass the exact requested count;
# see "Evaluator Requirements" below.)
my_sim <- function(theta, fidelity = "high", seed = NULL, ...) {
  # Toy counts, kept small so this first run finishes in seconds; the other
  # examples sync their fallbacks to fidelity_levels (see Evaluator
  # Requirements for the n_rep contract).
  n_rep <- switch(fidelity, low = 500, med = 2000, high = 5000)
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
fit$status           # why the run ended

# Get metrics for the best design (no fit$best_metrics field exists)
feasible <- fit$history[fit$history$feasible, ]
best_row <- feasible[which.min(feasible$objective), ]
best_row[, c("power", "type1", "EN")]
```

```
# Example output (abridged; a few seconds total on a laptop)
Initialising BATON calibration with seed = 42
Fidelity selection method: 'adaptive'
  Eval 001 (iter 00, low fidelity): EN = 66.0000
  Eval 002 (iter 00, low fidelity): EN = 72.0000 [feasible]
  ...
Iteration 1: fitting surrogates on 20 evaluations.
  -> max acquisition score: 26.618
  Eval 021 (iter 01, low fidelity): EN = 30.0000
  ...
  -> best observed feasible: 30.000 (change: 0.000)
Early stopping at iteration 9: no improvement > 0.1% for 2 consecutive iterations

> fit$best_theta
$threshold
[1] 0.3118077
$nmax
[1] 30.09349
> fit$best_objective
[1] 30
> fit$status
[1] "early_stopped"
> best_row[, c("power", "type1", "EN")]
  power type1    EN
1 0.832 0.052    30
```

## Minimal Template

To calibrate your own design, you replace exactly one function: the code
that simulates a single trial. The adapter below handles everything BATON
needs. Metric names are arbitrary (`power`/`type1`/`EN` are conventions of
this README's examples, not requirements); `objective` and `constraints`
just have to refer to names your evaluator returns.

```r
library(BATON)

# ---- 1. YOUR TRIAL LOGIC (replace this) -------------------------------
# Simulate ONE trial at design parameters theta; return named metrics.
sim_one_trial <- function(theta) {
  # ... your design logic here ...
  stop("REPLACE: simulate one replicate of YOUR design")
  # c(power = ..., type1 = ..., EN = ...)   # any names you like
}

# ---- 2. BATON ADAPTER (usually keep as-is) ----------------------------
my_evaluator <- function(theta, fidelity = c("low", "med", "high"),
                         seed = NULL, n_rep = NULL, ...) {
  fidelity <- match.arg(fidelity)
  if (is.null(n_rep)) {                      # fallback if BATON did not pass n_rep
    n_rep <- switch(fidelity, low = 2000, med = 4000, high = 10000)
  }
  if (!is.null(seed)) set.seed(seed)
  result <- welford_mean_var(function(i, theta) sim_one_trial(theta),
                             n_samples = n_rep, theta = theta)
  metrics <- result$mean
  attr(metrics, "variance") <- result$variance
  attr(metrics, "n_rep") <- result$n
  metrics
}

# ---- 3. CALIBRATE -----------------------------------------------------
fit <- bo_calibrate(
  sim_fun     = my_evaluator,
  bounds      = list(...),        # named list: one c(lower, upper) per parameter
  objective   = "EN",             # a metric your evaluator returns
  constraints = list(power = c("ge", 0.80), type1 = c("le", 0.10)),
  n_init = 20, budget = 60, seed = 1
)
```

Fill in `sim_one_trial()`, set `bounds` to your design parameters, and run.
The next section tells you how to size each setting.

## Choosing Your Settings

This table is the single source of truth for sizing a run. `d` is the number
of design parameters.

| Setting | Default | How to choose |
|---------|---------|---------------|
| `n_init` | 90 | Space-filling evaluations before BO starts. 10 x `d` for a working run; 4-5 x `d` for a quick pass (see profiles below). **Warning:** on an expensive simulator, omitting `n_init` means 90 initial evaluations before the first BO iteration. |
| `budget` | 300 | Total evaluations including `n_init`. 20-30 x `d` for well-behaved problems; add 50-100% when constraints are tight (narrow feasible region). Same warning: the default is 300. |
| `q` | 2 | Evaluations proposed per BO iteration. 2-4; larger batches give wall-clock parallelism at a small cost in statistical efficiency. |
| `bounds` | required | Bracket the plausible optimum generously; BATON rescales everything to the unit cube internally, so wide bounds do not break scaling (they just spend more budget exploring). If the best design lands on a bound, widen it and re-run. Bounds are independent boxes: for interdependent parameters (an interim size that must stay below the total, a stage-1 threshold below the final one), reparameterize so each parameter is free on its own range, e.g., calibrate `interim_fraction` in [0.3, 0.7] instead of `n_interim`, or a threshold *gap* instead of the second threshold. |
| `constraints` | required | Leave slack larger than the Monte Carlo error of your simulator. For a proportion, SE = sqrt(p(1-p)/n_rep): estimating type I error near 0.10 with n_rep = 2000 gives SE = 0.0067, so estimates within about +/- 0.013 (2 SE) of the threshold are indistinguishable from the boundary. Verify near-boundary designs at higher fidelity. |
| `fidelity_levels` | `c(low = 2000, med = 4000, high = 10000)` | Replications per fidelity tier. Size so the low tier is cheap enough for exploration (it absorbs most evaluations) and the high tier pushes the SE below the constraint slack you care about. |
| `integer_params` | `NULL` | Names of parameters (e.g., `"nmax"`) that BATON rounds to integers before calling your simulator. Alternative to rounding inside the simulator, and makes `best_theta` integer-valued for those parameters. |
| `seed` | 2025 | Run seed. BATON derives per-evaluation seeds from it and passes them to your simulator; your RNG stream is restored on exit. |

### Run Profiles: Quick, Balanced, Thorough

| Setting | Quick | Balanced | Thorough |
|---------|-------|----------|----------|
| `n_init` | 4-5 x d | 10 x d | 10 x d |
| `budget` | 10-15 x d | 20-30 x d | 40-50 x d or more |
| `fidelity_levels` | `c(low = 200, med = 500, high = 2000)` | default | raise high tier (e.g., `high = 50000`) |
| `q` | 2-4 | 2 | 2 |
| `early_stop` | default | default | `list(enabled = FALSE)` or tightened `threshold` |
| `multi_seed_verify` | `FALSE` | `TRUE` | `TRUE` with `multi_seed_n = 10` |
| Use when | Debugging your evaluator; first look at a new design space; checking feasibility is achievable at all; iterating on bounds and constraints | The run whose answer you intend to act on; overnight or lunch-scale runs | Manuscript or regulatory numbers; constraints within ~2 SE of a boundary (e.g., type I error near its threshold); d >= 6. Consider `options(BATON.cores = k)` on Unix/macOS. |
| Expect | Minutes. Operating characteristics good to roughly +/- 2-3 percentage points at the low tier (SE = sqrt(0.8 x 0.2 / 200) = 0.028) | Hours, simulator-dependent. SE at the high tier around 0.003-0.005 | Longest runs; SE at `high = 50000` around 0.001-0.002 |

Move between profiles as your problem firms up: use Quick while you iterate
on the problem formulation (bounds, constraints, metric definitions), switch
to Balanced once the formulation is stable, and reserve Thorough for the run
whose numbers will be quoted.

## Writing Your Own Evaluator

BATON works with **any** trial simulator -- whether the trial uses Bayesian
decision rules, frequentist group sequential boundaries, or hybrid designs.
You supply an evaluator function that simulates trials and returns their
operating characteristics; BATON handles the optimization. This section
walks through building one from scratch.

### Evaluator Requirements

The canonical signature is:

```
function(theta, fidelity, seed = NULL, n_rep = NULL, ...)
```

1. **`theta`**: named list of design parameters (e.g., `theta$alpha`,
   `theta$nmax`). Values arrive as-is unless listed in `integer_params`, in
   which case BATON rounds them before the call (the alternative to rounding
   inside your simulator).
2. **`fidelity`**: label from `c("low", "med", "high")` naming the requested
   fidelity tier.
3. **`n_rep`**: if your function declares an `n_rep` argument, BATON passes
   the exact requested replication count, `fidelity_levels[[fidelity]]`,
   including any dynamic escalation. A hardcoded
   `switch(fidelity, ...)` is only a fallback for calls made outside BATON;
   if you use one, keep it in sync with your `fidelity_levels` or the budget
   accounting will not match what your simulator actually did.
4. **`seed`**: always passed by BATON (derived from the run seed); apply it
   with `set.seed(seed)` when non-NULL for reproducibility.
5. **Return value:** named numeric vector of operating characteristics
   (e.g., `c(power = 0.83, type1 = 0.07, EN = 52)`). The names must match
   the names used in `objective` and `constraints`; beyond that they are
   arbitrary.
6. **Variance attributes (recommended):** attach `attr(result, "variance")`
   and `attr(result, "n_rep")` to enable the noise-aware GP surrogate, which
   improves efficiency by 30-50%.

### Complete Example: Group Sequential Design

This evaluator simulates a two-stage group sequential trial with a continuous
endpoint. It is entirely self-contained -- no external packages required beyond
base R.

```r
# Evaluator for a two-stage group sequential trial
# Design parameters: alpha_spend (Stage 1 alpha), nmax (max per arm)
# Fixed: effect size = 0.3, interim at 50% enrollment
gs_evaluator <- function(theta, fidelity = c("low", "med", "high"),
                         seed = NULL, n_rep = NULL, ...) {
  fidelity <- match.arg(fidelity)
  if (is.null(n_rep)) {  # fallback; keep in sync with your fidelity_levels
    n_rep <- switch(fidelity, low = 2000, med = 4000, high = 10000)
  }
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

**Fidelity:** Use the `fidelity` argument (or the passed `n_rep`) to set the
replication count; BATON's multi-fidelity engine schedules low fidelity for
exploration and high fidelity for verification automatically.

**Seeds:** BATON always passes `seed`; apply it with `set.seed(seed)` so
evaluations are reproducible.

**Welford variance:** Use `welford_mean_var()` rather than storing all
replicates: mean and variance of each metric in a single pass, with O(m)
memory regardless of the number of replications.

**Existing simulators:** If you already have a simulator that returns a
data.frame or list, write a thin wrapper that extracts the metrics into a
named vector and attaches variance attributes.

See `system.file("examples/simulator_with_variance.R", package = "BATON")`
for additional examples including parallel execution.

### Variance Attributes and Surrogate Selection

BATON uses noise-aware heteroskedastic Gaussian processes (hetGP), which let
the model trust precise observations more than noisy ones, when
per-evaluation variance estimates are available, and falls back to
constant-noise GPs (DiceKriging) when they are not. The surrogate choice is
made **per metric**:

| Scenario | Surrogate Used | How to Enable |
|----------|---------------|---------------|
| Evaluator attaches `attr(result, "variance")` for metric | **hetGP** (recommended) | Use `welford_mean_var()` in your evaluator |
| No variance attribute, metric in [0,1] (e.g., power, type I) | **hetGP** via binomial approximation | Automatic: BATON estimates `p(1-p)/n_rep` |
| No variance attribute, metric outside [0,1] (e.g., EN, ET) | **DiceKriging** (homoskedastic) | Provide variance attribute to upgrade |
| Analytical evaluator (no MC noise) | **DiceKriging** with small nugget | Appropriate: no noise to model |

**Why this matters:** Monte Carlo variance of expected sample size (EN) is
input-dependent: aggressive early stopping produces bimodal N distributions
(high variance), while designs that rarely stop early concentrate N near
N_max (low variance). A constant-noise GP oversmooths the former and
undersmooths the latter. In our benchmarks, providing EN variance improved
the best feasible design by approximately 15%.

**What to do if you cannot compute per-replication variance:**

If your simulator only returns aggregated means (e.g., from an external
program that reports `E[N]` but not individual trial sample sizes), BATON
still works. Constraint surrogates (power, type I error) automatically get
binomial variance estimates and use hetGP. Only the objective surrogate (EN)
falls back to DiceKriging. This is suboptimal but not incorrect. Options:

1. **Modify the simulator** to store per-replication values and compute
   `var(N_per_rep) / n_rep`. This is the recommended approach.
2. **Use `welford_mean_var()`** to wrap the per-trial simulation function if
   you have access to it. This is exactly the adapter in the Minimal
   Template above.
3. **Accept the fallback** for rapid prototyping: return the aggregated
   means as a plain named vector (no attributes) and increase budget by
   20-50% to compensate for the less efficient objective surrogate.

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
my_sim <- function(theta, fidelity = "high", seed = NULL, n_rep = NULL,
                   ncores = NULL, ...) {
  if (is.null(n_rep)) {  # fallback; keep in sync with your fidelity_levels
    n_rep <- switch(fidelity, low = 2000, med = 4000, high = 10000)
  }
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

# Use all available cores (autodetect); pass ncores = k to pin the count
fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, budget = 60, seed = 42
)
```

**Key points:**
- `parallel::detectCores(logical = FALSE)` returns physical cores; subtract
  1 to leave a core free for the BO loop itself.
- macOS/Linux fork via `mclapply` (no data copying); Windows has no fork,
  hence the `makeCluster`/`parLapply` socket branch.
- `pool_welford_results()` uses Chan's parallel variance algorithm to
  combine chunk means and variances with the same one-pass numerical
  stability as Welford.
- Each chunk gets `ceiling(n_rep / ncores)` replications, so the total may
  slightly exceed `n_rep`. This is harmless.

## Before You Trust the Answer

Three checks before a calibrated design goes anywhere important:

**1. Verify at high fidelity.** Most evaluations during a run happen at low
fidelity. Set `multi_seed_verify = TRUE` to re-evaluate the best design at
multiple seeds at high fidelity after the run:

```r
fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, budget = 60, seed = 42,
  multi_seed_verify = TRUE
)
fit$verdict              # "MULTI_SEED_PASS", "MULTI_SEED_WARN", or "MULTI_SEED_FAIL"
fit$multi_seed_summary   # cross-seed means and sds for each metric
```

`MULTI_SEED_PASS`: feasible at every verification seed. `MULTI_SEED_FAIL`:
at least one seed violated a constraint (default `multi_seed_strict = TRUE`).
`MULTI_SEED_WARN`: same failure with `multi_seed_strict = FALSE`.
Alternatively, verify manually: re-run your simulator at `fit$best_theta`
with a large `n_rep` and confirm the constraints hold.

**2. Round integer parameters.** `best_theta` comes back fractional (e.g.,
`nmax = 30.09`) unless you passed `integer_params`. Round to the design you
would actually run, then re-verify the operating characteristics at the
rounded design.

**3. Check `fit$status`.** It records why the run ended:

| `status` | Meaning |
|----------|---------|
| `budget_exhausted` | Normal completion: the evaluation budget was consumed |
| `early_stopped` | Converged: no meaningful improvement over the patience window |
| `acq_flatline` | Acquisition values collapsed; no candidate looks promising |
| `walltime` | The `max_walltime_s` cap tripped; result is whatever was found in time |
| `cancelled` | Your `callback` requested cancellation |
| `errored` | A failure was caught under `on_error = "return_partial"`; see `fit$error_message` |

## Service and Long-Run Controls (v0.7.0)

For unattended, scheduled, or cluster runs, `bo_calibrate()` has controls
for capping walltime, checkpointing, surviving simulator failures, and
cooperative cancellation. Details for each are in `?bo_calibrate`.

```r
options(BATON.cores = 4)  # parallel simulator evaluation (Unix/macOS; serial on Windows)

fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, n_init = 20, budget = 60, seed = 42,
  max_walltime_s = 6 * 3600,    # graceful stop before the scheduler kills the job
  checkpoint_fun = function(snap) saveRDS(snap, "checkpoint.rds"),
  checkpoint_every = 5,         # snapshot every 5 BO iterations
  on_error = "return_partial",  # keep completed evaluations if the simulator fails
  callback = function(info) !file.exists("STOP"),  # return TRUE to continue
  slim = TRUE                   # lightweight return: skips the final surrogate refit
)

# Resume a stopped run from its checkpoint
ck <- readRDS("checkpoint.rds")
fit <- bo_calibrate(
  sim_fun = my_sim, bounds = bounds, objective = "EN",
  constraints = constraints, initial_history = ck$history,
  n_init = nrow(ck$history), budget = nrow(ck$history) + 20, seed = 43
)
```

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

To sweep philosophies, vary `w_N`/`w_0`/`w_1` per the table above and re-run.
`bo_calibrate_philosophies()` orchestrates such a batch: it runs a named set
of philosophies in dependency order (donor philosophies first), warm-starts
recipients from donor solutions, and attaches a multi-seed verification
verdict to each fit. You still supply each philosophy's `objective` and
`constraints`; it does not construct weighted objectives for you.

## Multi-Stage Warm-Start Workflow

Recommended for `d >= 4`: run a broad Stage 1, then narrow bounds around the
best region and warm-start Stage 2 with `initial_history`.

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
| `status` | Why the run ended (see "Before You Trust the Answer") |
| `error_message` | Failure message when `status = "errored"`; otherwise `NULL` |
| `multi_seed_summary`, `multi_seed_runs`, `verdict` | Stage 4 verification results when `multi_seed_verify = TRUE` (`multi_seed_runs` is `NULL` under `slim = TRUE`) |
| `surrogates` | Fitted GP models for each metric (`NULL` under `slim = TRUE`) |
| `policies` | Acquisition, fidelity, and budget configuration |
| `diagnostics` | Posterior draws and sensitivity info (`NULL` under `slim = TRUE`) |
| `bounds` | Parameter bounds |
| `constraints` | Constraint specification |
| `constraint_tbl` | Parsed constraint tibble |

There is no `best_metrics` field. To retrieve metrics for the best design:

```r
feasible <- fit$history[fit$history$feasible, ]
best_row <- feasible[which.min(feasible$objective), ]
```

Convergence check: `check_convergence(fit$history)`.

## Benchmarking and Sensitivity Analysis

```r
# Compare BO vs grid vs random vs heuristic.
# bo_calibrate() settings go in bo_args, not top-level arguments.
benchmark_methods(sim_fun = my_sim, bounds = bounds,
                  objective = "EN", constraints = constraints,
                  bo_args = list(n_init = 20, budget = 100),
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
  bo_calibrate_philosophies.R  # Batch wrapper across design philosophies
  bo_parameter_importance.R
  bo_v04_helpers.R       # Warm-start donors and multi-seed verification
  bo_warmstart.R         # save/load state, refine_bounds, fix_parameters
  bounds.R               # Bound manipulation utilities
  case_study.R           # Summaries and diagnostics
  constraints.R          # Constraint parsing and feasibility
  convergence.R          # Convergence diagnostics
  init_stopping.R        # GP-based initialization early stopping
  reliability.R          # Constraint reliability estimation
  sensitivity.R          # Sobol, gradients, covariance analysis
  surrogates.R           # Heteroskedastic GP fitting
  utils.R                # Scaling, hashing, misc utilities
  welford.R              # Welford variance estimation and pooling
inst/examples/
  run_demo.R
  simulator_with_variance.R
vignettes/
  BATON-introduction.Rmd   # Getting started (self-contained)
  BATON-methods.Rmd        # How the optimizer works
  BATON-case-study.Rmd     # End-to-end workflow with the analysis helpers
  advanced-features.Rmd
  v04-hardening.Rmd
  variance-estimation.Rmd
tests/testthat/          # 32 test files
```

Plotting: 14 exported `plot_*` functions live alongside the analyses they
visualize (benchmark, sensitivity, case study, reliability, ablation,
initialization).

## Vignettes and Documentation

```r
vignette("BATON-introduction")  # Getting started: evaluator, settings, trust checklist
vignette("BATON-methods")       # How it works: GP surrogates, ECI, multi-fidelity, seeding
vignette("BATON-case-study")    # End to end: calibrate, verify, sensitivity, benchmark
vignette("advanced-features")   # Multi-stage, warm-start, fidelity control
vignette("v04-hardening")       # Cross-philosophy warm-start, multi-seed verification
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
@unpublished{young2026baton,
  title  = {{BATON}: Constrained {Bayesian} Optimization for Calibrating
            Adaptive Trials},
  author = {Young, A. M. and Li, D. and Hilsenbeck, S. G. and Tayob, N. and
            Chen, R. and Yuan, Y. and Bates, S. and Kelly, E. and Burns, R. and
            Jhaveri, K. and Tolaney, S. M. and Spears, P. A. and
            Goetz, M. P. and Davidson, N. E. and Norton, L. and
            Perou, C. M. and Krop, I. E. and Wolff, A. C. and Winer, E. P. and
            Carey, L. A. and Rashid, N. U.},
  note   = {Submitted. R package version 0.7.1,
            \url{https://github.com/rashidlab/BATON}},
  year   = {2026}
}
```

## License

Copyright (c) 2026, The University of North Carolina at Chapel Hill. For
not-for-profit research and educational use only; all other rights reserved.
See LICENSE. For commercial licensing, contact otc@unc.edu.

## Contact

**Naim Rashid** - Department of Biostatistics, UNC Chapel Hill; Lineberger Comprehensive Cancer Center
- Email: naim_rashid@unc.edu
- Issues: <https://github.com/rashidlab/BATON/issues>
