# Tests for candidate-pool sizing (v0.8 defect 2).
#
# The late-run +50 percent pool inflation must be reachable: it thresholds on
# the ACTUAL iteration budget (budget - n_init) / q, not budget / q. Under the
# package defaults (n_init = 90, budget = 300, q = 1) the loop runs at most 210
# iterations, so the old threshold 0.7 * (budget / q) = 210 could never be
# exceeded and the inflation was dead code.

test_that("late inflation fires in the last 30 percent of a default-shaped run", {
  d <- 2  # base pool = pmax(1000, pmin(5000, 500 * d)) = 1000
  # Defaults: n_init 90, budget 300, q 1 -> iteration budget 210.
  # Threshold 0.7 * 210 = 147: iterations 148..210 get the 1.5x pool.
  expect_equal(candidate_pool_size_for_iter(147, d = d, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 1000)
  expect_equal(candidate_pool_size_for_iter(148, d = d, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 1500)
  expect_equal(candidate_pool_size_for_iter(210, d = d, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 1500)
})

test_that("inflation threshold accounts for batch size q", {
  # n_init 90, budget 300, q 2 -> 105 iterations; threshold 73.5.
  expect_equal(candidate_pool_size_for_iter(73, d = 2, budget = 300,
                                            n_init = 90, q = 2,
                                            candidate_pool = 1), 1000)
  expect_equal(candidate_pool_size_for_iter(74, d = 2, budget = 300,
                                            n_init = 90, q = 2,
                                            candidate_pool = 1), 1500)
})

test_that("non-divisible budgets use the ceiling iteration count", {
  # budget - n_init = 11, q = 2: the final batch is partial, so the loop runs
  # ceiling(11 / 2) = 6 iterations. Threshold 0.7 * 6 = 4.2: inflation must
  # start at iteration 5, not at iteration 4 (raw division would give
  # threshold 0.7 * 5.5 = 3.85 and inflate one iteration early).
  expect_equal(candidate_pool_size_for_iter(4, d = 2, budget = 16,
                                            n_init = 5, q = 2,
                                            candidate_pool = 1), 1000)
  expect_equal(candidate_pool_size_for_iter(5, d = 2, budget = 16,
                                            n_init = 5, q = 2,
                                            candidate_pool = 1), 1500)
  expect_equal(candidate_pool_size_for_iter(6, d = 2, budget = 16,
                                            n_init = 5, q = 2,
                                            candidate_pool = 1), 1500)
})

test_that("base pool scales with dimension and is clamped", {
  expect_equal(candidate_pool_size_for_iter(1, d = 1, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 1000)
  expect_equal(candidate_pool_size_for_iter(1, d = 6, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 3000)
  expect_equal(candidate_pool_size_for_iter(1, d = 20, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 5000)
  # Inflated high-dimensional pool: 5000 * 1.5 = 7500 (below the 10000 cap).
  expect_equal(candidate_pool_size_for_iter(210, d = 20, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 1), 7500)
})

test_that("user candidate_pool floor is respected before and after inflation", {
  # Floor above both base (1000) and inflated (1500) sizes.
  expect_equal(candidate_pool_size_for_iter(1, d = 2, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 2000), 2000)
  expect_equal(candidate_pool_size_for_iter(210, d = 2, budget = 300,
                                            n_init = 90, q = 1,
                                            candidate_pool = 2000), 2000)
})

test_that("golden-run shape keeps its pool unchanged by the threshold fix", {
  # Golden master: n_init 6, budget 16, q 2, d 2, default candidate_pool 2000.
  # Iteration budget 5, threshold 3.5: iterations 4-5 now inflate 1000 -> 1500,
  # but the 2000 floor dominates both, so the realized pool (and therefore the
  # golden trajectory) is identical before and after the fix.
  for (it in 1:5) {
    expect_equal(candidate_pool_size_for_iter(it, d = 2, budget = 16,
                                              n_init = 6, q = 2,
                                              candidate_pool = 2000), 2000)
  }
})
