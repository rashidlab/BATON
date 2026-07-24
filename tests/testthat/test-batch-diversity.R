# Regression tests for local-penalization batch diversity (review item 1.2 / Task A3).
# Before the fix the penalty grew WITH distance, so batches clustered around the
# first pick (the opposite of Gonzalez et al. 2016).

test_that("local penalization spreads the batch (suppresses near-duplicates)", {
  # scores decreasing left to right; first pick is candidate 1 (x = 0.00).
  # Candidate 2 (x = 0.01) is a near-duplicate; candidate 3 (x = 0.99) is far.
  cand <- list(c(x = 0.00), c(x = 0.01), c(x = 0.99))
  scores <- c(1.0, 0.95, 0.80)
  sel <- BATON:::select_batch_local_penalization(cand, scores, q = 2, lipschitz = 10)
  # A diverse batch must take the FAR point (index 3), not the near-duplicate (2).
  expect_equal(sort(sel), c(1L, 3L))
})

test_that("penalty is strongest at the selected point and fades with distance", {
  # Three candidates at increasing distance from the top-scoring one; with equal
  # base scores the batch should walk outward, never re-pick the near neighbour.
  cand <- list(c(x = 0.50), c(x = 0.52), c(x = 0.95))
  scores <- c(1.0, 0.99, 0.98)
  sel <- BATON:::select_batch_local_penalization(cand, scores, q = 2, lipschitz = 10)
  expect_equal(sort(sel), c(1L, 3L))  # far point wins over the 0.02-away duplicate
})

test_that("q >= n returns all candidates unchanged", {
  cand <- list(c(x = 0.1), c(x = 0.2))
  expect_equal(BATON:::select_batch_local_penalization(cand, c(0.5, 0.4), q = 2),
               c(1L, 2L))
})
