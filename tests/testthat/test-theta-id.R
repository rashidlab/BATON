# Tests for theta_to_id stability (review item 1.12 / Task A6). Duplicate
# detection and cross-session warm-start matching depend on identical unit
# coordinates producing identical ids.

test_that("theta_to_id is stable at bound edges (no scientific-notation flip)", {
  a <- BATON:::theta_to_id(c(0.0, 0.5, 1.0))
  b <- BATON:::theta_to_id(c(1e-7, 0.5, 1.0))   # tiny coord must not reformat 0.5
  expect_true(grepl("0.500000", a, fixed = TRUE))
  expect_true(grepl("0.500000", b, fixed = TRUE))
})

test_that("theta_to_id is independent of session options (scipen)", {
  a <- BATON:::theta_to_id(c(0.0, 0.5, 1.0))
  old <- options(scipen = -10); on.exit(options(old))
  expect_identical(BATON:::theta_to_id(c(0.0, 0.5, 1.0)), a)
})

test_that("identical unit thetas map to identical ids; distinct ones differ", {
  expect_identical(BATON:::theta_to_id(c(0.123456, 0.25)),
                   BATON:::theta_to_id(c(0.123456, 0.25)))
  expect_false(identical(BATON:::theta_to_id(c(0.123456, 0.25)),
                         BATON:::theta_to_id(c(0.123457, 0.25))))
})
