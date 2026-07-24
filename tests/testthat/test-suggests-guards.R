# Task D9: analysis-only dependencies (ggplot2, randtoolbox) live in
# Suggests. Every exported plot function must fail with an informative
# install hint (via require_suggests) when ggplot2 is absent.

test_that("plot functions error informatively, not cryptically, without ggplot2", {
  # simulate absence via mocked requireNamespace is brittle; instead assert
  # every exported plot_* function guards ggplot2 with require_suggests().
  # plot_case_sobol and plot_case_gradient delegate to plot_sobol_indices /
  # plot_gradient_heatmap (which carry their own first-expression guard), so
  # the looser contains-check is enough for them; every other plot_* function
  # must have the guard as the FIRST expression of its body, before any work.
  delegators <- c("plot_case_sobol", "plot_case_gradient")
  for (fn_name in grep("^plot_", getNamespaceExports("BATON"), value = TRUE)) {
    fn_body <- body(get(fn_name, asNamespace("BATON")))
    if (fn_name %in% delegators) {
      body_txt <- paste(deparse(fn_body), collapse = " ")
      expect_match(body_txt, "require_suggests", info = fn_name)
    } else {
      first_expr <- fn_body[[2]]
      expect_true(
        is.call(first_expr) &&
          identical(first_expr[[1]], as.name("require_suggests")),
        info = fn_name
      )
    }
  }
})

test_that("require_suggests errors with an install hint for a missing package", {
  expect_error(
    BATON:::require_suggests("nonexistent.package.xyz"),
    "install.packages",
    fixed = TRUE
  )
  expect_silent(BATON:::require_suggests("stats"))
})

test_that("require_suggests names the calling function in its error", {
  needs_missing_pkg <- function() BATON:::require_suggests("nonexistent.package.xyz")
  expect_error(needs_missing_pkg(), "needs_missing_pkg\\(\\)")
})

test_that("init-stopping test grid falls back to random points with a message", {
  # randtoolbox is in Suggests, so the runif fallback in generate_test_grid
  # is live code on hosts without it. It must announce itself: two machines
  # with identical seeds would otherwise silently produce different
  # init-stopping decisions. The fallback branch is exercised directly via
  # the use_sobol seam (mocking base::requireNamespace is unreliable).
  X <- matrix(runif(20), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
  expect_message(
    grid <- BATON:::generate_test_grid(X, n = 50, use_sobol = FALSE),
    "randtoolbox not installed"
  )
  expect_equal(dim(grid), c(50L, 2L))
  expect_identical(colnames(grid), c("x1", "x2"))

  # With randtoolbox available, the Sobol path stays silent.
  skip_if_not_installed("randtoolbox")
  expect_silent(grid_sobol <- BATON:::generate_test_grid(X, n = 50, use_sobol = TRUE))
  expect_equal(dim(grid_sobol), c(50L, 2L))
})
