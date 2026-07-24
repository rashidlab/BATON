#' Internal helpers for parameter transformations
#' @keywords internal
scale_to_unit <- function(theta, bounds) {
  # Use base R instead of purrr::imap_dbl to avoid "In index: X" errors
  theta_names <- names(theta)
  result <- vapply(theta_names, function(name) {
    value <- as.numeric(theta[[name]])
    range <- bounds[[name]]
    if (length(range) != 2L) {
      stop(sprintf("Bounds for '%s' must have length 2.", name), call. = FALSE)
    }
    (value - range[[1]]) / (range[[2]] - range[[1]])
  }, FUN.VALUE = numeric(1))
  names(result) <- theta_names
  result
}

#' @keywords internal
scale_from_unit <- function(unit_x, bounds) {
  # Use base R instead of purrr::imap to avoid "In index: X" errors
  param_names <- names(unit_x)
  result <- lapply(param_names, function(name) {
    value <- as.numeric(unit_x[[name]])
    range <- bounds[[name]]
    range[[1]] + value * (range[[2]] - range[[1]])
  })
  names(result) <- param_names
  result
}

#' Sample an initial Latin hypercube design
#' @keywords internal
lhs_design <- function(n, bounds, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  d <- length(bounds)
  lhs <- lhs::randomLHS(n, d)
  theta_names <- names(bounds)
  # Use base R lapply instead of purrr::map
  lapply(seq_len(n), function(i) {
    unit_list <- as.list(lhs[i, , drop = TRUE])
    names(unit_list) <- theta_names
    scaled <- scale_from_unit(unit_list, bounds)
    names(scaled) <- theta_names
    scaled
  })
}

#' Convert theta list to named numeric vector
#' @keywords internal
theta_list_to_vec <- function(theta) {
  # Use base R vapply instead of purrr::map_dbl
  result <- vapply(theta, as.numeric, FUN.VALUE = numeric(1))
  names(result) <- names(theta)
  result
}

#' Create deterministic identifier for theta in unit space
#' @keywords internal
theta_to_id <- function(unit_theta) {
  # formatC(..., format = "f") formats each coordinate independently in fixed
  # notation. base::format() applies a common style to the whole vector, so a
  # single tiny coordinate (e.g. a bound at ~1e-7) flips the entire vector to
  # scientific notation, and the result depends on session options (scipen,
  # OutDec) - both break duplicate detection at bound edges and across sessions.
  paste(formatC(round(as.numeric(unit_theta), digits = 6), format = "f", digits = 6),
        collapse = "|")
}

#' Enforce integer parameters if requested
#' @keywords internal
coerce_theta_types <- function(theta, integer_params = NULL) {
  if (is.null(integer_params) || length(integer_params) == 0L) {
    return(theta)
  }
  # Use base R lapply instead of purrr::imap
  theta_names <- names(theta)
  result <- lapply(theta_names, function(name) {
    value <- theta[[name]]
    if (name %in% integer_params) {
      value <- round(as.numeric(value))
    }
    value
  })
  names(result) <- theta_names
  result
}

# NOTE: validate_bounds() moved to bounds.R with enhanced error messages
# and param_names parameter support

#' Require a suggested package at runtime (Task D9)
#'
#' Analysis-only dependencies (ggplot2, randtoolbox) live in Suggests so the
#' calibration core installs without them. Any exported function that needs
#' one at runtime calls this first, turning a cryptic "there is no package
#' called ..." into an actionable install hint that names the calling
#' function.
#'
#' @param pkg single package name.
#' @keywords internal
require_suggests <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    caller <- sys.call(-1)
    if (is.null(caller)) {
      # called from the top level (no caller to name)
      stop(sprintf(
        "Package '%s' is required for this function. Install it with install.packages(\"%s\").",
        pkg, pkg
      ), call. = FALSE)
    }
    stop(sprintf(
      "Package '%s' is required by %s(). Install it with install.packages(\"%s\").",
      pkg, paste(deparse(caller[[1]]), collapse = ""), pkg
    ), call. = FALSE)
  }
}

# Null-coalescing helper (internal; not exported, no .Rd generated).
# Rd's \name{} macro cannot contain |, so we suppress roxygen docs for the
# infix operator to avoid `checkRd` WARNINGs on R CMD check.
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

#' Effective parallel worker count for forked evaluation
#'
#' Central decision for how many workers `options(BATON.cores = k)` actually
#' buys: forking requires Unix and the parallel package, so anywhere else the
#' effective count is 1 (serial `lapply`). Shared by `evaluate_points()`,
#' `fit_surrogates()`, and the walltime chunk sizing in `bo_calibrate()` so
#' the parallel guard and the walltime-check granularity cannot drift apart
#' (a serial platform gets chunk size 1, i.e. a cap check per evaluation).
#'
#' @return A single integer >= 1.
#' @keywords internal
effective_cores <- function() {
  .cores <- as.integer(getOption("BATON.cores", 1L))
  if (.cores > 1L && .Platform$OS.type == "unix" &&
      requireNamespace("parallel", quietly = TRUE)) {
    .cores
  } else {
    1L
  }
}
