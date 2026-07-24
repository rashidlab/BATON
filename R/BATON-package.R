# Package-level imports and non-standard-evaluation declarations.
#
# The importFrom directives cover base-adjacent generics used unqualified
# across the package (R CMD check "no visible global function definition");
# the globalVariables call covers dplyr/ggplot2 column references
# ("no visible binding for global variable").

#' @importFrom methods as
#' @importFrom stats median predict runif sd setNames
NULL

utils::globalVariables(c(
  "eval_id", "feasible", "best_metrics", "n_points", "max_var"
))
