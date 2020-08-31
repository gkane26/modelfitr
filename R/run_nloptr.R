#' run_nloptr
#'
#' runs optimization using nloptr package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param algorithm string; method to use. See details. Default = "neldermead"
#' @param opts list; list of options. See nloptr::nloptr.print.options() for details.
#' @param hessian logical; if TRUE, find the hessian at the optimum. Default = T
#' @param ... further arguments passed to nloptr method and objective
#'
#' @details
#' Algorithm options =
#' "bobyqa",
#' "ccsaq",
#' "cobyla",
#' "crs2lm",
#' "direct",
#' "isres",
#' "lbfgs",
#' "mma",
#' "neldermead",
#' "newuoa",
#' "sbplx",
#' "slsqp",
#' "stogo",
#' "tnewton",
#' "varmetric"
#'
#' @export
run_nloptr <- function(objective, start, lower = NULL, upper = NULL, algorithm = "neldermead", opts = list(), hessian = T, ...) {
  check_nloptr <- require(nloptr)
  if (!check_nloptr) {
    stop("nloptr package must be installed to use package = \"nloptr\" in modelfitr")
  }

  all_algos <- matrix(c(
    "neldermead", "NLOPT_LN_NELDERMEAD",
    "lbfgs", "NLOPT_LD_LBFGS",
    "bobyqa", "NLOPT_LN_BOBYQA",
    "ccsaq", "NLOPT_LD_CCSAQ",
    "cobyla", "NLOPT_LN_COBYLA",
    "crs2lm", "NLOPT_GN_CRS2_LM",
    "direct", "NLOPT_GN_DIRECT",
    "direct_orig", "NLOPT_GN_ORIG_DIRECT",
    "direct_noscal", "NLOPT_GN_DIRECT_NOSCAL",
    "isres", "NLOPT_GN_ISRES",
    "mma", "NLOPT_LD_MMA",
    "newuoa", "NLOPT_LN_NEWUOA",
    "sbplx", "NLOPT_LN_SBPLX",
    "slsqp", "NLOPT_LD_SLSQP",
    "stogo", "NLOPT_GD_STOGO",
    "stogo_rand", "NLOPT_GD_STOGO_RAND",
    "tnewton", "NLOPT_LD_TNEWTON",
    "tnewton_restart", "NLOPT_LD_TNEWTON_RESTART",
    "varmetric", "NLOPT_LD_VAR1",
    "varmetric_rank2", "NLOPT_LD_VAR2"
  ),
  ncol = 2, byrow = T
  )

  if (substring(algorithm, 1, 5) != "NLOPT") {
    if (!(algorithm %in% all_algos[, 1])) {
      stop(paste("nloptr algorithm =", algorithm, "is not supported"))
    } else {
      algo_index <- which(all_algos[, 1] == algorithm)
      algorithm <- all_algos[algo_index, 2]
    }
  }

  opts <- c(algorithm = algorithm, opts)

  if ("lower" %in% names(list(...))) {
    list(...)$lb <- list(...)$lower
    list(...)$lower <- NULL
  }
  if ("upper" %in% names(list(...))) {
    list(...)$ub <- list(...)$upper
    list(...)$upper <- NULL
  }

  fit <- nloptr(start, objective, opts = opts, ...)

  fit$pars <- fit$solution
  fit$solution <- NULL
  fit$value <- fit$objective
  fit$objective <- NULL
  if (hessian) {
    fit$hess <- numDeriv::hessian(objective, fit$pars, ...)
  } else {
    fit$hess <- NULL
  }
  fit$code <- fit$status
  fit$status <- NULL

  return(fit)
}
