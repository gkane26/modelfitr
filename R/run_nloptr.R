#' run_nloptr
#'
#' runs optimization using nloptr package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param method string; method to use. See details. Default = "neldermead"
#' @param control list; list of options. See nloptr::nloptr.print.options() for details.
#' @param hessian logical; if TRUE, find the hessian at the optimum. Default = T
#' @param ... further arguments passed to nloptr method and objective
#'
#' @details
#' method options =
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
run_nloptr <- function(objective,
                       start,
                       lower = NULL,
                       upper = NULL,
                       hessian = FALSE,
                       method = "neldermead",
                       control = list(),
                       ...) {
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

  if (substring(method, 1, 5) != "NLOPT") {
    if (!(method %in% all_algos[, 1])) {
      stop(paste("nloptr method =", method, "is not supported"))
    } else {
      algo_index <- which(all_algos[, 1] == method)
      method <- all_algos[algo_index, 2]
    }
  }

  control <- c(algorithm = method, control)

  fit <- nloptr(start, objective, lb = lower, ub = upper, opts = control, ...)

  fit_pars <- fit$solution
  names(fit_pars) <- names(start)
  fit_val <- fit$objective
  fit_hess <- ifelse(hessian, numDeriv::hessian(objective, fit_pars, ...), NA)
  fit_conv <- ifelse(is.na(fit_hess), NA, matrixcalc::is.positive.definite(fit_hess))
  fit_code <- fit$status

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess,
    convergence = fit_conv,
    code = fit_code
  )

  return(list(res = res, fit = fit))
}
