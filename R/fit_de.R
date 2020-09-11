#' fit_de
#'
#' runs optimization using RcppDE package (Rcpp version of DEoptim)
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds on parameters
#' @param upper numeric vector; upper bounds on parameters
#' @param hessian logical; if TRUE, hessian is calculated at solution. Default = FALSE
#' @param sigma numeric vector; standard deviation to create random initial population
#' @param n_pop numeric vector; size of population
#' @param opt_args list; list of arguments passed to DEoptim, see RcppDE::DEoptim for details
#' @param obj_args list; list of arguments to pass on to objective function
#' @param ... further arguments passed to objective
#'
#' @export
fit_de <- function(objective,
                   start = NULL,
                   lower = NULL,
                   upper = NULL,
                   hessian = FALSE,
                   sigma = 0.1,
                   n_pop = 50L,
                   opt_args = list(),
                   obj_args = list(),
                   ...) {
  if (is.null(lower) | is.null(upper)) {
    stop("Bounds must be specified for DEoptim/RcppDE")
  }

  pass_through <- list(...)
  if ("control" %in% names(pass_through)) {
    if (!("control" %in% names(opt_args))) {
      opt_args$control <- pass_through$control
    } else {
      warning("control specified in opt_args and as a pass through argument -- ignoring the pass through control.")
    }
    pass_through$control <- NULL
  }

  if (!is.null(start)) {
    if ("control" %in% names(opt_args)) {
      if (!("initialpop" %in% names(opt_args$control))) {
        n_pop <- ifelse("NP" %in% names(opt_args$control), opt_args$control$NP, n_pop)
        opt_args$control$NP <- n_pop
        opt_args$control$initialpop <- sapply(start, function(x) rnorm(n_pop, x, sigma * abs(x)))
      }
    } else {
      opt_args$control <- list()
      opt_args$control$NP <- n_pop
      opt_args$control$initialpop <- sapply(start, function(x) rnorm(n_pop, x, sigma * abs(x)))
    }
  }

  fit <- do.call(RcppDE::DEoptim, c(
    objective,
    list(lower = lower),
    list(upper = upper),
    opt_args,
    obj_args,
    pass_through
  ))

  fit_pars <- fit$optim$bestmem
  fit_val <- fit$optim$bestval
  fit_hess <- ifelse(hessian, numDeriv::hessian(objective, fit_pars, ...), NA)
  fit_conv <- ifelse(is.na(fit_hess), NA, matrixcalc::is.positive.definite(fit_hess))
  fit_code <- NA

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess,
    convergence = fit_conv,
    code = fit_code
  )

  return(list(fit = fit, res = res))
}
