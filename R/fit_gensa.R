#' fit_gensa
#'
#' runs optimization using GenSA package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds on parameters
#' @param upper numeric vector; upper bounds on parameters
#' @param control list; list of control parameters, see GenSA function for details
#' @param ... further arguments passed to objective
#'
#' @export
fit_gensa <- function(objective, start, lower, upper, control = list(), hessian = T, ...) {
  fit <- GenSA::GenSA(start,
    objective,
    lower = lower,
    upper = upper,
    control = control,
    ...
  )

  fit$pars <- fit$par
  fit$par <- NULL

  if (hessian) {
    fit$hess <- numDeriv::hessian(objective, fit$pars, ...)
  }

  return(fit)
}
