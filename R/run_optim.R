#' run_optim
#'
#' runs optimization using optim package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds for parameters
#' @param upper numeric vector; upper bounds for parameters
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param method string; method to use. See optim for details
#' @param control list; a list of control parameters, see optim for details
#' @param ... further arguments passed to optim and objective
#'
#' @export
run_optim <- function(objective,
                      start,
                      lower = -Inf,
                      upper = Inf,
                      hessian = FALSE,
                      method = "Nelder-Mead",
                      control = list(),
                      ...) {

  if (any(is.null(lower))) lower <- -Inf
  if (any(is.null(upper))) upper <- Inf

  fit <- optim(
    start,
    objective,
    lower = lower,
    upper = upper,
    method = method,
    control = control,
    ...
  )

  fit_pars <- fit$par
  names(fit_pars) <- names(start)
  fit_val <- fit$value
  if (hessian) {
    fit_hess <- numDeriv::hessian(objective, fit_pars, ...)
    fit_conv <- matrixcalc::is.positive.definite(fit_hess)
  } else {
    fit_hess <- NA
    fit_conv <- NA
  }
  fit_code <- fit$convergence

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess,
    convergence = fit_conv,
    code = fit_code
  )

  return(list(res = res, fit = fit))
}
