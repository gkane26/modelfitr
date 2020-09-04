#' fit_gensa
#'
#' runs optimization using GenSA package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds on parameters
#' @param upper numeric vector; upper bounds on parameters
#' @param hessian logical; if TRUE, calculate hessian at solution. Default = FALSE
#' @param opt_args list; list of arguments passed on to GenSA, see GenSA function for details
#' @param obj_args list; list of arguments to pass on to objective function
#' @param ... further arguments passed to objective
#'
#' @export
fit_gensa <- function(objective,
                      start,
                      lower,
                      upper,
                      hessian = FALSE,
                      opt_args = list(),
                      obj_args = list(),
                      ...) {
  fit <- do.call(GenSA::GenSA, c(
    list(start),
    objective,
    list(lower = lower),
    list(upper = upper),
    opt_args,
    obj_args,
    list(...)
  ))

  fit_pars <- fit$par
  names(fit_pars) <- names(start)
  fit_val <- fit$value
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
