#' fit_rgenoud
#'
#' runs optimization using GenSA package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds on parameters
#' @param upper numeric vector; upper bounds on parameters
#' @param hessian logical; if TRUE, calculate hessian at solution. Default = FALSE
#' @param opt_args list; list of arguments to pass to genoud. see rgenoud::genoud for details
#' @param obj_args list; list of arguments to pass on to objective function
#' @param ... further arguments passed to objective
#'
#' @export
fit_rgenoud <- function(objective,
                        start = NULL,
                        lower = NULL,
                        upper = NULL,
                        hessian = FALSE,
                        opt_args = list(),
                        obj_args = list(),
                        ...) {
  if (!is.null(lower) & !is.null(upper)) {
    domains <- cbind(lower, upper)
    boundary.enforcement <- 2
  } else {
    domains <- NULL
    boundary.enforcement <- min(boundary.enforcement, 1)
  }

  if (!("Domains" %in% opt_args)) opt_args$Domains <- domains
  if (!("boundary.enforcement" %in% opt_args)) opt_args$boundary.enforcement <- boundary.enforcement

  if (!is.null(start)) {
    opt_args$nvars <- length(start)
  } else {
    if (!("nvars" %in% opt_args)) {
      stop("Must supply either a starting vector or n_vars for package = rgenoud.")
    }
  }

  opt_args$data.type.int <- ifelse(is.null(opt_args$data.type.int), FALSE, opt_args$data.type.int)

  fit <- do.call(rgenoud::genoud, c(
    objective,
    opt_args,
    obj_args,
    list(...)
  ))

  fit_pars <- fit$par
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
