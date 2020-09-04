#' run_optimx
#'
#' runs optimization using optimx package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds for parameters
#' @param upper numeric vector; upper bounds for parameters
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param method string; method to use. See optimx for details
#' @param control list; a list of control parameters, see optimx for details
#' @param ... further arguments passed to optimx and objective
#'
#' @export
run_optimx <- function(objective,
                       start,
                       lower = -Inf,
                       upper = Inf,
                       hessian = FALSE,
                       method = "nmkb",
                       control = list(),
                       ...) {
  if (!("starttests" %in% names(control))) control <- c(control, starttests = F)
  if (!("kkt" %in% names(control))) control <- c(control, kkt = F)

  if (any(is.null(lower))) lower <- -Inf
  if (any(is.null(upper))) upper <- Inf

  fit <- optimx::optimx(
    start,
    objective,
    lower = lower,
    upper = upper,
    method = method,
    control = control,
    ...
  )

  fit_pars <- as.numeric(fit[1, 1:length(start)])
  names(fit_pars) <- names(fit)[1:length(start)]
  fit_val <- fit[1, "value"]
  fit_hess <- ifelse(hessian, numDeriv::hessian(objective, fit_pars, ...), NA)
  fit_conv <- ifelse(is.na(fit_hess), NA, matrixcalc::is.positive.definite(fit_hess))
  fit_code <- fit[1, "convcode"]

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess,
    convergence = fit_conv,
    code = fit_code
  )

  return(list(res = res, fit = fit))
}
