#' run_optimx
#'
#' runs optimization using optimx package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param method string; method to use. See optimx for details
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param control list; a list of control parameters, see optimx for details
#' @param ... further arguments passed to optimx and objective
#'
#' @export
run_optimx <- function(objective, start, method = "nmkb", hessian = T, control = list(), failcode = 9999, ...) {
  if (!("starttests" %in% names(control))) control <- c(control, starttests = F)
  if (!("kkt" %in% names(control))) control <- c(control, kkt = F)

  fit <- optimx::optimx(start, objective, method = method, hessian = hessian, control = control, ...)

  fit_pars <- as.numeric(fit[1, 1:length(start)])
  names(fit_pars) <- names(fit)[1:length(start)]

  fit_val <- fit[1, "value"]
  fit_hess <- attr(fit, "details")$nhatend
  fit_code <- fit[1, "convcode"]

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess,
    code = fit_code
  )

  return(res)
}
