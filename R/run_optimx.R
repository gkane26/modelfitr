#' run_optimx
#'
#' runs optimization using optimx package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param control list; a list of control parameters, see optimx for details
#' @param ... further arguments passed to optimx and objective
#'
#' @export
run_optimx <- function(objective, start, hessian = T, control = list(), ...) {
  if (!("starttests" %in% names(control))) control <- c(control, starttests = F)
  if (!("kkt" %in% names(control))) control <- c(control, kkt = F)

  fit <- optimx::optimx(start, objective, hessian = hessian, control = control, ...)

  fit_pars <- as.numeric(fit[1, 1:length(start)])
  names(fit_pars) <- names(fit)[1:length(start)]

  fit_val <- fit[1, "value"]
  fit_hess <- attr(fit, "details")$nhatend

  res <- list(
    pars = fit_pars,
    value = fit_val,
    hess = fit_hess
  )

  return(res)
}
