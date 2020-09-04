#' run_pracma
#'
#' runs optimization using optimx package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param method string; method to use. See details. Default = "fminsearch".
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param ... further arguments passed to optimx and objective
#'
#' @details method options = "fminsearch", "fminunc", "fmincon", "fminbnd"
#'
#' @export
run_pracma <- function(objective,
                       start,
                       lower = NULL,
                       upper = NULL,
                       hessian = FALSE,
                       method = "fminsearch",
                       ...) {
  check_pracma <- require(pracma)
  if (!check_pracma) {
    stop("pracma package must be installed to use this method.")
  }

  if ("lower" %in% names(list(...))) {
    list(...)$lb <- list(...)$lower
    list(...)$lower <- NULL
  }
  if ("upper" %in% names(list(...))) {
    list(...)$ub <- list(...)$upper
    list(...)$upper <- NULL
  }

  if (method == "fminsearch") {
    fit <- pracma::fminsearch(objective, start, lower = lower, upper = upper, ...)
  } else if (method == "fmincon") {
    fit <- pracma::fmincon(start, objective, lb = lower, ub = upper, ...)
  } else if (method == "fminunc") {
    if (!is.null(lower) | !is.null(upper)) {
      warning("pracma::fminunc does not use bounds. Ignoring lower and upper.")
    }
    fit <- pracma::fminunc(start, objective, ...)
  } else {
    stop(paste("method =", method, "is not supported"))
  }

  fit_pars <- ifelse(method == "fminsearch", fit$xopt, fit$par)
  names(fit_pars) <- names(start)
  fit_val <- ifelse(method == "fminsearch", fit$min, fit$value)
  fit_hess <- ifelse(hessian, numDeriv::hessian(objective, fit_pars, ...), NA)
  fit_conv <- ifelse(is.na(fit_hess), NA, matrixcalc::is.positive.definite(fit_hess))
  fit_code <- fit$convergence

  return(list(res = res, fit = fit))
}
