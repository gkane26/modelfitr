#' fit_ga
#'
#' runs optimization using GA package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds on parameters
#' @param upper numeric vector; upper bounds on parameters
#' @param hessian logical; if TRUE, find the hessian at the optimum
#' @param sigma numeric; standard deviation to create random initial population
#' @param n_pop integer; population size
#' @param minimize logical; if TRUE, maximize -1 * objective
#' @param opt_args list; list of arguments to pass to ga, see GA::ga for details
#' @param obj_args list; list of arguments to pass on to objective function
#' @param ... further arguments passed to objective
#'
#' @export
fit_ga <- function(objective,
                   start = NULL,
                   lower = NULL,
                   upper = NULL,
                   hessian = T,
                   sigma = 0.1,
                   n_pop = 50L,
                   minimize = T,
                   opt_args = list(),
                   obj_args = list(),
                   ...) {
  if (is.null(lower) | is.null(upper)) {
    stop("Bounds are required for GA")
  }

  if (minimize) {
    obj <- function(...) -objective(...)
  } else {
    obj <- objective
  }

  opt_args$type <- ifelse(is.null(opt_args$type), "real-valued", opt_args$type)

  n_pop <- ifelse("popSize" %in% names(opt_args), opt_args$popSize, n_pop)
  if (!("suggestions" %in% names(opt_args))) {
    if (!(is.null(start))) {
      opt_args$popSize <- n_pop
      opt_args$suggestions <- sapply(start, function(x) rnorm(n_pop, x, sigma * abs(x)))
    }
  }

  fit <- do.call(GA::ga, c(
    fitness = obj,
    list(lower = lower),
    list(upper = upper),
    opt_args,
    obj_args,
    list(...)
  ))

  fit_pars <- as.numeric(fit@solution)
  names(fit_pars) <- names(start)
  fit_val <- ifelse(minimize, -fit@fitnessValue, fit@fitnessValue)
  if (hessian) {
    fit_hess <- numDeriv::hessian(objective, fit_pars, ...)
    fit_conv <- matrixcalc::is.positive.definite(fit_hess)
  } else {
    fit_hess <- NA
    fit_conv <- NA
  }
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
