#' fit_optimx
#'
#' runs optimization using optimx package
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param max_runs integer; maximum number of optimx runs
#' @param conescutive integer; stop after this number of consecutive runs without improvement
#' @param sigma numeric; standard deviation for random perturbation of starting values
#' @param tolerance numeric; new fit must change objective by 'tolerance' to be considered an improvement
#' @param verbose logical; if TRUE, print result summary after each iteration
#' @param ... further arguments passed to run_optimx, optimx, and objective
#'
#' @export
fit_optimx <- function(objective, start, max_runs = 5L, consecutive = 0L, sigma = .1, tol = 1e-5, verbose = T, ...) {
  convergence <- F
  cons <- 0
  it <- 0

  fit <- run_optimx(objective, start, ...)
  if (!is.null(fit$hess)) {
    if (matrixcalc::is.positive.definite(fit$hess)) {
      convergence <- T
    }
  }

  while ((cons < consecutive) & (!convergence)) {
    new_fit <- run_optimx(objective, start, ...)

    if (!is.null(new_fit$hess)) {
      if (matrixcalc::is.positive.definite(new_fit$hess)) {
        fit <- new_fit
        convergence <- T
      }
    }

    if (!convergence) {
      if (new_fit$value < fit$value - tol) {
        fit <- new_fit
        start <- new_fit$pars
        cons <- 0
      } else {
        cons <- cons + 1
        start <- rnorm(length(start), start, abs(start) * sigma)
      }
    }

    if (verbose) {
      cat(sprintf("iteration = %d // cons = %d // obj = %0.3f // convergence = %s\n", it, cons, new_fit$value, convergence))
    }

    it <- it + 1
  }

  if (verbose) {
    cat(sprintf("FINAL :: iteration = %d // cons = %d // obj = %0.3f // convergence = %s\n", it, cons, fit$value, convergence))
  }

  return(fit)
}
