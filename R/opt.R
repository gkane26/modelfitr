#' opt
#'
#' nonlinear optimization: wrapper for optimx and nloptr
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param package string; package to use, see Details. Default = "optimx"
#' @param failcode integer vector; convergence code from optimizer to be considered failure
#' @param max_runs integer; maximum number of optimx runs
#' @param conescutive integer; stop after this number of consecutive runs without improvement
#' @param sigma numeric; standard deviation for random perturbation of starting values
#' @param tolerance numeric; new fit must change objective by 'tolerance' to be considered an improvement
#' @param verbose logical; if TRUE, print result summary after each iteration
#' @param ... further arguments passed to run_optimx, optimx, and objective
#'
#' @details Package options = "optimx", "nloptr"
#' @export
opt <- function(objective, start, package = "optimx", failcode = NULL, max_runs = 5L, consecutive = 0L, sigma = .1, tol = 1e-5, verbose = T, ...) {
  if (package == "optimx") {
    fit_fn <- run_optimx
  } else {
    fit_fn <- run_nloptr
  }

  par_names <- names(start)
  convergence <- F
  cons <- 0
  it <- 1

  fit <- fit_fn(objective, start, ...)
  start <- fit$pars

  if (!(fit$code %in% failcode)) {
    if (!is.null(fit$hess)) {
      if (matrixcalc::is.positive.definite(fit$hess)) {
        convergence <- T
      }
    }
  }

  if (verbose) {
    cat(sprintf("iteration = %d // obj = %0.3f // code = %d // convergence = %s\n", it, fit$value, fit$code, convergence))
  }


  while ((cons < consecutive) & (!convergence)) {
    it <- it + 1

    new_fit <- fit_fn(objective, start, ...)

    if (!(new_fit$code) %in% failcode) {
      if (!is.null(new_fit$hess)) {
        if (matrixcalc::is.positive.definite(new_fit$hess)) {
          fit <- new_fit
          convergence <- T
        }
      }
    }

    if (!convergence) {
      if (new_fit$value < fit$value - tol) {
        fit <- new_fit
        start <- new_fit$pars
        cons <- 0
      } else {
        if (!(new_fit$code %in% failcode)) {
          cons <- cons + 1
        }
        start <- rnorm(length(start), start, abs(start) * sigma)
        names(start) <- par_names
      }
    }

    if (verbose) {
      cat(sprintf("iteration = %d // cons = %d // obj = %0.3f // code = %d // convergence = %s\n", it, cons, new_fit$value, new_fit$code, convergence))
    }
  }

  if ((verbose) & (it > 1)) {
    cat(sprintf("FINAL :: iteration = %d // cons = %d // obj = %0.3f // convergence = %s\n", it, cons, fit$value, convergence))
  }

  return(fit)
}
