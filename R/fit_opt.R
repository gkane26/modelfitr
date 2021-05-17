#' fit_opt
#'
#' nonlinear optimization: wrapper for optimx and nloptr
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; starting parameters
#' @param lower numeric vector; lower bounds for parameters
#' @param upper numeric vector; upper bounds for parameters
#' @param hessian logical; if TRUE, calculate hessian at solution
#' @param package string; package to use, see Details. Default = "optimx"
#' @param method string; method to use from package
#' @param restart logical; if TRUE, restart optimization until no more improvement. Default = FALSE
#' @param conescutive integer; stop after this number of consecutive runs without improvement. Default = 0
#' @param max_runs integer; maximum number of optimx runs. Default = Inf
#' @param sigma numeric; standard deviation for random perturbation of starting values
#' @param tolerance numeric; new fit must change objective by 'tolerance' to be considered an improvement
#' @param verbose logical; if TRUE, print result summary after each iteration.Default=FALSE
#' @param opt_args list; list of arguments to be passed to optimx
#' @param obj_args list; list of arguments to be passed to objective function
#' @param ... further arguments objective
#'
#' @details Package options = "optimx", "optim", "nloptr", "pracma"
#' @export
fit_opt <- function(objective,
                    start,
                    lower = NULL,
                    upper = NULL,
                    hessian = FALSE,
                    package = "optimx",
                    method = NULL,
                    restart = FALSE,
                    consecutive = 0L,
                    max_runs = Inf,
                    sigma = .1,
                    tol = 1e-5,
                    verbose = FALSE,
                    opt_args = list(),
                    obj_args = list(),
                    ...) {


  if (length(start) == 1) {
    warn_1 = (package != "optim")
    if (!is.null(method)) {
      warn_1 = (method != "Brent")
    }

    if (warn_1) {
      warning("For single parameter optimization, using package = \"optim\" & method = \"Brent\".")
    }

    package = "optim"
    method = "Brent"
  }

  if (package == "optim") {
    fit_fn <- run_optim
  } else if (package == "optimx") {
    fit_fn <- run_optimx
  } else if (package == "nloptr") {
    fit_fn <- run_nloptr
  } else if (package == "pracma") {
    fit_fn <- run_pracma
  } else {
    stop(paste("package =", package, "is not currently supported."))
  }

  if (is.null(method)) {
    if (package == "optim") {
      method = "Nelder-Mead"
    } else if (package == "optimx") {
      method <- "nmkb"
    } else if (package == "nloptr") {
      method <- "neldermead"
    } else if (package == "pracma") {
      method <- "fminsearch"
    }
    warning(sprintf("method not provided, using default method = %s", method))
  }

  par_names <- names(start)
  convergence <- F

  ### iteration 1

  it <- 1

  if (verbose) {
    cat(sprintf("\niteration = %d\n", it))
  }

  fit <- .opt_chain(fit_fn,
    objective,
    start,
    lower = lower,
    upper = upper,
    hessian = hessian,
    method = method,
    restart = restart,
    verbose = verbose,
    opt_args = opt_args,
    obj_args = obj_args,
    tol=tol,
    ...
  )
  convergence <- ifelse(is.na(fit$res$convergence), FALSE, fit$res$convergence)

  ### if first run didn't converge and consecutive runs > 0, continue...

  for (s in 1:length(sigma)) {
    cons <- 0

    while ((!convergence) & (cons < consecutive) & (it < max_runs)) {
      it <- it + 1

      if (verbose) {
        cat(sprintf("\niteration = %d\n", it))
      }

      # randomize starting values
      start <- rnorm(length(fit$res$pars), fit$res$pars, ifelse(fit$res$pars == 0, sigma[s], sigma[s] * abs(fit$res$pars)))
      names(start) <- names(fit$res$pars)

      # refit model
      new_fit <- .opt_chain(fit_fn,
        objective,
        start,
        lower = lower,
        upper = upper,
        hessian = hessian,
        method = method,
        restart = restart,
        verbose = verbose,
        opt_args = opt_args,
        obj_args = obj_args,
        tol=tol,
        ...
      )
      convergence <- ifelse(is.na(new_fit$res$convergence), FALSE, new_fit$res$convergence)

      if (new_fit$res$value < fit$res$value - tol) {
        fit <- new_fit
        cons <- 0
      } else {
        cons <- cons + 1
      }

      if (verbose) {
        cat(sprintf("best obj = %0.3f // sigma = %0.3f // iterations without improvement = %d\n", fit$res$value, sigma[s], cons))
      }
    }
  }

  if ((verbose) & (it > 1)) {
    cat(sprintf("\nFINAL :: obj = %0.3f // convergence = %s // cons = %d\n", fit$res$value, fit$res$convergence, cons))
  }

  return(fit)
}

.opt_chain <- function(fit_fn,
                       objective,
                       start,
                       lower = NULL,
                       upper = NULL,
                       hessian = FALSE,
                       method = NULL,
                       restart = FALSE,
                       verbose = FALSE,
                       opt_args = list(),
                       obj_args = list(),
                       tol=1e-5,
                       ...) {
  it <- 1
  target <- do.call(objective, c(
    list(start),
    obj_args,
    list(...)
  ))

  if (verbose) {
    cat(sprintf("START :: obj = %0.3f\n", target))
  }

  fit <- do.call(fit_fn, c(
    objective = objective,
    list(start = start),
    list(lower = lower),
    list(upper = upper),
    hessian = hessian,
    method = method,
    opt_args,
    obj_args,
    list(...)
  ))
  res <- fit$res
  convergence <- ifelse(is.na(res$convergence), FALSE, res$convergence)

  if (verbose & !convergence & restart) {
    cat(sprintf("%d :: obj = %0.3f // code = %d // convergence = %s\n", it, res$value, res$code, res$convergence))
  }

  while ((restart) & (!convergence) & (res$value < (target - tol))) {
    it <- it + 1
    target <- res$value

    fit <- do.call(fit_fn, c(
      objective = objective,
      list(start = res$pars),
      list(lower = lower),
      list(upper = upper),
      hessian = hessian,
      method = method,
      opt_args,
      obj_args,
      list(...)
    ))
    res <- fit$res
    convergence <- ifelse(is.na(res$convergence), FALSE, res$convergence)

    if (verbose) {
      cat(sprintf("%d :: obj = %0.3f // code = %d // convergence = %s\n", it, res$value, res$code, res$convergence))
    }
  }

  if (verbose) {
    cat(sprintf("END :: obj = %0.3f // code = %d // convergence = %s\n", res$value, res$code, res$convergence))
  }


  return(fit)
}
