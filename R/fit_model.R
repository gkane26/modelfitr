#' fit_model
#'
#' run optimization algorithms from a variety of packages
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; parameter vector. Default = NULL
#' @param lower numeric vector; lower bounds on parameters. Default = NULL
#' @param upper numeric vector; upper bounds on parameters. Default = NULL
#' @param hessian logical; if TRUE, calculate hessian at final solution; Default = FALSE.
#' @param package string; package to use for optimization, see Details. Default = "optimx"
#' @param method string; method to use from the package (if applicable, see Details for options). Default = NULL (use method default for package)
#' @param restart logical; if TRUE, restart optimization until no more improvement. Only applicable for package = optimx, nloptr, & pracma. Default = FALSE
#' @param conescutive integer; stop after this number of consecutive runs without improvement. Only applicable for package = optimx, nloptr, & pracma. Default = 0
#' @param max_runs integer; maximum number of optimx runs. Only applicable for package = optimx, nloptr, & pracma. Default = Inf
#' @param sigma numeric vector; standard deviation to perturb starting values. For packages = optimx, nloptr, pracma, if consecutive > 0, will perturb start values before each successive run. For packages = DEoptim, RcppDE, GA, used to create random initial population (with specified starting values). Default = 0.1
#' @param n_pop numeric vector; size of population. Only applicable for package = RcppDE, DEoptim, GA. Default = 50, recommended to be 10 * n parameters.
#' @param opt_args list; further arguments passed to fitting method
#' @param obj_args list; further arguments passed to objective function
#' @param aic logical; if TRUE, calculate Akaike Information Criterion. Default = FALSE
#' @param bic logical; if TRUE, calculate Bayes Information Criterion. Default = FALSE
#' @param n_obs integer; the number of observations (used to calculate BIC). Default = NULL
#' @param return_df logical; if TRUE, return results as a data frame
#' @param return_all logical; if TRUE, return list of results from modelfitr AND results from the package used to fit the model
#' @param ... further arguments passed to objective function (similar to obj_args)
#'
#' @details
#' Still under developement... \cr
#'
#' Package options = "optimx", "nloptr", "pracma", GenSA", "rgenoud", "GA", "DEoptim", "RcppDE". \cr
#' Package = "optimx", "nloptr", "pracma" should work. \cr
#' Others are implemented but no guarantees. \cr
#' Further documentation coming soon. \cr
#'
#' @return
#'
#' if return_df = FALSE and return_all = FALSE (default),
#' returns a list containing the following fields:
#'
#' \item{pars}{parameter values}
#' \item{value}{the value of the objective function at the solution}
#' \item{hess}{if hessian=TRUE, the hessian matrix at the solution, otherwise NA}
#' \item{convergence}{if hessian=TRUE, convergence = TRUE if the hessian is positive definite, otherwise FALSE. if hessian = FALSE, convergence = NA}
#' \item{code}{the message or code returned by the optimization method if applicable, otherwise NA}
#'
#' if aic = TRUE and/or bic = TRUE, list also contains:
#' \item{aic}{the aic for the model}
#' \item{bic}{the bic for the model}
#'
#' if return_df = TRUE, returns a data frame instead of a list. This data frame has 1 row and columns for:
#' \itemize{
#' \item each parameter
#' \item value
#' \item covergence
#' \item code
#' \item aic - if applicable
#' \item bic - if applicable
#' }
#'
#' if return_all = TRUE, returns a list containing:
#' \itemize{
#' \item the list or data frame as described above
#' \item the direct output from the optimization method
#' }
#'
#'
#' @export
fit_model <- function(objective,
                      start = NULL,
                      lower = NULL,
                      upper = NULL,
                      hessian = FALSE,
                      package = "optimx",
                      method = NULL,
                      restart = FALSE,
                      consecutive = 0L,
                      max_runs = Inf,
                      sigma = 0.1,
                      n_pop = 50,
                      opt_args = list(),
                      obj_args = list(),
                      aic = FALSE,
                      bic = FALSE,
                      n_obs = NULL,
                      return_df = FALSE,
                      return_all = FALSE,
                      ...) {
  if (package %in% c("optimx", "nloptr", "pracma")) {
    if (is.null(start)) {
      stop(paste("Must provide starting values with package =", package))
    }

    fit <- fit_opt(
      objective,
      start,
      lower = lower,
      upper = upper,
      hessian = hessian,
      package = package,
      method = method,
      restart = restart,
      consecutive = consecutive,
      max_runs = max_runs,
      sigma = sigma,
      opt_args = opt_args,
      obj_args = obj_args,
      ...
    )
  } else if (package == "GenSA") {
    if (is.null(lower) | is.null(upper)) {
      stop("Lower and upper bounds are required for package = GenSA")
    }

    fit <- fit_gensa(
      objective,
      start,
      lower,
      upper,
      hessian = hessian,
      opt_args = opt_args,
      obj_args = obj_args,
      ...
    )
  } else if (package == "rgenoud") {
    fit <- fit_rgenoud(
      objective,
      start,
      lower,
      upper,
      hessian = hessian,
      opt_args = opt_args,
      obj_args = obj_args,
      ...
    )
  } else if (package == "GA") {
    fit <- fit_ga(
      objective,
      start = start,
      lower = lower,
      upper = upper,
      hessian = hessian,
      sigma = sigma,
      n_pop = n_pop,
      opt_args = opt_args,
      obj_args = obj_args,
      ...
    )
  } else if (package %in% c("DEoptim", "RcppDE")) {
    fit <- fit_de(
      objective,
      start = start,
      lower = lower,
      upper = upper,
      hessian = hessian,
      sigma = sigma,
      n_pop = n_pop,
      opt_args = opt_args,
      obj_args = obj_args,
      ...
    )
  } else {
    stop(paste("package =", package, "is not implemented!"))
  }

  if (aic) fit$res$aic <- 2 * length(fit$res$pars) - 2 * log(fit$res$value)
  if (bic) {
    if (is.null(n_obs)) {
      warning("Not calculating BIC! Must provide number of observations to calculate BIC.")
    } else {
      fit$res$bic <- length(fit$res$pars) * log(n_obs) - 2 * log(fit$res$value)
    }
  }

  if (return_df) {
    fit_data <- data.frame(as.list(fit$res$pars), value = fit$res$value, convergence = fit$res$convergence)
    if ("aic" %in% names(fit$res)) fit_data$aic <- fit$res$aic
    if ("bic" %in% names(fit$res)) fit_data$bic <- fit$res$bic
    fit$res <- fit_data
  }

  if (return_all) {
    return(fit)
  } else {
    return(fit$res)
  }
}
