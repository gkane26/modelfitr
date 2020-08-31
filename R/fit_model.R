#' fit_model
#'
#' run optimization algorithm on given objective function
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; parameter vector
#' @param package string; package to use for optimization, see Details. Default = "optimx".
#' @param method string; fitting method. Default = Nelder-Mead
#' @param aic logical; if TRUE, calculate Akaike Information Criterion. Default = FALSE
#' @param bic logical; if TRUE, calculate Bayes Information Criterion. Default = FALSE
#' @param n_obs integer; the number of observations (used to calculate BIC)
#' @param ... further arguments passed to fitting method  and objective function and objective function
#'
#' @details
#' Package options = "optimx", "nloptr", "GenSA"
#'
#' @export
fit_model <- function(objective,
                      start = NULL,
                      lower = NULL,
                      upper = NULL,
                      package = "optimx",
                      aic = F,
                      bic = F,
                      n_obs = NULL,
                      as_data_frame = F, ...) {
  if (package %in% c("optimx", "nloptr")) {
    if (is.null(start)) {
      stop(paste("Must provide starting values with package =", package))
    }

    if (package == "optimx") {
      if (is.null(lower)) lower <- -Inf
      if (is.null(upper)) upper <- Inf
    }

    fit <- opt(objective, start, lower = lower, upper = upper, package = package, ...)
  } else if (package == "GenSA") {
    if (is.null(lower) | is.null(upper)) {
      stop("Lower and upper bounds are required for GenSA")
    }

    fit <- fit_gensa(objective, start, lower, upper, ...)
  } else {
    stop(paste("package =", package, "is not implemented!"))
  }

  if (aic) fit$aic <- 2 * length(fit$pars) - 2 * log(fit$value)
  if (bic) {
    if (is.null(n_obs)) {
      warning("Not calculating BIC! Must provide number of observations to calculate BIC.")
    } else {
      fit$bic <- length(fit$pars) * log(n_obs) - 2 * log(fit$value)
    }
  }

  if (as_data_frame) {
    fit_data <- data.frame(as.list(fit$pars), value = fit$value)
    if ("hess" %in% names(fit)) {
      if (!is.null(fit$hess)) {
        fit_data$hess_is_pdm <- matrixcalc::is.positive.definite(fit$hess)
      }
    }
    if ("aic" %in% names(fit)) fit_data$aic <- fit$aic
    if ("bic" %in% names(fit)) fit_data$bic <- fit$bic
    return(fit_data)
  } else {
    return(fit)
  }
}
