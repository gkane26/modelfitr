#' fit_model
#'
#' run optimization algorithm on given objective function
#'
#' @param objective function; the objective function to minimize
#' @param start numeric vector; parameter vector
#' @param method string; fitting method; by default, Nelder-Mead; see Details
#' @param aic logical; if TRUE, calculate Akaike Information Criterion
#' @param bic logical; if TRUE, calculate Bayes Information Criterion
#' @param n_obs integer; the number of observations (used to calculate BIC)
#' @param ... further arguments passed to fit_optimx, run_optimx, optimx, and objective
#'
#' @export
fit_model <- function(objective, start = NULL, method = "Nelder-Mead", aic = T, bic = F, n_obs = NULL, as_data_frame = F, ...) {
  optimx_methods <- c(
    "Nelder-Mead",
    "BFGS",
    "CG",
    "L-BFGS-B",
    "nlm",
    "nlminb",
    "spg",
    "ucminf",
    "newuoa",
    "bobyqa",
    "nmkb",
    "hjkb",
    "Rcgmin",
    "Rvmmin"
  )

  if (method %in% optimx_methods) {
    if (is.null(start)) {
      stop(paste("Must provide starting values with method =", method))
    }

    fit <- fit_optimx(objective, start, method = method, ...)
  } else {
    stop(paste("method =", method, "is not implemented!"))
  }

  if (aic) fit$aic <- 2 * length(fit$pars) - 2 * log(fit$value)
  if (bic) {
    if (is.null(n_obs)) {
      warning("Not calculating BIC! Must provide number of observations to calculate BIC.")
    } else {
      fit$bic <- length(fit$pars) * log(n_obs) - 2 * log(fit_value)
    }
  }

  if (as_data_frame) {
    fit_data <- data.frame(as.list(fit$pars), value = fit$value)
    if ("aic" %in% names(fit)) fit_data$aic <- fit$aic
    if ("bic" %in% names(fit)) fit_data$bic <- fit$bic
    return(fit_data)
  } else {
    return(fit)
  }
}
