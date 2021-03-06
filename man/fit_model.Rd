% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model.R
\name{fit_model}
\alias{fit_model}
\title{fit_model}
\usage{
fit_model(
  objective,
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
  verbose = FALSE,
  n_pop = 50,
  opt_args = list(),
  obj_args = list(),
  aic = FALSE,
  bic = FALSE,
  n_obs = NULL,
  return_df = FALSE,
  return_all = FALSE,
  ...
)
}
\arguments{
\item{objective}{function; the objective function to minimize}

\item{start}{numeric vector; parameter vector. Default = NULL}

\item{lower}{numeric vector; lower bounds on parameters. Default = NULL}

\item{upper}{numeric vector; upper bounds on parameters. Default = NULL}

\item{hessian}{logical; if TRUE, calculate hessian at final solution; Default = FALSE.}

\item{package}{string; package to use for optimization, see Details. Default = "optimx"}

\item{method}{string; method to use from the package (if applicable, see Details for options). Default = NULL (use method default for package)}

\item{restart}{logical; if TRUE, restart optimization until no more improvement. Only applicable for package = optimx, nloptr, & pracma. Default = FALSE}

\item{max_runs}{integer; maximum number of optimx runs. Only applicable for package = optimx, nloptr, & pracma. Default = Inf}

\item{sigma}{numeric vector; standard deviation to perturb starting values. For packages = optimx, nloptr, pracma, if consecutive > 0, will perturb start values before each successive run. For packages = DEoptim, RcppDE, GA, used to create random initial population (with specified starting values). Default = 0.1}

\item{verbose}{logical; if TRUE, print progress statements. Only applicable for packages = optimx, nloptrm, pracma. Default=FALSE}

\item{n_pop}{numeric vector; size of population. Only applicable for package = RcppDE, DEoptim, GA. Default = 50, recommended to be 10 * n parameters.}

\item{opt_args}{list; further arguments passed to fitting method}

\item{obj_args}{list; further arguments passed to objective function}

\item{aic}{logical; if TRUE, calculate Akaike Information Criterion. Default = FALSE}

\item{bic}{logical; if TRUE, calculate Bayes Information Criterion. Default = FALSE}

\item{n_obs}{integer; the number of observations (used to calculate BIC). Default = NULL}

\item{return_df}{logical; if TRUE, return results as a data frame}

\item{return_all}{logical; if TRUE, return list of results from modelfitr AND results from the package used to fit the model}

\item{...}{further arguments passed to objective function (similar to obj_args)}

\item{conescutive}{integer; stop after this number of consecutive runs without improvement. Only applicable for package = optimx, nloptr, & pracma. Default = 0}
}
\value{
if return_df = FALSE and return_all = FALSE (default),
returns a list containing the following fields:

\item{pars}{parameter values}
\item{value}{the value of the objective function at the solution}
\item{hess}{if hessian=TRUE, the hessian matrix at the solution, otherwise NA}
\item{convergence}{if hessian=TRUE, convergence = TRUE if the hessian is positive definite, otherwise FALSE. if hessian = FALSE, convergence = NA}
\item{code}{the message or code returned by the optimization method if applicable, otherwise NA}

if aic = TRUE and/or bic = TRUE, list also contains:
\item{aic}{the aic for the model}
\item{bic}{the bic for the model}

if return_df = TRUE, returns a data frame instead of a list. This data frame has 1 row and columns for:
\itemize{
\item each parameter
\item value
\item covergence
\item code
\item aic - if applicable
\item bic - if applicable
}

if return_all = TRUE, returns a list containing:
\itemize{
\item the list or data frame as described above
\item the direct output from the optimization method
}
}
\description{
run optimization algorithms from a variety of packages
}
\details{
Still under developement... \cr

Package options = "optimx", "nloptr", "pracma", "GenSA", "rgenoud", "GA", "DEoptim", "RcppDE". \cr
Others are implemented but no guarantees. \cr
Further documentation coming soon. \cr
}
