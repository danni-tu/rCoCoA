% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator_reml.R
\name{loglik_reml}
\alias{loglik_reml}
\title{Restricted multivariate normal likelihood of CoCoA model parameters}
\usage{
loglik_reml(params, dat_xyz)
}
\arguments{
\item{params}{a numeric vector of the form \code{(mu_x, mu_y, sd_x, sd_y, alpha, beta, ...)} in that order}

\item{dat_xyz}{a data frame containing the columns X, Y, Z, T, ... in that order}
}
\value{
The restricted log-likelihood value multiplied by -1.
}
\description{
Restricted multivariate normal likelihood of CoCoA model parameters
}
\details{
Given \eqn{n} observations of Z and T, X and Y are bivariate normal with
means \code{mu_x, mu_y}, standard deviations \code{sd_x, sd_y}, and correlation
\deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
Currently, only the \code{tanh} link is implemented. Apart from the X and Y columns, \code{dat_xyz} may include
multiple covariates Z, T, ... for the conditional correlation model. The number of these covariates should be
equal to the length of \code{params} minus 4.
}
