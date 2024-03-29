% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator_gee2.R
\name{get_params_gee}
\alias{get_params_gee}
\title{Second order generalized estimating equations (GEE2) estimation of CoCoA model parameters}
\usage{
get_params_gee(
  dat_xyz,
  jack_se = FALSE,
  init_mle = TRUE,
  no_meanmod = FALSE,
  no_meanscalemod = FALSE,
  rescale_xy = FALSE,
  rescale_zt = FALSE,
  cor_link = "tanh",
  sca_link = "log"
)
}
\arguments{
\item{dat_xyz}{a data frame containing the columns X, Y, Z, T, ... in that order}

\item{jack_se}{a logical value indicating if the jackknife standard errors should be computed}

\item{init_mle}{a logical value indicating if the parameter estimates should be initialized at the MLE}

\item{no_meanmod}{a logical value indicating if the mean model should be estimated (see Details)}

\item{no_meanscalemod}{a logical value indicating if the mean and scale models should be estimated (see Details)}

\item{rescale_xy}{a logical value indicating if X and Y should be rescaled (z-scored)}

\item{rescale_zt}{a logical value indicating if the conditional correlation model covariates Z, T, ... should be rescaled (z-scored)}

\item{cor_link}{a character string indicating the correlation model link function; passed on to \code{geese()}}

\item{sca_link}{a character string indicating the scale model link function; passed on to \code{geese()}}
}
\value{
A list containing \code{params}, a numeric matrix of the estimated parameters, \code{out2}, a table of the estimated parameters,
\code{runtime}, the runtime of the optimization routine.
}
\description{
Estimate CoCoA model parameters using second order generalized estimating equations (GEE2).
}
\details{
Given \eqn{n} observations of Z and T, the conditional mean, scale, and correlation of X and Y
is described using second order generalized estimating equations (GEE2), introduced
in Yan and Fine (2004). In particular, the GEE2 model includes the correlation model
\deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
This function is a wrapper around the \code{geese()} function from geepack, and is an
expanded version of \code{CNM.full()} from the LiquidAssociation package.

If \code{no_meanmod = TRUE}, then only an intercept term will be estimated in the mean model (i.e.,
all participants will have the same predicted mean of X and the same predicted mean of Y).
If \code{no_meanscalemod = TRUE}, then only an intercept term will be estimated in both
the mean and scale models.
}
\examples{
T = as.numeric(purrr::rbernoulli(n = 100, p = 0.5))
Z = rep(NA, times = 100)
for (i in 1:100){
Z[i] <- if (T[i] == 1){
# Z|T=1 ~ Bernoulli(0.3)
as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
} else {
# Z|T=0 ~ Bernoulli(0.9)
as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
}
}
dat_zt0 <- data.frame(Z = Z, T = T, ZT = Z*T)
XY = sim_xy2_mvn(theta0 = c(0, 0, 0, 0), dat_zt0 = dat_zt0,
mu_x = 1, mu_y = 1, sd_x = 1, sd_y = 2, cor_link = "tanh")
dat_xyz = data.frame(X = XY$x, Y = XY$y)
dat_xyz = cbind(dat_xyz, dat_zt0)
get_params_gee(dat_xyz)

}
