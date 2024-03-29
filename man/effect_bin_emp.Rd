% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/association_size.R
\name{effect_bin_emp}
\alias{effect_bin_emp}
\title{Association Size of Binary Variable T}
\usage{
effect_bin_emp(params, dat5)
}
\arguments{
\item{params}{a numeric vector containing the estimated correlation parameters in the CoCoA model}

\item{dat5}{a data frame containing only the covariates in the correlation model (i.e., excluding X and Y)}
}
\value{
A numeric vector containing the marginal and conditional association sizes.
}
\description{
Association Size of Binary Variable T
}
\details{
The data frame \code{dat5} may contain a column with name "1" consisting of 1's
corresponding to the intercept term in the correlation model. If not, this column
will be created. Currently, only the \code{tanh} link is implemented.
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
effect_bin_emp(params = c(0.5, 0.2, 0.5, 0.4), dat5 = dat_zt0)

}
