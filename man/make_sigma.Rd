% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimator_reml.R
\name{make_sigma}
\alias{make_sigma}
\title{Construct a 2 by 2 covariance matrix}
\usage{
make_sigma(sd_x, sd_y, rho)
}
\arguments{
\item{sd_x}{the standard deviation of X}

\item{sd_y}{the standard deviation of X}

\item{rho}{the correlation between X and Y}
}
\value{
A 2 by 2 numeric matrix.
}
\description{
Construct a 2 by 2 covariance matrix
}
\details{
Given variables X and Y, construct the covariance matrix of (X,Y) using
their standard deviations \code{sd_x, sd_y} and correlation \code{rho}.
}
\keyword{internal}
