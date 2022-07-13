
#' @title Bound values in [-1,1] away from -1, 0, and 1
#'
#' @param x a numeric vector with elements in [-1,1].
#'
#' @return The bounded value.
#'
#' @examples
#' linkfun_tanh(1) #Inf
#' linkfun_tanh(bound_away(1))
#'
#' @keywords internal
#'
bound_away <- function(x){
  ind_pos = which(x >= 0)
  ind_neg = which(x < 0)
  out = x

  # x in [0,1]
  out[ind_pos] = pmax(x[ind_pos], .Machine$double.eps)
  out[ind_pos]= pmin(out[ind_pos], 1 - .Machine$double.eps)
  # x in [-1, 0)
  out[ind_neg] = pmin(x[ind_neg], -1*.Machine$double.eps)
  out[ind_neg]= pmax(out[ind_neg], -1 + .Machine$double.eps)

  return(out)
}

#' @title Beta log-likelihood
#'
#' @param params a numeric vector of length 2 for the Beta shape parameters
#' @param dat_x a numeric vector of observations
#'
#' @return The Beta log-likelihood given the data \code{dat_x} and parameters \code{params}.
#'
#' @keywords internal
#'
loglik_beta <- function(params = c(1,1), dat_x){
  sum(log(dbeta(dat_x, shape1 = params[1], shape2 = params[2])))
}

#' @title Beta maximum likelihood estimator
#'
#' @param dat_x a numeric vector of observations
#'
#' @return The maximum likelihood estimates of the Beta parameters given the data \code{dat_x}.
#'
#' @keywords internal
#'
get_beta_params <- function(dat_x){
  opt = optim(par = c(1,1),
              fn = loglik_beta,
              dat_x = dat_x,
              method = "L-BFGS-B",
              lower = c(.Machine$double.eps, .Machine$double.eps),
              upper = c(10,10),
              control = list(fnscale = -1))

  mle = c(shape1 = opt$par[1], shape2 = opt$par[2])
  return(mle)
}

#' @title Gamma log-likelihood
#'
#' @param params a numeric vector of length 2 for the Gamma shape parameters
#' @param dat_x a numeric vector of observations
#'
#' @return The Gamma log-likelihood given the data \code{dat_x} and parameters \code{params}.
#'
#' @keywords internal
#'
loglik_gamma <- function(params = c(1,1), dat_x){
  sum(log(dgamma(dat_x, shape = params[1], rate = params[2])))
}

#' @title Gamma maximum likelihood estimator
#'
#' @param dat_x a numeric vector of observations
#'
#' @return The maximum likelihood estimates of the Gamma parameters given the data \code{dat_x}.
#'
#' @keywords internal
#'
get_gamma_params <- function(dat_x){

  opt = optim(par = c(1,1),
              fn = loglik_gamma,
              dat_x = dat_x,
              method = "L-BFGS-B",
              lower = c(.Machine$double.eps, .Machine$double.eps),
              upper = c(10,10),
              control = list(fnscale = -1))

  mle = c(shape = opt$par[1], rate = opt$par[2])
  return(mle)
}

#' @title Transform (X,Y) to a Gaussian copula with N(0,1) marginals
#'
#' @description Transform (X,Y) to a Gaussian copula with N(0,1) marginals using a quantile transform.
#'
#' @details Given \eqn{n} observations of X, Y, and their distribution functions \code{F_x, F_y}, use a quantile
#' transform so that the transformed (X,Y) follows a Gaussian copula with marginal N(0,1) distributions. If \code{F_x, F_y}
#' are not supplied, then the empirical distribution function is used.
#'
#' @param dat_xyz_u a data frame where the first column is X (untransformed) and the second column is Y (untransformed)
#' @param F_x the distribution function of X
#' @param F_y the distribution function of Y
#'
#' @return The original data frame with the columns of X and Y replaced by their transformed values.
#'
#' @examples
#' T = as.numeric(rbernoulli(n = n, p = 0.5))
#' Z = rep(NA, times = n)
#' for (i in 1:n){
#' Z[i] <- if (T[i] == 1){
#' # Z|T=1 ~ Bernoulli(0.3)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
#' } else {
#' # Z|T=0 ~ Bernoulli(0.9)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
#' }
#' }
#' # Simulate X,Y given Z,T, and Z*T using a normal copula with gamma and beta marginals.
#' dat_zt0 <- data.frame(Z = Z, T = T, ZT = Z*T)
#' XY = sim_xy2_mvn(theta0 = c(0, 0, 0, 0), dat_zt0 = dat_zt0,
#' mu_x = 0, mu_y = 0, sd_x = 1, sd_y = 1)
#' X_t = qgamma(bound_away(pnorm(XY$x, mean = 0, sd = 1)), shape = 2,  rate = 1) # Gamma(2,1)
#' Y_t = qbeta(bound_away(pnorm(XY$y, mean = 0, sd = 1)), shape1 = 2, shape2 = 2) # Beta(2,2)
#' dat_xyz = data.frame(X = X_t, Y = Y_t)
#' dat_xyz = cbind(dat_xyz, dat_zt0)
#'
#' # If the marginal distributions are knonw, but arams unknown, estimate with MLE
#' mle_x = get_gamma_params(dat_x = X_t)
#' F_x = function(p) pgamma(p, shape=mle_x[1],  rate=mle_x[2])
#' mle_y = get_beta_params(dat_x = Y_t)
#' F_y = function(p) pbeta(p,shape1 = mle_y[1], shape2 = mle_y[1])
#' dat_xyz = gcop_transform2(dat_xyz, F_x = F_x, F_y = F_y)
#'
#' @export
#'
gcop_transform2 <- function(dat_xyz_u, F_x = NULL, F_y = NULL){
  # Transform X and Y to uniform random variables
  x <- dat_xyz_u[,1]
  y <- dat_xyz_u[,2]

  # If F_x is not supplied, use empirical CDF
  if (is.null(F_x)){
    F_x <- ecdf(x)
  }
  u_x <- bound_away(F_x(x))
  # If F_y is not supplied, use empirical CDF
  if (is.null(F_y)){
    F_y <- ecdf(y)
  }
  u_y <- bound_away(F_y(y))

  # Transform u_x and u_y to x_t and y_t using N(0,1) quantile function
  x_t <- qnorm(p = u_x, mean = 0, sd = 1)
  y_t <- qnorm(p = u_y, mean = 0, sd = 1)

  # Output the data with the transformed columns
  dat_xyz_t <- dat_xyz_u
  dat_xyz_t[,1] <- x_t
  dat_xyz_t[,2] <- y_t

  return(dat_xyz_t)
}
