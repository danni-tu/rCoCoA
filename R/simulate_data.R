#' @title Simulate multivariate normal data in CoCoA model
#'
#' @description Simulate multivariate normal (X,Y) given observations of (Z,T) and the multivariate Gaussian CoCoA model.
#'
#' @details Given \eqn{n} observations of Z and T, generate X and Y according to the model \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...)}
#' where \eqn{(X,Y)} follow a multivariate normal distribution with means \code{mu_x, mu_y} and standard deviations \code{sd_x, sd_y}. The inverse link function \eqn{g^{-1}}
#' allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
#'
#' @param theta0 numeric vector of true parameters
#' @param dat_zt0 data.frame containing observations of Z and T
#' @param mu_x mean of X
#' @param mu_y means of Y
#' @param sd_x standard deviation of X
#' @param sd_y standard deviation of Y
#' @param cor_link character vector corresponding to the correlation link function (default = "tanh")
#'
#' @return Generated X and Y values based on the multivariate Gaussian model.
#'
#' @examples
#' T = as.numeric(purrr::rbernoulli(n = 100, p = 0.5))
#' Z = rep(NA, times = 100)
#' for (i in 1:100){
#' Z[i] <- if (T[i] == 1){
#' # Z|T=1 ~ Bernoulli(0.3)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
#' } else {
#' # Z|T=0 ~ Bernoulli(0.9)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
#' }
#' }
#' dat_zt0 <- data.frame(Z = Z, T = T, ZT = Z*T)
#' XY = sim_xy2_mvn(theta0 = c(0, 0, 0, 0), dat_zt0 = dat_zt0,
#' mu_x = 1, mu_y = 1, sd_x = 1, sd_y = 2, cor_link = "tanh")
#'
#'
#' @export
#'
sim_xy2_mvn <- function(theta0, dat_zt0, mu_x = 0, mu_y = 0, sd_x = 1, sd_y = 1, cor_link = "tanh"){

  n = nrow(dat_zt0)

  # Ensure number of parameters = number of variables - 1
  stopifnot(length(theta0) == (1 + ncol(dat_zt0)))

  # Append column of 1's for the intercept
  dat_zt <- as.matrix(cbind(1,dat_zt0))

  theta = matrix(theta0, ncol = 1)

  # Vector of conditional correlations
  cor_link_inv_fun = switch(cor_link,
                            "tanh" = linkinv_tanh,
                            "identity" = linkinv_tid,
                            "logit" = linkinv_logit,
                            "probit" = linkinv_probit,
                            "cloglog" = linkinv_cloglog)
  rhos = cor_link_inv_fun(dat_zt %*% theta)

  # Simulate X,Y based on MVN data
  x = rep(NA, times = n)
  y = rep(NA, times = n)

  for (i in 1:n){
    Sigma = matrix(c(sd_x^2, sd_x*sd_y*rhos[i],
                     sd_x*sd_y*rhos[i], sd_y^2),
                   byrow = TRUE, nrow = 2)

    samp = MASS::mvrnorm(n = 1, mu = c(mu_x,mu_y), Sigma = Sigma)
    x[i] = samp[1]
    y[i] = samp[2]
  }

  return(list(x = x, y = y))
}


#' @title Simulate copula-valued data in CoCoA model
#'
#' @description Simulate (X,Y) from a Gaussian copula with Beta and Gamma margins given observations of (Z,T) and the CoCoA model.
#'
#' @details Given \eqn{n} observations of Z and T, generate X and Y according to the model \deqn{Corr_s(X,Y) = g^{-1}(\alpha_s + \beta_s Z + \gamma_s T + ...)}
#' where \eqn{Corr_s(X,Y)} is the Spearman correlation, X follows a Gamma(2,1) marginal, and Y follows a Beta(2,2) marginal.
#'
#' @param theta0 numeric vector of true parameters
#' @param dat_zt0 data.frame containing observations of Z and T
#' @param cor_link character vector corresponding to the correlation link function (default = "tanh")
#'
#' @return Generated X and Y values based on the copula model.
#'
#' @examples
#' T = as.numeric(purrr::rbernoulli(n = 100, p = 0.5))
#' Z = rep(NA, times = 100)
#' for (i in 1:100){
#' Z[i] <- if (T[i] == 1){
#' # Z|T=1 ~ Bernoulli(0.3)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
#' } else {
#' # Z|T=0 ~ Bernoulli(0.9)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
#' }
#' }
#' dat_zt0 <- data.frame(Z = Z, T = T, ZT = Z*T)
#' XY = sim_xy2_cop(theta0 = c(0, 0, 0, 0), dat_zt0 = dat_zt0, cor_link = "tanh")
#'
#'
#' @export
#'
sim_xy2_cop <- function(theta0, dat_zt0, cor_link = "tanh"){

  n = nrow(dat_zt0)

  # Ensure parameters = number of variables - 1
  stopifnot(length(theta0) == (1 + ncol(dat_zt0)))

  # Append column of 1's for the intercept
  dat_zt <- as.matrix(cbind(1, dat_zt0))

  theta = matrix(theta0, ncol = 1)

  # Vector of conditional correlations--here, spearmans
  # rhos = tanh(dat_zt %*% theta)
  cor_link_inv_fun = switch(cor_link,
                            "tanh" = linkinv_tanh,
                            "identity" = linkinv_tid,
                            "logit" = linkinv_logit,
                            "probit" = linkinv_probit,
                            "cloglog" = linkinv_cloglog)
  rhos = cor_link_inv_fun(dat_zt %*% theta)

  # Simulate X,Y based on MVN data
  x = rep(NA, times = n)
  y = rep(NA, times = n)

  for (i in 1:n){

    # Elliptical copula
    myCop <- copula::normalCopula(param = rhos[i], dim = 2, dispstr = "un")

    # Create multivariate distribution using myCop
    myMvd <- copula::mvdc(copula = myCop, margins = c("gamma", "beta"),
                  paramMargins=list(list(shape=2,  scale=1),
                                    list(shape1 = 2, shape2 = 2)))
    samp <- copula::rMvdc(n = 1, myMvd)

    x[i] = samp[1]
    y[i] = samp[2]
  }

  return(list(x = x, y = y))
}


#' @title Obtain Pearson correlation parameters from copula-based CoCoA model (internal function)
#'
#' @description Obtain Pearson correlation parameters from copula-based CoCoA model (internal function)
#'
#' @details Given \eqn{n} observations of binary Z and T, and X and Y generated according to the model \deqn{Corr_s(X,Y) = g^{-1}(\alpha_s + \beta_s Z + \gamma_s T + \delta_s Z*T)}
#' where \eqn{Corr_s(X,Y)} is the Spearman correlation, X follows a Gamma(2,1) marginal, and Y follows a Beta(2,2) marginal, recover the parameters
#' \eqn{(\alpha, \beta, \gamma, \delta)} corresponding to the *Pearson* correlation model \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + \delta Z*T).}
#'
#' @param alpha numeric parameter in copula model (see details)
#' @param beta numeric parameter in copula model (see details)
#' @param gamma numeric parameter in copula model (see details)
#' @param delta numeric parameter in copula model (see details)
#' @param cor_link character vector corresponding to the correlation link function (default = "tanh")
#'
#' @return A numeric matrix containing the model parameters corresponding to the Pearson correlation model.
#'
#' @examples
#' theta_true = true_params_cop(0, 0, 0, 0, cor_link = "tanh")
#'
true_params_cop <- function(alpha = 0.5, beta = 0.8, gamma = 0.1, delta = 0.5, cor_link = "tanh"){
  n = 100000

  # Simulate Z,T
  T = as.numeric(purrr::rbernoulli(n = n, p = 0.5))
  Z = rep(NA, times = n)

  # Z,T
  for (i in 1:n){
    Z[i] <- if (T[i] == 1){
      # Z|T=1 ~ Bernoulli(0.3)
      as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
    } else {
      # Z|T=0 ~ Bernoulli(0.9)
      as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
    }
  }

  # Simulate X,Y
  dat_zt0 <- data.frame(Z = Z, T = T, ZT = Z*T)
  XY = sim_xy2_cop(theta0 = c(alpha, beta, gamma, delta), dat_zt0 = dat_zt0, cor_link = cor_link)
  dat_xy = data.frame(X = XY$x, Y = XY$y)

  # Estimate pearson correlation
  r1 = cor(dat_xy[which(Z == 1 & T == 1),])[1,2] # alpha + beta + gamma + delta
  r2 = cor(dat_xy[which(Z == 1 & T == 0),])[1,2] # alpha + beta
  r3 = cor(dat_xy[which(Z == 0 & T == 1),])[1,2] # alpha + gamma
  r4 = cor(dat_xy[which(Z == 0 & T == 0),])[1,2] # alpha

  cor_link_fun = switch(cor_link,
                        "tanh" = linkfun_tanh,
                        "identity" = linkfun_tid,
                        "logit" = linkfun_logit,
                        "probit" = linkfun_probit,
                        "cloglog" = linkfun_cloglog)
  b = matrix(cor_link_fun(c(r1, r2, r3, r4)), ncol = 1)

  # Find params
  a = matrix(c(1, 1, 1, 1,
               1, 1, 0, 0,
               1, 0, 1, 0,
               1, 0, 0, 0), nrow = 4)
  theta = solve(a = a, b = b)
  return(theta)

}
