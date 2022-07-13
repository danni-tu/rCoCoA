
#' @title Multivariate normal likelihood of CoCoA model parameters
#'
#' @description Multivariate normal likelihood of CoCoA model parameters
#'
#' @details Given \eqn{n} observations of Z and T, X and Y are bivariate normal with
#' means \code{mu_x, mu_y}, standard deviations \code{sd_x, sd_y}, and correlation
#' \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
#' allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
#' Currently, only the \code{tanh} link is implemented. Apart from the X and Y columns, \code{dat_xyz} may include
#' multiple covariates Z, T, ... for the conditional correlation model. The number of these covariates should be
#' equal to the length of \code{params} minus 4.
#'
#' @param params a numeric vector of the form \code{(mu_x, mu_y, sd_x, sd_y, alpha, beta, ...)} in that order
#' @param dat_xyz a data frame containing the columns X, Y, Z, T, ... in that order
#'
#' @return The log-likelihood value multiplied by -1.
#'
#'
loglik <- function(params = c(0, 0, 1, 1, 0.1, 0.1, 0.1), dat_xyz){

  n = nrow(dat_xyz)

  # Convert X,Y, and conditional correlation covariates to matrices
  x = as.matrix(dat_xyz[,1])
  y = as.matrix(dat_xyz[,2])
  # Matrix of z's, t's, etc. as a design matrix
  dat_zs = as.matrix(dat_xyz[,-(1:2)])
  dat_zs <- cbind(1, dat_zs)

  # Parameters provided
  n_params = ncol(dat_zs) + 4
  if (length(params) != n_params){
    warning("Warning: Not enough parameters provided, defaulting to 0.1")
    params = rep(0.1, n_params)
    params = as.matrix(params, ncol = 1)
  }

  # Split parameters:
  params_m = as.matrix(params[1:2]) # Means
  params_s = as.matrix(params[3:4]) # Variances
  params_c = as.matrix(params[-(1:4)]) # Conditional correlations

  # Conditional correlation vector
  rho = linkinv_tanh(dat_zs %*% params_c)

  # Log-likelihood
  # ll = 0
  # for (i in 1:n){
  #   # Conditional correlation matrix
  #   Sigma = c(params_s[1]^2,
  #             params_s[1]*params_s[2]*rho[i],
  #             params_s[1]*params_s[2]*rho[i],
  #             params_s[2]^2) %>%
  #     matrix(ncol = 2)
  #
  #   # Add to log-likelihood
  #   ll = ll + mvtnorm::dmvnorm(c(x[i],y[i]),  mean = params[1:2],
  #                              sigma = Sigma, log = TRUE, checkSymmetry = TRUE)
  # }
  # logl = -1*ll

  # Calculate log-likelihood
  ll1 = n*log(2*pi) + (1/2)*sum(log(1-rho^2)) + n*log(params_s[1]) + n*log(params_s[2])
  ll2_inner = ((x - params_m[1])/params_s[1])^2 +
    ((y - params_m[2])/params_s[2])^2 -
    2*rho*((x - params_m[1])/params_s[1])*((y - params_m[2])/params_s[2])
  ll2 = (1/2)*sum((1/(1-rho^2))*ll2_inner)
  logl = ll1 + ll2

  return(logl)
}


#' @title Maximum likelihood estimation of CoCoA model parameters
#'
#' @description Estimate CoCoA model parameters using the maximum likelihood estimator and the multivariate
#' normal likelihood.
#'
#' @details Given \eqn{n} observations of Z and T, X and Y are bivariate normal with
#' means \code{mu_x, mu_y}, standard deviations \code{sd_x, sd_y}, and correlation
#' \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
#' allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
#' Currently, only the \code{tanh} link is implemented.
#'
#' @param dat_xyz a data frame containing the columns X, Y, Z, T, ... in that order
#' @param rescale_zt a logical value indicating if the conditional correlation model covariates Z, T, ... should be rescaled (z-scored)
#'
#' @importFrom dplyr mutate
#'
#' @return A list containing \code{params}, a numeric matrix of the estimated parameters, \code{out2}, a table of the estimated parameters,
#' and \code{runtime}, the runtime of the optimization routine.
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
#' dat_xyz = data.frame(X = XY$x, Y = XY$y)
#' dat_xyz = cbind(dat_xyz, dat_zt0)
#' get_params_mle(dat_xyz)
#'
#' @export
#'
get_params_mle <- function(dat_xyz, rescale_zt = FALSE){

  n = nrow(dat_xyz) # Number of observations
  p = ncol(dat_xyz) - 1 # Number of params to estimate

  dat_xyz2 = as.data.frame(dat_xyz)

  if (rescale_zt){
    # Z-score the covariates, but save the means and SDs
    xy_ind = which(colnames(dat_xyz2) %in% c("X", "Y"))
    means_zt = apply(dat_xyz2[,-xy_ind], 2, mean)
    sds_zt = apply(dat_xyz2[,-xy_ind], 2, sd)

    dat_xyz2[,-xy_ind] <- apply(dat_xyz2[,-xy_ind], 2, scale2)
  }

  # Optimization without constraints
  # Initialize means/scales at the empirical means/sds
  start_vals = c(mean(dat_xyz2$X),
                 mean(dat_xyz2$Y),
                 sd(dat_xyz2$X),
                 sd(dat_xyz2$Y),
                 # Initialize conditional correlation params at 0.001
                 rep(0.001, p))

  t0 = Sys.time()
  opt = optim(par = start_vals, f = loglik, gr = NULL,
              control = list(reltol = 1e-15),
              dat_xyz = data.matrix(dat_xyz2))
  t1 = Sys.time() - t0

  # MLEs of means, variance components, and conditional correlation parameters
  pars = opt$par
  # Log-likelihood
  ll = -1*(opt$value)

  # Output
  param_names <- rCoCoA:::grk_letters[1:(length(pars)-4)]
  var_names <- c("1", colnames(dat_xyz)[-(1:2)])

  # Predicted correlations
  # rho_hats = as.matrix(cbind(1, dat_xyz2[,-(1:2)]))  %*%  as.matrix(pars[-(1:4)], ncol = 1)

  # Back-transform if rescaling Z,T
  if (rescale_zt){
    alpha0 = pars[5]
    beta_gamma0 = pars[-(1:5)]

    alpha2 = alpha0 - sum((means_zt*beta_gamma0)/sds_zt)
    beta_gamma2 = beta_gamma0/sds_zt

    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", "sigma_y", param_names),
                      Variable = c("X", "Y", "X", "Y", var_names),
                      estimate = c(pars[1:4], alpha2, beta_gamma2))
    out2 = dplyr::mutate(out2,
                         san.se = as.numeric(NA),
             wald = as.numeric(NA),
             p = as.numeric(NA),
             Model =c(rep("mean", 2), rep("scale", 2), rep("correlation", length(param_names))))

  } else {
    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", "sigma_y", param_names),
                      Variable = c("X", "Y", "X", "Y", var_names),
                      estimate = pars)
    out2 = dplyr::mutate(out2,
                  san.se = as.numeric(NA),
                  wald = as.numeric(NA),
                  p = as.numeric(NA),
                  Model =c(rep("mean", 2), rep("scale", 2), rep("correlation", length(param_names))))
  }

  return(list(params = as.matrix(pars), # CAUTION: params are Untransformed parameters, if rescale_zt = TRUE.
              out2 = out2,
              runtime = t1))
}
