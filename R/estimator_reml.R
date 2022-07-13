
#' @title Construct a 2 by 2 covariance matrix
#'
#' @details Given variables X and Y, construct the covariance matrix of (X,Y) using
#' their standard deviations \code{sd_x, sd_y} and correlation \code{rho}.
#'
#' @param sd_x the standard deviation of X
#' @param sd_y the standard deviation of X
#' @param rho the correlation between X and Y
#'
#' @return A 2 by 2 numeric matrix.
#'
#' @keywords internal
#'
make_sigma <- function(sd_x, sd_y, rho){
  matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),
         nrow = 2, byrow = TRUE)
}

#' @title Format data into longer format
#'
#' @details Given a data frame \code{dat_xyz}, reformat it so that there is 1 row for each X and 1 row for each Y observation per person.
#' Note: this function is an adapted version of the \code{form} function from the LiquidAssociation package; see
#' \url{https://www.bioconductor.org/packages/release/bioc/html/LiquidAssociation.html}.
#'
#' @param dat_xyz a data frame with columns containing observations of X, Y, ... (in that order)
#'
#' @keywords internal
#'
format_predictors <- function(dat_xyz) {
  # Set number of samples = number of rows
  pts <- nrow(dat_xyz)

  # Use groupid to make 2 entries per sample (1 row for X, 1 row for Y)
  groupid <- rep(1:pts, each = 2)

  # Unravel the X,Y data to match groupid (X1,Y1,X2,Y2,...)
  xy <- matrix(t(dat_xyz[, 1:2]), ncol = 1, nrow = 2 * pts)
  xy <- as.data.frame(xy)
  colnames(xy) <- "xy"

  # Indicator variables for x,y
  i1 <- rep(c(1, 0), pts) # X indicator
  i2 <- rep(c(0, 1), pts) # Y indicator

  # Conditioning variables--if only one, then it's T
  if (ncol(dat_xyz) == 3){
    cnames_z = colnames(dat_xyz)[-(1:2)]
    dat_z <- data.frame(T = rep(dat_xyz[,3], each = 2))

  } else {
    # Otherwise, use the
    cnames_z <- colnames(dat_xyz[,-(1:2)])
    dat_z <- dat_xyz[,-(1:2)]
    dat_z <- lapply(dat_z, rep, each = 2)
    dat_z <- Reduce(cbind, dat_z)
  }

  dat_z1 <- sweep(dat_z, MARGIN=1, i1, `*`)
  colnames(dat_z1) <- paste0(cnames_z, "1")

  dat_z2 <- sweep(dat_z, MARGIN=1, i2, `*`)
  colnames(dat_z2) <- paste0(cnames_z, "2")

  visit <- rep(c(1, 2), pts)

  dat <- as.data.frame(cbind(groupid, xy, i1, i2, dat_z1, dat_z2,  visit))
  return(dat)
}


#' @title Restricted multivariate normal likelihood of CoCoA model parameters
#'
#' @description Restricted multivariate normal likelihood of CoCoA model parameters
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
#'
#' @return The restricted log-likelihood value multiplied by -1.
#'
#' @export
#'
loglik_reml <- function(params, dat_xyz){

  # REML likelihood
  n = nrow(dat_xyz)

  # Input: initial params and covariates for the mean model
  dat_mean = matrix(c(rep(c(1, 0), n), # X indicator
                      rep(c(0, 1), n)), ncol = 2)  # Y indicator

  # input: initial params for the scale/correlation
  sd_x = params[1]
  sd_y = params[2]

  # Split data into outcome/covariates
  # Correlation model predictors
  dat_zs = as.matrix(dat_xyz[,-(1:2)])
  dat_zs = cbind(1, dat_zs)

  # Outcomes stacked
  dat_xyz_long = format_predictors(dat_xyz)
  matrix_xy = as.matrix(dat_xyz_long$xy)

  # Conditional correlations
  params_c = as.matrix(params[-(1:2)])
  rhos = linkinv_tanh(dat_zs %*% params_c)

  # Step 1: calculate beta-hat using the GLS estimator
  list_sigmas = base::lapply(X = rhos, FUN = make_sigma, sd_x = sd_x, sd_y = sd_y) # Subject-specific covariance matrices

  # Reciprocal Condition numbers
  list_sigmas_rcn = base::lapply(X = list_sigmas, FUN = rcond, norm = "1")
  list_sigmas_rcn = unlist(list_sigmas_rcn)
  if (any(list_sigmas_rcn <.Machine$double.eps)){
    message("Subject-specific covariate matrices invertible.")
    # Return lowest possible likelihood value.
    return(-1/.Machine$double.eps)
  }
  # Inverse of block-diagonal matrix can be done block by block
  list_sigmas_inv = lapply(X = list_sigmas, FUN = Matrix::solve) # Subject-specific inverse covariance matrices

  # Stack covariance matrices into big matrix
  matrix_sigma = Matrix::bdiag(list_sigmas)
  matrix_sigma_inv = templateICAr:::bdiag_m(list_sigmas_inv)

  beta_hat = solve(t(dat_mean) %*% matrix_sigma_inv %*% dat_mean) %*% (t(dat_mean)) %*% matrix_sigma_inv %*% matrix_xy
  beta_hat = as.matrix(beta_hat)

  # Step 2: calculate reml likelihood

  # 2.1 calculate log det (sigma)
  log_det_sigma = sum(log(unlist(map(list_sigmas, det))))

  # 2.2 calculate log det (X' sigma_inv X)
  x_sigma_inv_x = t(dat_mean) %*% matrix_sigma_inv %*% dat_mean
  log_det_x_sigmainv_x = log(det(as.matrix(x_sigma_inv_x)))

  # 2.3 calculate exponential term (y-y_hat)' sigma_inv (y - y_hat)
  exp_part = t(matrix_xy - dat_mean %*% beta_hat) %*% matrix_sigma_inv %*% (matrix_xy - dat_mean %*% beta_hat)

  ll = (-1/2)*log_det_sigma + (-1/2)*log_det_x_sigmainv_x - (1/2)*exp_part
  ll = as.numeric(ll)

  return(ll)

}


#' @title Estimate mean parameters from the REML likelihood
#'
#' @details Given an estimate of the variance and correlation parameters, estimate the mean parameters using
#' restricted maximum likelihood.
#'
#' @param pars a numeric vector containing the standard deviation of X and Y, followed by the conditional correlation parameters
#' @param dat_xyz a data frame containing observations of X, Y, Z, ... in that order
#'
#' @return A numeric vector of the REML estimates of the mean parameters.
#'
#' @keywords internal
#'
reml_get_beta <- function(pars, dat_xyz){
  # REML likelihood
  n = nrow(dat_xyz)

  # input: initial params and covariates for the mean model
  dat_mean = matrix(c(rep(c(1, 0), n), # X indicator
                      rep(c(0, 1), n)), ncol = 2)  # Y indiator

  # input: initial params for the scale/correlation
  sd_x = pars[1]
  sd_y = pars[2]

  # Split data into outcome/covariates
  dat_xyz_long = format_predictors(dat_xyz)

  # Correlation model predictors
  dat_zs = as.matrix(dat_xyz[,-(1:2)])
  dat_zs = cbind(1, dat_zs)

  # Outcomes stacked
  matrix_xy = as.matrix(dat_xyz_long$xy)

  # Conditional correlations
  params_c = as.matrix(pars[-(1:2)])
  rhos = linkinv_tanh(dat_zs %*% params_c)

  # Step 1: calculate beta-hat using the GLS estimator
  list_sigmas = purrr::map(rhos, make_sigma, sd_x = sd_x, sd_y = sd_y) # Subject-specific covariance matrices
  list_sigmas_inv = purrr::map(list_sigmas, solve) # Subject-specific inverse covariance matrices

  # Stack covariance matrices into big matrix
  matrix_sigma = templateICAr:::bdiag_m(list_sigmas)
  # Inverse of block-diagonal matrix can be done block by block
  matrix_sigma_inv = templateICAr:::bdiag_m(list_sigmas_inv)

  beta_hat = solve(t(dat_mean) %*% matrix_sigma_inv %*% dat_mean) %*% (t(dat_mean)) %*% matrix_sigma_inv %*% matrix_xy

  return(as.numeric(beta_hat))
}

#' @title Restricted maximum likelihood (REML) estimation of CoCoA model parameters
#'
#' @description Estimate CoCoA model parameters using the restricted maximum likelihood (REML) estimator and the multivariate
#' normal likelihood.
#'
#' @details Given \eqn{n} observations of Z and T, X and Y are bivariate normal with
#' means \code{mu_x, mu_y}, standard deviations \code{sd_x, sd_y}, and correlation
#' \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
#' allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
#' Currently, only the \code{tanh} link is implemented.
#'
#' @param dat_xyz a data frame containing the columns X, Y, Z, T, ... in that order
#' @param init_mle a logical value indicating if the parameter estimates should be initialized at the MLE
#' @param rescale_zt a logical value indicating if the conditional correlation model covariates Z, T, ... should be rescaled (z-scored)
#' @param optim_type a character value. If \code{optim_type = "optim"}, the usual \code{optim} optimization routine is used.
#' Otherwise the "Rvmmin" method from the optimr package is used.
#'
#' @return A list containing \code{params}, a numeric matrix of the estimated parameters, \code{out2}, a table of the estimated parameters,
#' \code{runtime}, the runtime of the optimization routine, and \code{loglik}, the value of the final restricted log-likelihood.
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
#' get_params_reml(dat_xyz)
#'
#' @export
#'
get_params_reml <- function(dat_xyz, init_mle = TRUE, rescale_zt = FALSE, optim_type = "optim"){
  n = nrow(dat_xyz)
  p = ncol(dat_xyz) - 1 # Number of params to estimate

  if (init_mle){
    # Get estimates from MLE
    params_mle = get_params_mle(dat_xyz)$params
    start_vals = params_mle[-(1:2)] # sd_x, sd_y, cc params
  } else {
    # Otherwise, initialize at empirical SD
    start_vals = c(
      sd(dat_xyz$X),
      sd(dat_xyz$Y),
      # Initial estimates of the CC params
      rep(0.1, p))
  }

  dat_xyz2 <- dat_xyz

  if (rescale_zt){
    xy_ind = which(colnames(dat_xyz2) %in% c("X", "Y"))
    means_zt = apply(dat_xyz2[,-xy_ind], 2, mean)
    sds_zt = apply(dat_xyz2[,-xy_ind], 2, sd)

    dat_xyz2[,-xy_ind] <- apply(dat_xyz2[,-xy_ind], 2, scale2)
  }

  # Maximize the log-likelihood
  t0 = Sys.time()
  if (optim_type == "optim"){
    opt = optim(par=start_vals, loglik_reml, dat_xyz = dat_xyz2,
                # Letting bounds be 0 gives convergence errors
                method="L-BFGS-B", lower=c(1e-3, 1e-3, rep(-5, p)), upper=c(Inf, Inf, rep(5, p)), hessian = TRUE,
                control = list(fnscale = -1))

    # Log-likelihood
    ll = 1*(opt$value)
  } else {
    # Use optimr, which is faster but may be less good at finding the global maximum.
    opt = optimr::optimr(par=start_vals, loglik_reml, dat_xyz = dat_xyz2,
                 # Letting bounds be 0 gives convergence errors
                 method="Rvmmin", lower=c(1e-3, 1e-3, rep(-5, p)), upper=c(Inf, Inf, rep(5, p)), hessian = TRUE,
                 control = list(maximize = TRUE))

    # Log-likelihood
    ll = -1*(opt$value)

  }
  t1 = Sys.time() - t0

  # MLEs of alpha, beta
  pars = opt$par

  # Calculate beta from the V estimates
  final_beta_hat = reml_get_beta(pars, dat_xyz)

  # Output
  param_names <- grk_letters[1:(length(pars)-2)]
  var_names <- c("1", colnames(dat_xyz)[-(1:2)])

  if (rescale_zt){
    alpha0 = pars[3]
    beta_gamma0 = pars[-(1:3)]

    alpha2 = alpha0 - sum((means_zt*beta_gamma0)/sds_zt)
    beta_gamma2 = beta_gamma0/sds_zt

    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", "sigma_y", param_names),
                      Variable = c("X", "Y", "X", "Y", var_names),
                      estimate = c(final_beta_hat, pars[1:2], alpha2, beta_gamma2))
    out2 = dplyr::mutate(out2,
                         san.se = as.numeric(NA),
                         wald = as.numeric(NA),
                         p = as.numeric(NA),
                         Model =c(rep("mean", 2), rep("scale", 2), rep("correlation", length(param_names))))

  } else {
    out2 = data.frame(Coefficient = c("mu_x", "mu_y", "sigma_x", "sigma_y", param_names),
                      Variable = c("X", "Y", "X", "Y", var_names),
                      estimate = c(final_beta_hat, pars))
    out2 = dplyr::mutate(out2,
                         san.se = as.numeric(NA),
                         wald = as.numeric(NA),
                         p = as.numeric(NA),
                         Model =c(rep("mean", 2), rep("scale", 2), rep("correlation", length(param_names))))
  }


  return(list(params = as.matrix(pars),
              out2 = out2,
              runtime = t1,
              loglik = ll))
}
