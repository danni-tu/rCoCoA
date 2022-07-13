
#' @title Second order generalized estimating equations (GEE2) estimation of CoCoA model parameters
#'
#' @description Estimate CoCoA model parameters using second order generalized estimating equations (GEE2).
#'
#' @details Given \eqn{n} observations of Z and T, the conditional mean, scale, and correlation of X and Y
#' is described using second order generalized estimating equations (GEE2), introduced
#' in Yan and Fine (2004). In particular, the GEE2 model includes the correlation model
#' \deqn{Corr(X,Y) = g^{-1}(\alpha + \beta Z + \gamma T + ...).} The inverse link function \eqn{g^{-1}}
#' allows for unconstrained parameter estimation while ensuring that the predicted correlation is bounded.
#' This function is a wrapper around the \code{geese()} function from geepack, and is an
#' expanded version of \code{CNM.full()} from the LiquidAssociation package.
#'
#' If \code{no_meanmod = TRUE}, then only an intercept term will be estimated in the mean model (i.e.,
#' all participants will have the same predicted mean of X and the same predicted mean of Y).
#' If \code{no_meanscalemod = TRUE}, then only an intercept term will be estimated in both
#' the mean and scale models.
#'
#' @param dat_xyz a data frame containing the columns X, Y, Z, T, ... in that order
#' @param jack_se a logical value indicating if the jackknife standard errors should be computed
#' @param init_mle a logical value indicating if the parameter estimates should be initialized at the MLE
#' @param no_meanmod a logical value indicating if the mean model should be estimated (see Details)
#' @param no_meanscalemod a logical value indicating if the mean and scale models should be estimated (see Details)
#' @param rescale_xy a logical value indicating if X and Y should be rescaled (z-scored)
#' @param rescale_zt a logical value indicating if the conditional correlation model covariates Z, T, ... should be rescaled (z-scored)
#' @param cor_link a character string indicating the correlation model link function; passed on to \code{geese()}
#' @param sca_link a character string indicating the scale model link function; passed on to \code{geese()}
#'
#' @return A list containing \code{params}, a numeric matrix of the estimated parameters, \code{out2}, a table of the estimated parameters,
#' \code{runtime}, the runtime of the optimization routine.
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
#' get_params_gee(dat_xyz)
#'
#' @export
#'
get_params_gee <- function(dat_xyz, jack_se = FALSE, init_mle = TRUE,
                           no_meanmod = FALSE, no_meanscalemod = FALSE,
                           rescale_xy = FALSE, rescale_zt = FALSE,
                           cor_link = "tanh", sca_link = "log"){

  # Tanh link is the fisherz link for geepack
  if (cor_link == "tanh") {cor_link <- "fisherz"}
  stopifnot(cor_link %in% c("fisherz", "identity", "logit", "probit", "cloglog"))
  stopifnot(sca_link %in% c("identity", "log"))

  # Filter to non-missing rows
  dat_xyz2 <- dplyr::filter(dat_xyz, complete.cases(dat_xyz))

  # Colnames of the non-X/Y variables
  cnames_zt = colnames(dat_xyz2)[-(1:2)]

  # Rescale X,Y?
  xy_ind = which(colnames(dat_xyz2) %in% c("X", "Y"))
  if (rescale_xy){
    means_xy = apply(dat_xyz2[,xy_ind], 2, mean)
    sds_xy = apply(dat_xyz2[,xy_ind], 2, sd)

    dat_xyz2[,xy_ind] <- apply(dat_xyz2[,xy_ind], 2, scale2)
  }

  # Rescale Z,T?
  if (rescale_zt){
    means_zt = apply(dat_xyz2[,-xy_ind], 2, mean)
    sds_zt = apply(dat_xyz2[,-xy_ind], 2, sd)

    dat_xyz2[,-xy_ind] <- apply(dat_xyz2[,-xy_ind], 2, scale2)
  }

  # Turn the data into long format
  dat_mod <- format_predictors(dat_xyz2)

  # Columns corresponding to the conditioning vars
  z_cnames1 <- paste0(cnames_zt, 1)
  z_cnames2 <- paste0(cnames_zt, 2)
  z_cnames <- c(z_cnames1, z_cnames2)

  # Linear predictor for the correlation
  zcor <- cbind(1, dat_xyz2[,-xy_ind])
  colnames(zcor) <- c("1", cnames_zt)

  # GEE model
  # Create formulas for mean and scale link functions
  form_mean = paste0("xy ~ ", paste0(z_cnames, collapse = " + "), " - 1")
  form_scale = paste0(" ~ i1 + i2 + ", paste0(z_cnames, collapse = " + "), " - 1")

  # if (init_mle){
  #   # Get estimates from MLE
  #   params_mle = get_params_mle(dat_xyz2)$params
  #   b_init = params_mle[1:2] # Mean
  #   # gm_init = params_mle[3:4] # Scale
  #   alpha_init = params_mle[-(1:4)] # Correlation
  # } else {
  #   alpha_init = rep(0, ncol(zcor))
  # }

  t0 = Sys.time()
  if (no_meanmod){
    # Only include Z,T in estimation of the scale and correlation parameters
    # Fit just the intercept for X/Y mean
    fit <- geepack::geese(formula = xy ~ i1 + i2 - 1, id = groupid, data = dat_mod,
                 sformula = as.formula(form_scale), waves = visit,
                 family = gaussian,
                 sca.link = sca_link, cor.link = cor_link, corstr = "userdefined",
                 # Estimate jackkknife SEs if specified
                 control = geepack::geese.control(epsilon = 1e-15, maxit = 35,
                                         trace = FALSE, fij = jack_se,
                                         jack = jack_se, j1s = jack_se),
                 # Initialize values if specified
                 # b = b_init, #gm = gm_init,
                 # alpha = alpha_init,
                 # Linear predictor for the correlation component
                 zcor = zcor)

  } else {
    # Only include Z,T in estimation of the correlation parameters
    # Fit just the intercept for X/Y mean/scale
    if (no_meanscalemod){
      # No scale mod?
      fit <- geepack::geese(formula = xy ~ i1 + i2 - 1, id = groupid, data = dat_mod,
                   sformula = ~ i1 + i2 - 1, waves = visit,
                   sca.link = sca_link, cor.link = cor_link, corstr = "userdefined",
                   # Estimate jackkknife SEs if specified
                   control = geepack::geese.control(epsilon = 1e-15, maxit = 35,
                                           trace = FALSE, fij = jack_se,
                                           jack = jack_se, j1s = jack_se),
                   # Initialize values if specified
                   # b = b_init, #gm = gm_init,
                   # alpha = alpha_init,
                   # Linear predictor for the correlation component
                   zcor = zcor)
    } else {
      # Otherwise, include Z,T in the mean estimating equation
      fit <- geepack::geese(formula = as.formula(form_mean), id = groupid, data = dat_mod,
                   sformula = as.formula(form_scale), waves = visit,
                   sca.link = sca_link, cor.link = cor_link, corstr = "userdefined",
                   # Estimate jackkknife SEs if specified
                   control = geepack::geese.control(epsilon = 1e-15, maxit = 35,
                                           trace = FALSE, fij = jack_se,
                                           jack = jack_se, j1s = jack_se),
                   # Initialize values if specified
                   #b = b_init, #gm = gm_init,
                   # alpha = alpha_init,
                   # Linear predictor for the correlation component
                   zcor = zcor)
    }
  }
  t1 = Sys.time() - t0

  fit_coefs <- summary(fit)

  # Coefficients
  out1 = tibble::rownames_to_column(fit_coefs[["mean"]], "Variable")
  out1$Model <- "mean"

  out2 = tibble::rownames_to_column(fit_coefs[["correlation"]], "Variable")
  out2$Model <- "correlation"

  out3 = tibble::rownames_to_column(fit_coefs[["scale"]], "Variable")
  out3$Model <- "scale"

  out = rbind(out1, out2, out3)
  out = dplyr::mutate_if(out, is.numeric, round, 3) # Round numeric columns

  # Output model name
  Model <- paste("Model: Corr(", colnames(dat_xyz)[1], ",",
                 colnames(dat_xyz)[2], "|",
                 paste0(colnames(dat_xyz)[-(1:2)], collapse = ","), ")",
                 sep = "")

  # If Z,T were re-scaled, convert the coefficients back to original scale
  if (rescale_zt){
    alpha0 = out2$estimate[out2$Variable == "1"]
    beta_gamma0 = out2$estimate[out2$Variable != "1"]

    alpha2 = alpha0 - sum((means_zt*beta_gamma0)/sds_zt)
    beta_gamma2 = beta_gamma0/sds_zt

    # Only correlation model output
    out2 = dplyr::mutate(out2,
        Coefficient = grk_letters[1:nrow(out2)],
        estimate = c(alpha2, beta_gamma2),
        Rescaled = "Fit on scaled Z",
        san.se = as.numeric(NA))
    out2 = dplyr::select(out2, Coefficient, everything())
    out2 = dplyr::mutate_if(out2, is.numeric, round, 3) # Round numeric columns

  } else {
    # Only correlation model output
    out2 = dplyr::mutate(out2,
        Coefficient = grk_letters[1:nrow(out2)],
        Rescaled = "Fit on orig Z")
    out2 = dplyr::select(out2, Coefficient, everything())
    out2 = dplyr::mutate_if(out2, is.numeric, round, 3) # Round numeric columns

  }

  # Output covariance matrix of the estimated correlation parameters (sandwich)
  v_corr_params = fit$valpha

  return(list(Model = Model,
              params = as.matrix(out2$estimate),
              coefs = out,
              out2 = out2,
              v_corr_params = v_corr_params,
              runtime = t1,
              cor_link = cor_link,
              sca_link = sca_link))
}
