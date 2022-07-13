
#' @title Association Size of Binary Variable T
#'
#' @details The data frame \code{dat5} may contain a column with name "1" consisting of 1's
#' corresponding to the intercept term in the correlation model. If not, this column
#' will be created. Currently, only the \code{tanh} link is implemented.
#'
#' @param params a numeric vector containing the estimated correlation parameters in the CoCoA model
#' @param dat5 a data frame containing only the covariates in the correlation model (i.e., excluding X and Y)
#'
#' @return A numeric vector containing the marginal and conditional association sizes.
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
#' effect_bin_emp(params = c(0.5, 0.2, 0.5, 0.4), dat5 = dat_zt0)
#'
#' @export
#'
effect_bin_emp <- function(params, dat5){

  n_p = length(params)
  params2 = matrix(params, ncol = 1)
  stopifnot("T" %in% colnames(dat5))

  # Add intercept if not in the data
  if (!("1" %in% colnames(dat5))){
    dat_zt <- cbind(1, dat5)
    colnames(dat_zt) <- c("1", colnames(dat5))
  } else {
    dat_zt <- dat5
  }
  stopifnot(ncol(dat_zt) == n_p)

  # Marginal effect
  dat_zt_t0 = as.matrix(dplyr::mutate(dat_zt, T = 0))
  dat_zt_t1 = as.matrix(dplyr::mutate(dat_zt, T = 1))
  V1 = mean(linkinv_tanh(dat_zt_t1 %*% params2)- linkinv_tanh(dat_zt_t0 %*% params2))

  # Conditional effect
  dat_zt_c0 = as.matrix(dplyr::filter(dat_zt, T == 0))
  dat_zt_c1 = as.matrix(dplyr::filter(dat_zt, T == 1))
  V2 = mean(linkinv_tanh(dat_zt_c1 %*% params2)) - mean(linkinv_tanh(dat_zt_c0 %*% params2))

  return(c(Marg_Effect = V1, Cond_Effect = V2))
}


#' @title Association Size of Continuous Variable T
#'
#' @details The data frame \code{dat5} may contain a column with name "1" consisting of 1's
#' corresponding to the intercept term in the correlation model. If not, this column
#' will be created. Currently, only the \code{tanh} link is implemented.
#'
#' @param params a numeric vector containing the estimated correlation parameters in the CoCoA model
#' @param dat5 a data frame containing only the covariates in the correlation model (i.e., excluding X and Y)
#'
#' @return A numeric data frame containing the marginal and conditional association sizes for each
#' unique value of T, as well as the estimated parameters.
#'
#' @examples
#' T = sample(x = seq(from = 0, to = 2, by = 0.25), size = 100, replace = TRUE)
#' Z = rep(NA, times = 100)
#' for (i in 1:100){
#' Z[i] <- if (T[i] < 1){
#' # Z|T=1 ~ Bernoulli(0.3)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.3))
#' } else {
#' # Z|T=0 ~ Bernoulli(0.9)
#' as.numeric(purrr::rbernoulli(n = 1, p = 0.9))
#' }
#' }
#' dat_zt0 <- data.frame(Z = Z, T = T)
#' effect_cont_emp(params = c(0.5, 0.2, 0.5), dat5 = dat_zt0)
#'
#' @export
#'
effect_cont_emp <- function(params, dat5){
  n_p = length(params)
  params2 = matrix(params, ncol = 1)
  stopifnot("T" %in% colnames(dat5))

  # Add intercept if not in the data
  if (!("1" %in% colnames(dat5))){
    dat_zt <- cbind(1, dat5)
    colnames(dat_zt) <- c("1", colnames(dat5))
  } else {
    dat_zt <- dat5
  }

  if ("data.frame" %in% class(dat_zt)){
    dat_zt <- as.matrix(dat_zt)
  }

  # Ensure that dat_zt is in order by T
  dat_zt <- dat_zt[order(dat_zt[,"T"]),]

  # Unique values of T, range of T
  T_unique = unique(dat_zt[,"T"])
  rho_t_m = vector("numeric", length = length(T_unique))
  mean_z_m = vector("numeric", length = length(T_unique))

  # Marginal effect: how does rho change with T?
  for (i in seq_along(T_unique)){
    t_i = T_unique[i]
    dat_zt_i = dat_zt
    dat_zt_i[,"T"] <- t_i

    # What is the average rho(t), averaged over all values of Z?
    rho_t_m[i] = mean(linkinv_tanh(dat_zt_i %*% params2))
    mean_z_m[i] = ifelse(is.null(dim(dat_zt_i)),dat_zt_i["Z"] ,mean(dat_zt_i[,"Z"]))
  }
  effect_marg = data.frame(Effect = "Marginal", T = T_unique, rho_t = rho_t_m, mean_z = mean_z_m)

  # Conditional effect: how does rho change with T?
  # Windows of T based on nb_rad_n
  rho_t_c = vector("numeric", length = length(T_unique))
  mean_z_c = vector("numeric", length = length(T_unique))
  for (i in seq_along(T_unique)){
    t_i = T_unique[i]
    inds_i = which(dat_zt[,"T"] == t_i)
    dat_zt_i = dat_zt[inds_i,]

    # What is the average rho(t), averaged over all values of Z?
    rho_t_c[i] = mean(linkinv_tanh(dat_zt_i %*% params2))
    mean_z_c[i] = ifelse(is.null(dim(dat_zt_i)),dat_zt_i["Z"] ,mean(dat_zt_i[,"Z"]))
  }

  effect_cond = data.frame(Effect = "Conditional", T = T_unique, rho_t = rho_t_c, mean_z = mean_z_c)

  out = rbind(effect_marg, effect_cond)

  param_cnames = paste0(grk_letters[1:length(params)], "_hat")
  out[, param_cnames] <- as.data.frame(matrix(params, nrow = nrow(out), ncol=length(params), byrow=TRUE))

  return(out)
}
