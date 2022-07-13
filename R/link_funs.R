##' @name link_funs
##' @rdname link_funs
##'
##' @title Link and inverse link functions for the correlation model.
##'
##' @details Link functions are defined according to the C++ source code in geepack \url{https://github.com/cran/geepack/blob/master/src/famstr.cc}.
##' Both \code{linkfun_tid()} and \code{linkinv_tid()} correspond to the truncated identity link.
##'
##' @param eta a numeric value corresponding to the linear predictor, generally taking values in (-Inf, Inf).
##' @param p a numeric value corresponding to the outcome (correlation), taking values in (0,1) for \code{logit}, \code{probit}; or (-1,1) for \code{tanh}.
##'
NULL

##' @rdname link_funs
##' @examples
##' linkfun_tanh(0.5)
##' @export
linkfun_tanh <- function(p){
  log((1+p)/(1-p))
}

##' @rdname link_funs
##' @examples
##' linkinv_tanh(3)
##' @export
linkinv_tanh <- function(eta){
  (exp(eta) - 1)/(exp(eta) + 1)
}

##' @rdname link_funs
##' @export
linkfun_identity <- function(p){p}

##' @rdname link_funs
##' @export
linkinv_identity <- function(eta){eta}

##' @rdname link_funs
##' @export
linkfun_tid <- function(p){
  p <- pmax(p, 0)
  p <- pmin(p, 1)
  return(p)
}

##' @rdname link_funs
##' @export
linkinv_tid <- function(eta){
  eta <- pmax(eta, 0)
  eta <- pmin(eta, 1)
  return(eta)
}


##' @rdname link_funs
##' @export
linkfun_logit <- function(p){
  log(p/(1-p))
}

##' @rdname link_funs
##' @export
linkinv_logit <- function(eta){
  exp(eta)/(1+exp(eta))
}

##' @rdname link_funs
##' @export
linkfun_probit <- function(p){
  qnorm(p, mean = 0, sd = 1)
}
##' @rdname link_funs
##' @export
linkinv_probit <- function(eta){
  pnorm(eta, mean = 0, sd = 1)
}

##' @rdname link_funs
##' @export
linkfun_cloglog <- function(p){
  log(-log(1-p))
}

##' @rdname link_funs
##' @export
linkinv_cloglog <- function(eta){
  ans = 1 - exp(- exp(eta))
  # Ensure answer is bounded away from 0 and 1.
  ans = pmin(1 - .Machine$double.eps, ans)
  return(pmax(.Machine$double.eps, ans))
}

#' Link function label
#'
#' For plots/tables, get the name of the link function.
#'
#' @param str character of the link function
#'
#' @examples
#' ptab(cars)
#'
get_link_name <- function(str){
  str = gsub(x = str, pattern = "linkinv_", replacement = "")
  str = gsub(x = str, pattern = "linkfun_", replacement = "")

  out = dplyr::recode_factor(str,
                      "cloglog" = "Comp. Log-Log",
                      "probit" = "Probit",
                      "logit" = "Logit",
                      "fisherz" = "Tanh",
                      "identity" = "Identity")

  return(out)
}
