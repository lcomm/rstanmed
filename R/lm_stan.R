#' Bayesian linear regression with Stan
#'
#' @param x vector of inputs
#' @param y vector of outputs
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
lm_stan <- function(x, y) {
  out <- rstan::sampling(stanmodels$lm, data=list(x=x, y=y, N=length(y)))
  return(out)
}
