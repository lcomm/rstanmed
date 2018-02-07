#' Binary (M) - binary (Y) sensitivity analysis 
#'
#' @param outcomes List of named regression outcomes (named m and y)
#' @param designs List of named design matrices (named x_m, x_y, x_u)
#' @param prior List of named prior means and variance-covariance matrices
#' @param ... Additional arguments passed to \code{rstan::sampling}
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
bin_bin_sens_stan <- function(outcomes, designs, prior, ...) {
  out <- rstan::sampling(stanmodels$bin_bin_sens, 
                        data=list(N = length(outcomes$y), 
                                  P_m = ncol(designs$x_m) + 1,
                                  P_y = ncol(designs$x_y) + 1,
                                  P_u = ncol(designs$x_u),
                                  y = outcomes$y,
                                  m = outcomes$m,
                                  x_m = designs$x_m,
                                  x_y = designs$x_y,
                                  x_u = designs$x_u,
                                  prior_mean_beta = prior[["beta"]][["mean"]],
                                  prior_mean_alpha = prior[["alpha"]][["mean"]],
                                  prior_mean_gamma = prior[["gamma"]][["mean"]],
                                  prior_vcov_beta = prior[["beta"]][["vcov"]],
                                  prior_vcov_alpha = prior[["alpha"]][["vcov"]],
                                  prior_vcov_gamma = prior[["gamma"]][["vcov"]]),
                        ... )
  return(out)
}
