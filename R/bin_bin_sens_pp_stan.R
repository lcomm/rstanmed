#' Binary (M) - binary (Y) sensitivity analysis with power prior
#'
#' (Does collapsing of design matrices internally)
#' In power prior, \code{prior[[parameter]][[mean]]} and 
#' \code{prior[[parameter]][[vcov]]} are NOT used to construct a prior, but just
#' to have reasonable default scaling and centering of parameters for 
#' sampling and initialization.
#'
#' @param outcomes List of named regression outcomes (named m and y)
#' @param designs List of named design matrices (named x_m, x_y, x_u) - 
#' with exposure named a
#' @param prior List of named prior means and variance-covariance matrices
#' (only used for, not for setting prior), and a 
#' @param small_outcomes List of named regression outcomes (named m and y)
#' @param small_designs List of named design matrices (named x_m, x_y, x_u) - 
#' with exposure named a, from external data
#' @param pp_alpha Alpha tuning parameter for power prior 
#' @param u_ei 0/1 Whether u is assumed to be exposure-induced
#' @param am_intx 0/1 Whether exposure-mediation interaction is included in 
#' outcome model
#' @param mean_only 0/1 Whether to sample outcomes or use means for mediation qtys
#' @param ... Additional arguments passed to \code{rstan::sampling}
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
bin_bin_sens_pp_stan <- function(outcomes, designs,
                                 prior,
                                 pp_alpha,
                                 u_ei = NULL,
                                 am_intx = NULL,
                                 mean_only = 1,
                                 ...) {
  
  # Exposure is last column of mediator design
  exposure <- designs$x_m[, ncol(designs$x_m)]
  
  # Infer exposure-induced if exposure in last column of unmeasured design
  last_u <- designs$x_u[, ncol(designs$x_u)]
  last_u_eq_exp <- all(last_u == exposure) * 1
  if (is.null(u_ei)) {
    u_ei <- last_u_eq_exp
  } else {
    if ((u_ei == 1) && (last_u_eq_exp == 0)) {
      stop("u_ei = 1 but exposure is not in final column of designs$x_u")
    } else if ((u_ei == 0) && (last_u_eq_exp == 1)) {
      stop("u_ei = 0 but exposure is in final column of designs$x_u")
    }
  }
  
  # Infer exposure-mediator interaction in outcome model
  intx <- exposure * outcomes$m
  last_y <- designs$x_y[, ncol(designs$x_y)]
  second_last_y <- designs$x_y[, ncol(designs$x_y) - 1]
  if (is.null(am_intx)) {
    am_intx <- all(last_y == intx) * 1
  }
  if (am_intx == 0) {
    if (any(last_y != outcomes$m)) {
      stop("am_intx = 0 but final column of designs$x_y is not mediator")  
    } else if (any(second_last_y != outcomes$a)) {
      stop("am_intx = 0 but second from final column of designs$x_y is not exposure")
    }
  } else if (am_intx == 1) {
    third_last_y <- designs$x_y[, ncol(designs$x_y) - 2]
    if (any(last_y != outcomes$intx)) {
      stop("am_intx = 1 but final column of designs$x_y is not interaction")  
    } else if (any(second_last_y != outcomes$m)) {
      stop("am_intx = 1 but second from final column of designs$x_y is not mediator")
    } else if (any(third_last_y != outcomes$a)) {
      stop("am_intx = 1 but third from final column of designs$x_y is not exposure")
    }
  }
  
  # Collapsed version of main data set design matrices and outcomes
  cdo <- collapse_do_list(do_list = list(designs = designs, 
                                         outcomes = outcomes))

  out <- rstan::sampling(stanmodels$bin_bin_sens_pp, 
                         data=list(N = length(cdo$outcomes$y), 
                                   P_m = ncol(cdo$designs$x_m) + 1,
                                   P_y = ncol(cdo$designs$x_y) + 1,
                                   P_u = ncol(cdo$designs$x_u),
                                   y = cdo$outcomes$y,
                                   m = cdo$outcomes$m,
                                   x_m = cdo$designs$x_m,
                                   x_y = cdo$designs$x_y,
                                   x_u = cdo$designs$x_u,
                                   location_beta = prior[["beta"]][["mean"]],
                                   location_alpha = prior[["alpha"]][["mean"]],
                                   location_gamma = prior[["gamma"]][["mean"]],
                                   scale_V_beta = prior[["beta"]][["vcov"]],
                                   scale_V_alpha = prior[["alpha"]][["vcov"]],
                                   scale_V_gamma = prior[["gamma"]][["vcov"]],
                                   N_s = length(prior$small_cdo$outcomes$y),
                                   small_y = prior$small_cdo$outcomes$y,
                                   small_m = prior$small_cdo$outcomes$m,
                                   small_u = prior$small_cdo$outcomes$u,
                                   small_x_m = prior$small_cdo$designs$x_m,
                                   small_x_y = prior$small_cdo$designs$x_y,
                                   small_x_u = prior$small_cdo$designs$x_u,
                                   u_ei = u_ei,
                                   am_intx = am_intx,
                                   mean_only = mean_only,
                                   w = cdo$w,
                                   small_w = prior$small_cdo$w,
                                   pp_alpha = pp_alpha),
                        ...)
  return(out)
}
