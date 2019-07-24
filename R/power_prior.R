#' Make small data list to be used in power prior
#' 
#' @param small_dl Named list of lists containing small external data 
#' frame (named \code{df}, with column name \code{u}), outcomes (named 
#' \code{outcomes}, with vector list elements \code{y} and \code{m}), and 
#' design matrices (named \code{designs}, vector list elements \code{x_y}, 
#' \code{x_m}, \code{x_u})
#' @return Named list of collapsed outcomes (including \code{u}) and design 
#' matrices, plus a weight vector \code{w} with the number of times a data 
#' pattern is observed in the small data set
#' @export
repackage_small_dl <- function(small_dl) {
  small_dl$outcomes$u <- small_dl$df$u
  small_dl$designs$x_m <- cbind(small_dl$designs$x_m, u = small_dl$df$u)
  small_dl$designs$x_y <- cbind(small_dl$designs$x_y, u = small_dl$df$u)
  return(collapse_do_list(do_list = small_dl))
}



#' Run a single replicate of the binary-binary mediation sensitivity analysis
#' with a data power prior specification
#' 
#' @param n Number of observations in data set
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @param am_intx 0/1 flag for exposure-mediator interaction in outcome model
#' @param yu_strength Log-OR of U in outcome model
#' @param mu_strength Log-OR of U in mediator model
#' @param params List of data generating parameters for big data set (will be 
#' created by \code{\link{return_dgp_parameters}} if not specified)
#' @param small_yu_strength Log-OR of U in outcome model for small data
#' @param small_mu_strength Log-OR of U in mediator model for small data
#' @param small_params List of data generating parameters for small data set 
#' @param dd_control TODO(LCOMM)
#' @param pp_alpha Alpha parameter for power prior (default is 1, which is 
#' equivalent to data pooling)
#' (will be created by \code{\link{return_dgp_parameters}} if not specified)
#' @param result_type Whether to return full object ("raw") or only selected
#' information about NDER
#' @param n_ratio Ratio of big:small sample sizes
#' @param ... Additional parameters to be passed to bin_bin_sens_stan
#' @return Stan model fit object
#' @export
run_powerprior_replicate <- function(n, u_ei, am_intx, 
                                     yu_strength, mu_strength, params = NULL,
                                     small_yu_strength, small_mu_strength, 
                                     small_params = NULL, dd_control = NULL,
                                     pp_alpha = 1,
                                     n_ratio = 10,
                                     result_type = c("raw", "processed"),
                                     ...) {
  
  # Basic checks
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei, am_intx, yu_strength, mu_strength)
  }
  params$prior <- prior_type <- "powerprior"
  
  if (is.null(small_params)) {
    small_params <- return_dgp_parameters(u_ei, am_intx, 
                                          small_yu_strength, small_mu_strength)
  }
  
  # Simulate data
  dl <- simulate_data(n = n, params = params)
  
  #TODO(LCOMM): add prior 
  # Prior 
  #TODO(LCOMM): implement dd_control as function argument better
  if (is.null(dd_control)) {
    dd_control = list(small_n = floor(n / n_ratio),
                      params = small_params,
                      partial_vague = FALSE,
                      inflate_factor = 1)
  }
  prior <- make_prior(params = small_params, prior_type = prior_type,
                      dd_control = dd_control)
  # browser()
  # Run Stan model and return fit
  sf <- bin_bin_sens_pp_stan(dl$outcomes, dl$designs, prior, u_ei, am_intx, ...)
  
  res <- list(stan_fit = sf, params = params, data_list = dl)
  truth_nder <- calculate_nder(params = params,
                               u_ei = u_ei, am_intx = am_intx,
                               mean = TRUE)
  
  # Post-process
  gc <- extract_nder_gcomp(res$stan_fit)
  gf <- extract_nder_gform(res$stan_fit, dl$df, u_ei = u_ei, am_intx = am_intx)
  nder_gc <- gc[1]
  nder_gf <- gf[1]
  ci_gc <- gc[2:3]
  ci_gf <- gf[2:3]
  cov_gc  <- check_coverage(ci = ci_gc, truth = truth_nder)
  cov_gf  <- check_coverage(ci = ci_gf, truth = truth_nder)
  width_gc <- get_ci_width(ci_gc)
  width_gf <- get_ci_width(ci_gf)
  
  if (result_type == "raw") {
    return(list(stan_fit = sf, params = params, data_list = dl,
                truth_nder = truth_nder, 
                estimate = c(gc = nder_gc, gf = nder_gf),
                bias = c(gc = nder_gc - truth_nder,
                         gf = nder_gf - truth_nder),
                ci = c(gc = ci_gc, gf = ci_gf),
                ci_cov = c(gc = cov_gc, gf = cov_gf),
                ci_width = c(gc = width_gc, gf = width_gf)))
    
  } else {
    return(list(truth_nder = truth_nder, 
                estimate = c(gc = nder_gc, gf = nder_gf),
                bias = c(gc = nder_gc - truth_nder,
                         gf = nder_gf - truth_nder),
                ci = c(gc = ci_gc, gf = ci_gf),
                ci_cov = c(gc = cov_gc, gf = cov_gf),
                ci_width = c(gc = width_gc, gf = width_gf)))
  }
}

