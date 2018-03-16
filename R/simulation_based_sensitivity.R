#' Function to do simulation-based (g-computation) on stan sensitivity fit
#' 
#' @param sens_res Result from \code{\link{run_sensitivity()}}
#' @param mean_only 0/1 flag for whether to use only mean of final counterfactual
#' @param seed Seed for random number generation
#' @return Length-B vector of marginal NDER based on simulated counterfactuals
#' @export
calculate_simnder_stan <- function(sens_res, mean_only = 1) {
  dl <- sens_res$data_list
  u_ei <- sens_res$params$u_ei
  as <- extract(sens_res$stan_fit, par = "alpha")$alpha
  bs <- extract(sens_res$stan_fit, par = "beta")$beta
  gs <- extract(sens_res$stan_fit, par = "gamma")$gamma
  # TODO(LCOMM): Fix this to expose correctly
  nders <- rstanmed:::sim_nder_rng(alpha = as, beta = bs, gamma = gs, u_ei = u_ei,
                        x_y = dl$designs$x_y, 
                        x_m = dl$designs$x_m, 
                        x_u = dl$designs$x_u,
                        mean_only = mean_only, seed = get_seed(sens_res$stan_fit))
  mnders <- colMeans(nders)
  return(mnders)
}

