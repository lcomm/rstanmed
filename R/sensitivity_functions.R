
#' Calculate bias of posterior mean in Stan parameter samples
#' 
#' @param samples SxP Matrix with each column corresponding to a parameter
#' @param truth Length-P vector of true parameter values
#' @return Length-P vector of estimated bias
#' @export
calculate_bias <- function(samples, truth) {
  post_mean <- colMeans(samples)
  bias <- post_mean - truth
}

#' Evaluate Y/N credible interval coverage for Stan samples
#' 
#' @param samples SxP Matrix with each column corresponding to a parameter
#' @param truth Length-P vector of true parameter values
#' @return Length-P logical vector of interval coverage
#' @export
check_coverage <- function(samples, truth) {
  ci <- matrixStats::colQuantiles(samples, probs = c(0.025, 0.975))
  coverage <- (ci[, 1] < truth) & (ci[, 2] > truth)
  return(coverage)
}

#' Calculate parameter bias and coverage from Stan model fit
#' 
#' @param stan_fit Stan fit object
#' @param parameter_name String name for parameter extraction
#' @param truth True values of parameter
#' @return Length-2 list of length-P bias vector and length-P coverage
#' @export
get_bias_coverage <- function(stan_fit, parameter_name, truth) {
  samples <- extract(stan_fit, pars = parameter_name)[[parameter_name]]
  bias <- calculate_bias(samples, truth)
  coverage <- check_coverage(samples, truth)
  names(bias) <- names(coverage) <- paste0(parameter_name, 1:length(truth))
  return(list(bias = bias, coverage = coverage))
}

#' Return true DGP for binary-binary sensitivity simulations
#' Note: does not have A -> U!
#' @return Named list of parameter values
#' @export
return_dgp_parameters <- function() {
  list(pz1   = 0.5, 
       pz2   = 0.5, 
       gamma = c(-1, 0.2, 0.4),
       beta  = c(-2, 0.3, 0.2, 0.7, 0.8),
       alpha = c(-3, 0.3, 0.2, 1, 0.8, 0.3),
       a_params = c(-1, 0.5, 0.7))
}

#' Construct all 4 possible combinations of binary Z1 and Z2
#' @return Named data frame of the baseline covariate types
#' @export
make_bl_types <- function() {
  expand.grid(z1 = 0:1, z2 = 0:1)
}

#' Calculate population prevalance of baseline covariate types
#' 
#' @param bl_types Data frame containing the 4 unique covariate types
#' @param pz1 P(Z1 = 1)
#' @param pz2 P(Z2 = 1), assumed to be independent of Z1
#' @export
make_bl_weights <- function(bl_types, pz1, pz2) {
  stopifnot(NROW(bl_types) == 4)
  bl_weights <- rep(NA, 4)
  bl_weights[bl_types$z1 == 1 & bl_types$z2 == 1] <- pz1 * pz2
  bl_weights[bl_types$z1 == 1 & bl_types$z2 == 0] <- pz1 * (1 - pz2)
  bl_weights[bl_types$z1 == 0 & bl_types$z2 == 1] <- (1 - pz1) * pz2
  bl_weights[bl_types$z1 == 0 & bl_types$z2 == 0] <- (1 - pz1) * (1 - pz2)
  return(bl_weights)
}

#' Construct design matrix for unmeasured confounder regression model
#' Note that it does not depend on A yet!
#' 
#' @param z1 Length-N vector for baseline confounder Z1
#' @param z2 Length-N vector for baseline confounder Z2
#' @param a Length-N vector for observed or hypothetical A
#' @return NxP design matrix for the U outcome model
#' @export
make_x_u <- function(z1, z2, a) {
  cbind(1, z1, z2)
}

#' Construct design matrix for mediator regression model
#' 
#' @param z1 Length-N vector for baseline confounder Z1
#' @param z2 Length-N vector for baseline confounder Z2
#' @param a Length-N vector for observed or hypothetical A
#' @param u Length-N vector for observed (true) or hypothetical U
#' @param type Whether to make naive matrix excluding U ("observed") or 
#' full matrix that includes U ("full")
#' @return NxP design matrix for the M outcome model
#' @export
make_x_m <- function(z1, z2, a, u = NULL, type = c("observed", "full")) {
  if (type == "full") {
    return(cbind(1, z1, z2, a, u))
  } else if (type == "observed") {
    return(cbind(1, z1, z2, a))
  }
}

#' Construct design matrix for outcoe regression model
#' 
#' @param z1 Length-N vector for baseline confounder Z1
#' @param z2 Length-N vector for baseline confounder Z2
#' @param a Length-N vector for observed or hypothetical A
#' @param m Length-N vector for observed or hypothetical M
#' @param u Length-N vector for observed (true) or hypothetical U
#' @param type Whether to make naive matrix excluding U ("observed") or 
#' full matrix that includes U ("full")
#' @return NxP design matrix for the Y outcome model
#' @export
make_x_y <- function(z1, z2, a, m, u = NULL, type = c("observed", "full")) {
  if (type == "full") {
    return(cbind(1, z1, z2, a, m, u))
  } else if (type == "observed") {
    return(cbind(1, z1, z2, a, m))
  }
}

#' Simulate data set
#' 
#' @param n Number of observations to simulate
#' @param params (Optional) List of data-generating parameters
#' see \code{\link{return_dgp_parameters}} for details on structure
#' @return List containing: df (data frame with n rows),
#' outcomes (named list of m and y outcomes), and designs (named list of 
#' design matrices for u, m, and y regression models)
#' @export
simulate_data <- function(n, params = NULL) {
  if (is.null(params)) {
    params <- return_dgp_parameters()
  }
  
  # Baseline and treatment
  z1  <- rbinom(n, size = 1, prob = params$pz1)
  z2  <- rbinom(n, size = 1, prob = params$pz2)
  a   <- rbinom(n, size = 1, prob = plogis(cbind(1, z1, z2) %*% params$a_params))
  
  # Regression model data simulation
  x_u <- make_x_u(z1, z2, a)
  u   <- rbinom(n, size = 1, prob = plogis(x_u %*% params$gamma))
  x_m_full <- make_x_m(z1, z2, a, u, type = "full")
  x_m <- make_x_m(z1, z2, a, u, type = "observed")
  m   <- rbinom(n, size = 1, prob = plogis(x_m_full %*% params$beta))
  x_y_full <- make_x_y(z1, z2, a, m, u, type = "full")
  x_y <- make_x_y(z1, z2, a, m, u, type = "observed")
  y   <- rbinom(n, size = 1, prob = plogis(x_y_full %*% params$alpha))
  
  # Return
  return(list(df = data.frame(z1 = z1, z2 = z2, a = a, u = u, m = m, y = y),
              outcomes = list(y = y, m = m),
              designs = list(x_u = x_u, x_m = x_m, x_y = x_y)))
}

#' Construct prior list based on type of sensitivity analysis
#' 
#' @param alpha_mean Y regression parameter mean
#' @param beta_mean M regression parameter mean
#' @param gamma_mean U regression parameter mean
#' @param prior_type Variance type: unit variance ("unit"), informative priors on
#' non-identifiable parameters ("partial"), or strong prior information 
#' ("strict")
#' @return Named list of prior means and variance-covariance matrices
#' @export
make_prior <- function(alpha_mean, beta_mean, gamma_mean, prior_type) {
  
  # Parse
  stopifnot(prior_type %in% c("unit", "partial", "strict"))
  P_y <- length(alpha_mean)
  P_m <- length(beta_mean)
  P_u <- length(gamma_mean)
  
  # Initialize container
  prior <- list()
  prior[["alpha"]] <- list(mean = alpha_mean, vcov = NA)
  prior[["beta"]] <- list(mean = beta_mean, vcov = NA)
  prior[["gamma"]] <- list(mean = gamma_mean, vcov = NA)
  
  # Make variances
  if (prior_type == "unit") {
    prior[["alpha"]][["vcov"]] <- diag(P_y)
    prior[["beta"]][["vcov"]] <- diag(P_m)
    prior[["gamma"]][["vcov"]] <- diag(P_u)
  } else if (prior_type == "partial") {
    prior[["alpha"]][["vcov"]] <- diag(c(rep(1, P_y - 1), 0.001))
    prior[["beta"]][["vcov"]] <- diag(c(rep(1, P_m - 1), 0.001))
    prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
  } else if (prior_type == "strict") {
    prior[["alpha"]][["vcov"]] <- diag(P_y) * 0.001
    prior[["beta"]][["vcov"]] <- diag(P_m) * 0.001
    prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
  }
  
  # Return
  return(prior)
}


#' Run a single replicate of the binary-binary mediation sensitivity analysis
#' 
#' @param n Number of observations in data set
#' @param params List of data generating parameters (will be created by 
#' \code{\link{return_dgp_parameters}} if not specified)
#' @param prior_type See \code{\link{make_prior}} for details
#' @param ... Additional parameters to be passed to bin_bin_sens_stan
#' @return Stan model fit object
#' @export
run_sensitivity <- function(n, prior_type, params = NULL, ...) {
  
  # Basic checks
  stopifnot(prior_type %in% c("unit", "partial", "strict"))
  if (is.null(params)) {
    params <- return_dgp_parameters()
  }
  
  # Simulate data
  dl <- simulate_data(n = n, params = params)
  
  # Prior variances
  prior <- make_prior(alpha_mean = params$alpha, beta_mean = params$beta, 
                      gamma_mean = params$gamma, prior_type = prior_type)
  
  # Run Stan model and return fit
  sf <- bin_bin_sens_stan(dl$outcomes, dl$designs, prior, ...)
  return(sf)
  
}

#' Calculate bias and coverage from a mediation fit
#' 
#' @param stan_fit The Stan fit object (see \code{\link{run_sensitivity}})
#' Should have parameter vectors labeled alpha, beta, and gamma
#' @param params List of true parameter values, labeled alpha, beta, and gamma
#' @return List of bias and coverage vectors for all parameters
#' @export
get_parameter_bias_coverage <- function(stan_fit, params) {
  
  labs <- c(paste0("alpha[", 1:length(params$alpha), "]"),
            paste0("beta[",  1:length(params$beta), "]"),
            paste0("gamma[", 1:length(params$gamma), "]"))
  
  alpha_res <- get_bias_coverage(stan_fit, "alpha", params$alpha)
  beta_res  <- get_bias_coverage(stan_fit, "beta",  params$beta)
  gamma_res <- get_bias_coverage(stan_fit, "gamma", params$gamma)
  
  bias <- c(alpha_res$bias, beta_res$bias, gamma_res$bias)
  coverage <- c(alpha_res$coverage, beta_res$coverage, gamma_res$coverage)
  names(bias) <- names(coverage) <- labs
  
  return(list(bias = bias, coverage = coverage))
  
}

