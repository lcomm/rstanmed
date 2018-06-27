
#' Calculate bias of posterior mean in Stan parameter samples
#' 
#' @param estimate Length-P vector with each column corresponding to an estimate
#' @param truth Length-P vector of true parameter values
#' @return Length-P vector of estimated bias
#' @export
calculate_bias <- function(estimate, truth) {
  bias <- estimate - truth
  return(bias)
}

#' Calculate posterior mean in Stan parameter samples
#' 
#' @param samples SxP Matrix with each column corresponding to a parameter
#' @return Length-P vector of posterior means
#' @export
calculate_post_mean <- function(samples) {
  if (is.vector(samples)) {
    post_mean <- mean(samples)
  } else {
    post_mean <- colMeans(as.matrix(samples))
  }
  return(post_mean)
}

#' Make credible interval from Stan samples
#' 
#' @param samples SxP Matrix with each column corresponding to a parameter
#' @return Length-2 vector or Px2 matrix of credible interval bounds
#' @export
make_ci <- function(samples) {
  if (is.vector(samples) || (length(dim(samples)) == 1)) {
    ci <- quantile(samples, probs = c(0.025, 0.975))
  } else {
    ci <- matrixStats::colQuantiles(as.matrix(samples), probs = c(0.025, 0.975))  
  }
  names(ci) <- NULL
  return(ci)
}

#' Function to calculate width of credible intervals
#' 
#' @param ci Length-2 vector or Px2
#' @return Scalar or length-P vector of CI widths
#' @export
get_ci_width <- function(ci) {
  if (is.vector(ci)) {
    ci_width <- ci[2] - ci[1]
  } else {
    ci_width <- ci[,2] - ci[,1]
  }
  return(ci_width)
}


#' Evaluate Y/N credible interval coverage
#' 
#' @param ci Length-2 vector or Px2 matrix of CI lower and upper bounds
#' @param truth Length-P vector of true parameter values
#' @return Length-P logical vector of interval coverage
#' @export
check_coverage <- function(ci, truth) {
  if (is.vector(ci)) {
    coverage <- (ci[1] < truth) & (ci[2] > truth)
  } else {
    coverage <- (ci[, 1] < truth) & (ci[, 2] > truth)
  }
  return(coverage)
}

#' Calculate parameter bias and coverage from Stan model fit
#' 
#' @param stan_fit Stan fit object
#' @param parameter_name String name for parameter extraction
#' @param truth True values of parameter
#' @return Length-2 list of length-P bias, coverage, and CI width
#' export
# get_bias_coverage <- function(stan_fit, parameter_name, truth) {
#   if (length(parameter_name) > 1) {
#     samples <- as.matrix(extract(stan_fit, pars = parameter_name)[[parameter_name]])
#   } else {
#     samples <- extract(stan_fit, pars = parameter_name)[[parameter_name]]  
#   }
#   estimate <- calculate_post_mean(samples)
#   bias <- calculate_bias(estimate, truth)
#   ci <- make_ci(samples)
#   coverage <- check_coverage(ci, truth)
#   ci_width <- get_ci_width(ci)
#   if (is.vector(ci)) {
#     stopifnot(length(ci) == 2)
#     ci <- as.matrix(ci, ncol = 2)
#   }
#   return(list(truth = truth, estimate = estimate, bias = bias, 
#               coverage = coverage, ci_lb = ci[,1], ci_ub = ci[,2],
#               ci_width = ci_width))
# }


#' Calculate bias from a stan fit
#' 
#' @param stan_fit Stan model fit object.
#' @param pars Parameter names. See \code{\link[rstan]{extract}} for parsing details.
#' @param truth True values for parameters.
#' @return Vector of bias values.
#' @export
get_bias_stan <- function(stan_fit, pars, truth) {
  estimate <- summary(stan_fit, pars = pars)$summary[, "mean"]
  stopifnot(length(truth) == length(estimate))
  return(estimate - truth)
}

#' Calculate 95\% credible interval coverage from a stan fit
#' 
#' @param stan_fit Stan model fit object.
#' @param pars Parameter name(s). See \code{\link[rstan]{extract}} for parsing details.
#' @param truth True values for parameter(s).
#' @return Vector of logical(s) indicating whether credible interval(s) cover truth.
#' @export
get_ci_coverage_stan <- function(stan_fit, pars, truth) {
  s <- summary(stan_fit, pars = pars)$summary[, c("2.5%", "97.5%")]
  ci <- matrix(unname(s), ncol = 2)
  stopifnot(length(truth) == nrow(ci))
  covers <- (truth > ci[, 1]) & (truth < ci[, 2])
  return(covers)
}

#' Calculate 95\% credible interval coverage from a stan fit
#' 
#' @param stan_fit Stan model fit object.
#' @param pars Parameter name(s). See \code{\link[rstan]{extract}} for parsing details.
#' @param truth True values for parameter(s).
#' @return Vector of credible interval width(s).
#' @export
get_ci_width_stan <- function(stan_fit, pars) {
  s <- summary(stan_fit, pars = pars)$summary[, c("2.5%", "97.5%")]
  ci <- matrix(unname(s), ncol = 2)
  width <- ci[, 2] - ci[, 1]
  return(width)
}

#' Return true DGP for binary-binary sensitivity simulations
#' 
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @param am_intx 0/1 flag for including A*M interaction in outcome model
#' @param yu_strength Strength of U -> Y relationship on log-odds scale
#' @param mu_strength Strength of U -> M relationship on log-odds scale
#' @return Named list of parameter values
#' @export
return_dgp_parameters <- function(u_ei, am_intx, yu_strength, mu_strength) {
  stopifnot(c(u_ei, am_intx) %in% 0:1, is.numeric(c(yu_strength, mu_strength)))
  
  yu_strength <- unname(as.numeric(yu_strength))
  mu_strength <- unname(as.numeric(mu_strength))
  
  if (u_ei == 0) {
    gamma <- setNames(c(-2, 0, 0), c("(i)", "z1", "z2"))
  } else if (u_ei == 1) {
    gamma <- setNames(c(-2, 0, 0, 0.35), c("(i)", "z1", "z2", "a"))
  }
  
  if (am_intx == 0) {
    alpha <- setNames(c(-3, 0.3, 0.2, 1, 0.8, yu_strength),
                      c("(i)", "z1", "z2", "a", "m", "u"))
  } else if (am_intx == 1) {
    alpha <- setNames(c(-3, 0.3, 0.2, 1, 0.8, 0.1, yu_strength),
                      c("(i)", "z1", "z2", "a", "m", "a:m", "u"))
  }
  return(list(u_ei     = u_ei,
              am_intx  = am_intx,
              pz1      = 0.5, 
              pz2      = 0.5, 
              gamma    = gamma,
              beta     = setNames(c(-3.7, 0.3, 0.2, 0.7, mu_strength),
                                  c("(i)", "z1", "z2", "a", "u")),
              alpha    = alpha,
              a_params = setNames(c(-0.2, 0.5, 0.7),
                                  c("(i)", "z1", "z2"))))
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
#' Can depend on A or not
#' 
#' @param z1 Length-N vector for baseline confounder Z1
#' @param z2 Length-N vector for baseline confounder Z2
#' @param a Length-N vector for observed or hypothetical A
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @return NxP design matrix for the U outcome model
#' @export
make_x_u <- function(z1, z2, a, u_ei) {
  if (u_ei == 1) {
    return(cbind(`(i)` = 1, z1, z2, a))
  } else if (u_ei == 0) {
    return(cbind(`(i)` = 1, z1, z2))
  }
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
    return(cbind(`(i)` = 1, z1, z2, a, u))
  } else if (type == "observed") {
    return(cbind(`(i)` = 1, z1, z2, a))
  }
}

#' Construct design matrix for outcome regression model
#' 
#' @param z1 Length-N vector for baseline confounder Z1
#' @param z2 Length-N vector for baseline confounder Z2
#' @param a Length-N vector for observed or hypothetical A
#' @param m Length-N vector for observed or hypothetical M
#' @param u Length-N vector for observed (true) or hypothetical U
#' @param type Whether to make naive matrix excluding U ("observed") or 
#' full matrix that includes U ("full")
#' @param am_intx 0/1 Whether to include A*M interaction
#' @return NxP design matrix for the Y outcome model
#' @export
make_x_y <- function(z1, z2, a, m, u = NULL, type = c("observed", "full"),
                     am_intx) {
  base <- cbind(`(i)` = 1, z1, z2, a, m)
  if (am_intx == 1) {
    base <- cbind(base, a * m)
    colnames(base)[ncol(base)] <- "a:m"
  }
  if (type == "full") {
    return(cbind(base, u))
  } else if (type == "observed") {
    return(base)
  }
  
}

#' Simulate data set
#' 
#' @param n Number of observations to simulate
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @param am_intx 0/1 flag for A*M interaction in outcome model
#' @param yu_strength Log-OR of U in outcome model
#' @param mu_strength Log-OR of U in mediator model
#' @param params (Optional) List of data-generating parameters
#' see \code{\link{return_dgp_parameters}} for details on structure
#' @return List containing: df (data frame with n rows),
#' outcomes (named list of m and y outcomes), and designs (named list of 
#' design matrices for u, m, and y regression models)
#' @export
simulate_data <- function(n, u_ei, am_intx, yu_strength, mu_strength, params = NULL) {
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei = u_ei, am_intx = am_intx,
                                    yu_strength = yu_strength, mu_strength = mu_strength)
  }
  
  # Baseline and treatment
  # Center baseline covariates at true prevalence
  # TODO(LCOMM): figure out if centering Z1 and Z2 is the problem
  z1  <- rbinom(n, size = 1, prob = params$pz1)# - params$pz1
  z2  <- rbinom(n, size = 1, prob = params$pz2)# - params$pz2
  a   <- rbinom(n, size = 1, prob = plogis(cbind(1, z1, z2) %*% params$a_params))
  
  # Regression model data simulation
  x_u <- make_x_u(z1, z2, a, params$u_ei)
  u   <- rbinom(n, size = 1, prob = plogis(x_u %*% params$gamma))
  x_m_full <- make_x_m(z1, z2, a, u, type = "full")
  x_m <- make_x_m(z1, z2, a, u, type = "observed")
  m   <- rbinom(n, size = 1, prob = plogis(x_m_full %*% params$beta))
  x_y_full <- make_x_y(z1, z2, a, m, u, type = "full", am_intx = params$am_intx)
  x_y <- make_x_y(z1, z2, a, m, u, type = "observed", am_intx = params$am_intx)
  y   <- rbinom(n, size = 1, prob = plogis(x_y_full %*% params$alpha))
  
  # Return
  return(list(df = data.frame(z1 = z1, z2 = z2, a = a, u = u, m = m, y = y),
              outcomes = list(y = y, m = m),
              designs = list(x_u = x_u, x_m = x_m, x_y = x_y)))
}

#' Make template for prior distribution list
#' @return Empty list with correct naming structure
make_prior_shell <- function() {
  prior <- list()
  prior[["gamma"]] <- prior[["beta"]] <- prior[["alpha"]] <- 
    list(mean = NA, vcov = NA)
  return(prior)
}

#' Make binary-binary data driven prior based on simulated or real external data
#' 
#' Make a list of prior mean and variance-covariance matrices using an external
#' data source. Assumes binary mediator and binary outcome (for now).
#' 
#' @param small_data Data list with data set on which to base prior information 
#' (see \code{\link{simulate_data}} for format)
#' @param partial_vague Whether or not to inflate prior variance for more 
#' identified parameters
#' @param inflate_factor Variance inflation factor (only used if 
#' partial_vague = TRUE)
#' @return Named list with prior means and variances
#' @export
make_dd_prior <- function(small_data, partial_vague, inflate_factor = NULL) {
  
  #Fit maximum likelihood models to small data set
  fit_m <- glm(small_data$df$m ~ -1 + small_data$designs$x_m + small_data$df$u, 
               family = binomial(link = "logit"))
  fit_y <- glm(small_data$df$y ~ -1 + small_data$designs$x_y + small_data$df$u, 
               family = binomial(link = "logit"))
  fit_u <- glm(small_data$df$u ~ -1 + small_data$designs$x_u, 
               family = binomial(link = "logit"))
  
  # Set prior means and extract covariance matrices
  prior <- make_prior_shell()
  prior[["alpha"]][["mean"]] <- unname(coef(fit_y))
  prior[["beta"]][["mean"]]  <- unname(coef(fit_m))
  prior[["gamma"]][["mean"]] <- unname(coef(fit_u))
  a_vcov <- unname(vcov(fit_y))
  b_vcov <- unname(vcov(fit_m))
  g_vcov <- unname(vcov(fit_u))
  
  if (partial_vague) {
    # Inflate all variances and covariances by inflate_factor,
    # EXCEPT u-related parameters, which are all in last column
    # (and obviously gamma is u-related)
    sd_inflate <- sqrt(inflate_factor)
    a_vcov <- a_vcov * inflate_factor
    a_vcov[, ncol(a_vcov)] <- a_vcov[, ncol(a_vcov)] / sd_inflate
    a_vcov[nrow(a_vcov), ] <- a_vcov[nrow(a_vcov), ] / sd_inflate
    
    b_vcov <- b_vcov * inflate_factor
    b_vcov[, ncol(b_vcov)] <- b_vcov[, ncol(b_vcov)] / sd_inflate
    b_vcov[nrow(b_vcov), ] <- b_vcov[nrow(b_vcov), ] / sd_inflate
  }
  
  # Set prior covariance matrices
  prior[["alpha"]][["vcov"]] <- a_vcov
  prior[["beta"]][["vcov"]]  <- b_vcov
  prior[["gamma"]][["vcov"]] <- g_vcov
  
  # Apply names
  names(prior[["alpha"]][["mean"]]) <- colnames(small_data$designs$x_y)
  names(prior[["beta"]][["mean"]])  <- colnames(small_data$designs$x_m)
  names(prior[["gamma"]][["mean"]]) <- colnames(small_data$designs$x_u)
  
  return(prior)
}

#' Construct prior list based on type of sensitivity analysis
#' 
#' @param params List of parameters
#' @param prior_type Variance type: unit variance ("unit"), informative priors 
#' on non-identifiable parameters ("partial"), strong prior information 
#' ("strict"), or data-driven ("dd") based on external or simulated data set
#' @param dd_control Named list of data-driven options. If external data is to
#' be used, "small_data" needs to be formatted as output from 
#' \code{\link{simulate_data}}. If secondary data source is to be simulated, 
#' need a list of parameters as from \code{\link{return_dgp_parameters}}. 
#' @return Named list of prior means and variance-covariance matrices
#' @export
make_prior <- function(params, prior_type = c("unit", "partial", "strict", "dd"),
                       dd_control = NULL) {

  stopifnot(prior_type %in% c("unit", "partial", "strict", "dd"))
  
  if (prior_type %in% c("unit", "partial", "strict")) {
    if (is.null(params)) {
      stop("Need to supply underlying parameters for this prior_type")
    }
    P_y <- length(params$alpha)
    P_m <- length(params$beta)
    P_u <- length(params$gamma)
    
    # Initialize container
    prior <- make_prior_shell()
    
    # Set prior means to be truth
    prior[["alpha"]][["mean"]] <- params$alpha
    prior[["beta"]][["mean"]]  <- params$beta
    prior[["gamma"]][["mean"]] <- params$gamma
    
    # Set prior variances according to specification
    if (prior_type == "unit") {
      prior[["alpha"]][["vcov"]] <- diag(P_y)
      prior[["beta"]][["vcov"]]  <- diag(P_m)
      prior[["gamma"]][["vcov"]] <- diag(P_u)
    } else if (prior_type == "partial") {
      aa <- ifelse(params$am_intx == 1, 2, 1)
      prior[["alpha"]][["vcov"]] <- diag(c(rep(1, P_y - aa), rep(0.001, aa)))
      prior[["beta"]][["vcov"]]  <- diag(c(rep(1, P_m - 1), 0.001))
      prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
    } else if (prior_type == "strict") {
      prior[["alpha"]][["vcov"]] <- diag(P_y) * 0.001
      prior[["beta"]][["vcov"]]  <- diag(P_m) * 0.001
      prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
    }
    
    # Apply names
    names(prior[["alpha"]][["mean"]]) <- names(params$alpha)
    names(prior[["beta"]][["mean"]])  <- names(params$beta)
    names(prior[["gamma"]][["mean"]]) <- names(params$gamma)
    
  } else if (prior_type == "dd") {
    
    if (is.null(dd_control)) {
      stop("Need to specify input data set or simulation parameters for dd")
    } else {
      # Basic option checks
      valid_options <- c("params", "small_n", "small_data", 
                         "partial_vague", "inflate_factor")
      opts_given <- names(dd_control)
      stopifnot(opts_given %in% valid_options)
      if (c("params", "small_n") %in% opts_given &&
          !("small_data" %in% opts_given)) {
        dd_control$small_data <- simulate_data(n = dd_control$small_n, 
                                               params = dd_control$params)
      } else if (!("small_data" %in% opts_given)) {
        stop("Need to specify input data set or simulation parameters for dd")
      } else {
        stop("Can only specify small_data OR params and small_n")
      }
      
      # Defaults
      if (is.null(dd_control$partial_vague)) {
        dd_control$partial_vague <- TRUE
      } 
      if (is.null(dd_control$inflate_factor)) {
        dd_control$inflate_factor <- 100
      }
      
      prior <- make_dd_prior(small_data = dd_control$small_data, 
                             partial_vague = dd_control$partial_vague,
                             inflate_factor = dd_control$inflate_factor)
    }
  }
  
  # Return
  return(prior)
}


#' Run a single replicate of the binary-binary mediation sensitivity analysis
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
#' (will be created by \code{\link{return_dgp_parameters}} if not specified)
#' @param prior_type See \code{\link{make_prior}} for details
#' @param result_type Whether to return full object ("raw") or only bias,
#' coverage and width of NDER
#' @param ... Additional parameters to be passed to bin_bin_sens_stan
#' @return Stan model fit object
#' @export
run_bdf_replicate <- function(n, prior_type, u_ei, am_intx, 
                              yu_strength, mu_strength, params = NULL,
                              small_yu_strength, small_mu_strength, 
                              small_params = NULL, dd_control = NULL,
                              result_type = c("raw", "processed"),
                              ...) {
  
  # Basic checks
  stopifnot(prior_type %in% c("unit", "partial", "strict", "dd"))
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei, am_intx, yu_strength, mu_strength)
  }
  params$prior = prior_type
  
  if (is.null(small_params)) {
    small_params <- return_dgp_parameters(u_ei, am_intx, 
                                          small_yu_strength, small_mu_strength)
  }
  
  # Simulate data
  dl <- simulate_data(n = n, params = params)
  
  # Prior 
  #TODO(LCOMM): implement dd_control as function argument better
  #TODO(LCOMM): default to prior variance inflation later
  if (is.null(dd_control)) {
    dd_control = list(small_n = floor(n/10), 
                      params = small_params,
                      partial_vague = TRUE,
                      inflate_factor = 1)
  } else {
    if (!("small_n" %in% names(dd_control))) {
      dd_control$small_n <- floor(n/10)
    }
    if (!("params" %in% names(dd_control))) {
      dd_control$params <- small_params
    }
    if (!("partial_vague" %in% names(dd_control))) {
      dd_control$partial_vague <- TRUE
    }
    if (!("inflate_factor" %in% names(dd_control))) {
      dd_control$inflate_factor <- 10
    }
  }
  prior <- make_prior(params = small_params, prior_type = "dd",
                      dd_control = dd_control)
  
  # Run Stan model and return fit
  sf <- bin_bin_sens_stan(dl$outcomes, dl$designs, prior, u_ei, am_intx, ...)
  res <- list(stan_fit = sf, params = params, data_list = dl)
  truth_nder <- calculate_nder(u_ei = params$u_ei,
                               am_intx = params$am_intx,
                               yu_strength = yu_strength,
                               mu_strength = mu_strength,
                               mean = TRUE)
  estimate1_nder <- summary(res$stan_fit, pars = "meffects[1]")$summary[, "mean"]
  v2_KxR <- calculate_nder_stan(res$stan_fit, params$u_ei, params$am_intx)
  v2_samples <- colMeans(v2_KxR)
  v2_bootweights <- get_bayesian_bootweights(dl$df, R = length(v2_samples))
  v2_ci <- make_ci(v2_samples)
  v2_ci_boot <- make_ci(colSums(v2_bootweights * v2_KxR))
  estimate2_nder <- mean(v2_samples)
  bias1 <- calculate_bias(estimate1_nder, truth_nder)
  bias2 <- calculate_bias(estimate2_nder, truth_nder)
  coverage1 <- get_ci_coverage_stan(res$stan_fit, "meffects[1]", truth_nder)
  coverage2 <- check_coverage(v2_ci, truth_nder)
  coverage2_boot <- check_coverage(v2_ci_boot, truth_nder)
  width1 <- get_ci_width_stan(res$stan_fit, "meffects[1]")
  width2 <- v2_ci[2] - v2_ci[1]
  width2_boot <- v2_ci_boot[2] - v2_ci_boot[1]
  
  if (result_type == "raw") {
    return(list(stan_fit = sf, params = params, data_list = dl,
                truth_nder = truth_nder, 
                estimate1_nder = estimate1_nder,
                estimate2_nder = estimate2_nder,
                bias1 = bias1, bias2 = bias2,
                coverage1 = coverage1, coverage2 = coverage2, 
                coverage2_boot = coverage2_boot,
                width1 = width1, width2 = width2,
                width2_boot = width2_boot))
  } else {
    return(list(truth_nder = truth_nder, 
                estimate1_nder = estimate1_nder,
                estimate2_nder = estimate2_nder,
                bias1 = bias1, bias2 = bias2,
                coverage1 = coverage1, coverage2 = coverage2, 
                coverage2_boot = coverage2_boot,
                width1 = width1, width2 = width2,
                width2_boot = width2_boot))
  }
  
}

#' Calculate parameter bias and coverage from a mediation fit
#' 
#' @param stan_fit The Stan fit object (see \code{\link{run_bdf_replicate}})
#' Should have parameter vectors labeled alpha, beta, and gamma
#' @param params List of true parameter values, labeled alpha, beta, and gamma
#' @return List of bias and coverage vectors for all parameters
#' @export
get_parameter_bias_coverage <- function(stan_fit, params) {
  
  labels <- c(paste0("alpha[", 1:length(params$alpha), "]"),
              paste0("beta[",  1:length(params$beta), "]"),
              paste0("gamma[", 1:length(params$gamma), "]"))
  
  alpha_res <- get_bias_coverage(stan_fit, "alpha", params$alpha)
  beta_res  <- get_bias_coverage(stan_fit, "beta",  params$beta)
  gamma_res <- get_bias_coverage(stan_fit, "gamma", params$gamma)
  
  bias <- c(alpha_res$bias, beta_res$bias, gamma_res$bias)
  coverage <- c(alpha_res$coverage, beta_res$coverage, gamma_res$coverage)
  ci_width <- c(alpha_res$ci_width, beta_res$ci_width, gamma_res$ci_width)
  names(bias) <- names(coverage) <- names(ci_width) <- labels
  
  return(list(bias = bias, coverage = coverage, ci_width = ci_width))
  
}


#' Calculate randomized interventional analog to natural direct effect
#'
#' @param params List of parameters (optional)
#' @param alpha Vector of Y regression coefficients (optional) 
#' @param beta Vector of M regression coefficients (optional) 
#' @param gamma Vector of U regression coefficients (optional)
#' @param u_ei Whether U is exposure-induced (optional)
#' @param am_intx Whether exposure-M interaction in outcome (optional)
#' @param yu_strength U coefficient in Y regression (optional)
#' @param mu_strength U coefficient in M regression (optional)
#' @param z1 Vector containing levels of z1 covariate to evaluate conditional NDER
#' @param z2 Vector containing levels of z2 covariate to evaluate conditional NDER
#' @param mean Whether to average (defaults)
#' @param weights Vector of weights for z1 and z2 combinations (optional). Default is to 
#' give every combination equal weight
#' @return K x 1 Matrix of conditional NDERs, for each of K covariate patterns
#' @export
calculate_nder <- function(params = NULL, 
                           alpha = NULL, beta = NULL, gamma = NULL, 
                           u_ei = NULL, am_intx = NULL, 
                           yu_strength = NULL, mu_strength = NULL,
                           z1 = 0:1, z2 = 0:1, mean = FALSE, weights = NULL) {
  
  if ((is.null(u_ei) | is.null(am_intx) | is.null(yu_strength) | 
       is.null(mu_strength)) & (is.null(params))) {
    stop("Need to specify u_ei, am_intx, yu_strength, and mu_strength if not 
         supplying params")
  }
  
  if (is.null(params)) {
    params  <- return_dgp_parameters(u_ei, am_intx, yu_strength, mu_strength)
    alpha   <- params$alpha
    beta    <- params$beta
    gamma   <- params$gamma
    u_ei    <- params$u_ei
    am_intx <- params$am_intx
  }
  
  # Make all possible combinations for provided z1, z2
  df <- expand.grid(z1 = z1, z2 = z2, u_for_y = 0:1, m = 0:1)
  bl_types <- make_bl_types()
  bl_types$bl_type <- 1:NROW(bl_types)
  if (mean == TRUE & is.null(weights)) {
    bl_types$weights <- 1/NROW(bl_types)
  } else if (!is.null(weights)) {
    stopifnot(length(weights) == NROW(bl_types))
    bl_types$weights <- weights
  }
  df <- merge(df, bl_types, by = c("z1", "z2"))
  
  # for a = 0 world
  pu1_a0 <- plogis(make_x_u(df$z1, df$z2, a = 0, u_ei) %*% gamma)
  pu_for_y_a0 <- ifelse(df$u_for_y == 1, pu1_a0, 1 - pu1_a0)
  
  # for randomized intervention on m (only matters in a = 0 world)
  pm1_a0u0 <- plogis(make_x_m(df$z1, df$z2, a = 0, u = 0, "full") %*% beta)
  pm1_a0u1 <- plogis(make_x_m(df$z1, df$z2, a = 0, u = 1, "full") %*% beta)
  pm1_a0   <- pm1_a0u1 * pu1_a0 + pm1_a0u0 * (1 - pu1_a0)
  pm_a0    <- ifelse(df$m == 1, pm1_a0, 1 - pm1_a0)
  
  # for a = 1 world
  pu1_a1 <- plogis(make_x_u(df$z1, df$z2, a = 1, u_ei) %*% gamma)
  pu_for_y_a1 <- ifelse(df$u_for_y == 1, pu1_a1, 1 - pu1_a1)
  
  # mean y for covariate/m/u combinations
  ey_cf1 <- plogis(make_x_y(df$z1, df$z2, a = 1, df$m, df$u_for_y, "full", am_intx) 
                   %*% alpha)
  ey_cf2 <- plogis(make_x_y(df$z1, df$z2, a = 0, df$m, df$u_for_y, "full", am_intx)
                   %*% alpha)
  
  # CF1 : E[Y(a = 1, u = u(a = 1), m = g(a = 0))]
  w1 <- pm_a0 * pu_for_y_a1
  ECF1 <- sapply(bl_types$bl_type, function(type) {
    weighted.mean(ey_cf1[df$bl_type == type], w = w1[df$bl_type == type])
  })
  
  # CF2 : E[Y(a = 0, u = u(a = 0), m = g(a = 0))]
  w2 <- pm_a0 * pu_for_y_a0
  ECF2 <- sapply(bl_types$bl_type, function(type) {
    weighted.mean(ey_cf2[df$bl_type == type], w = w2[df$bl_type == type])
  })
  
  # Return
  if (mean == TRUE) {
    return(weighted.mean(ECF1 - ECF2, w = bl_types$weights))
  } else {
    return(ECF1 - ECF2)
  }
}

#' Calculate NDER from Stan sensitivity analysis fit
#' Not simulation-based!
#' Returns conditional on specific z1, z2 values
#' 
#' @param stan_fit Result list with element stan_fit containing the stan fit
#' @param u_ei Whether unmeasured confounder is exposure-induced
#' @param am_intx Whether outcome model includes exposure-mediator interaction
#' @param \dots Parameters to pass to \code{\link{calculate_nder}}
#' @return K x R matrix of conditional NDER values for the K covariate patterns requested
#' from R MCMC iterations (excludes warmup)
#' @export
calculate_nder_stan <- function(stan_fit, u_ei, am_intx, ...) {
  as <- extract(stan_fit, pars = "alpha")[["alpha"]]
  bs <- extract(stan_fit, pars = "beta")[["beta"]]
  gs <- extract(stan_fit, pars = "gamma")[["gamma"]]
  niter <- nrow(as)
  ncol_alpha <- ncol(as)
  ncol_beta <- ncol(bs)
  
  return(sapply(1:niter, function(r) {
    calculate_nder(params = NULL, 
                   alpha = as[r,], beta = bs[r,], gamma = gs[r,], 
                   u_ei = u_ei, am_intx = am_intx,
                   yu_strength = as[r, ncol_alpha], 
                   mu_strength = bs[r, ncol_beta],
                   ...)
  }))
}

#' Make matrix of covariate pattern weights for Bayesian bootstrap
#' 
#' @param dat Data frame containing z1 and z2
#' @param R Number of weight replicates to draw
#' @return K x R matrix with dirichlet weights in columns (K = 4 right now)
get_bayesian_bootweights <- function(dat, R = 1) {
  poss_patterns <- make_bl_types()
  dat$.one <- 1
  obs_patterns <- aggregate(.one ~ z1 + z2, FUN = NROW, data = dat)
  counts <- merge(poss_patterns, obs_patterns, all.x = TRUE)$.one
  counts[is.na(counts)] <- 0
  bootweights <- replicate(R, sapply(counts, FUN = function(x) 
                           rgamma(n = 1, shape = x, scale = 1)))
  return(bootweights %*% diag(1 / colSums(bootweights)))
}
