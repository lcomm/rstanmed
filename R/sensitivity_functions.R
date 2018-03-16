
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
#' 
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @return Named list of parameter values
#' @export
return_dgp_parameters <- function(u_ei) {
  if (u_ei == 0) {
    gamma <- setNames(c(-1, 0.2, 0.4), c("(i)", "z1", "z2"))
  } else if (u_ei == 1) {
    gamma <- setNames(c(-1, 0.2, 0.4, 0.5), c("(i)", "z1", "z2", "a"))
  }
  return(list(u_ei     = u_ei,
              pz1      = 0.5, 
              pz2      = 0.5, 
              gamma    = gamma,
              beta     = setNames(c(-2, 0.3, 0.2, 0.7, 0.8),
                                  c("(i)", "z1", "z2", "a", "u")),
              alpha    = setNames(c(-3, 0.3, 0.2, 1, 0.8, 0.3),
                                  c("(i)", "z1", "z2", "a", "m", "u")),
              a_params = setNames(c(-1, 0.5, 0.7),
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
    return(cbind(1, z1, z2, a))
  } else if (u_ei == 0) {
    return(cbind(1, z1, z2))
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
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @param params (Optional) List of data-generating parameters
#' see \code{\link{return_dgp_parameters}} for details on structure
#' @return List containing: df (data frame with n rows),
#' outcomes (named list of m and y outcomes), and designs (named list of 
#' design matrices for u, m, and y regression models)
#' @export
simulate_data <- function(n, u_ei, params = NULL) {
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei)
  }
  
  # Baseline and treatment
  z1  <- rbinom(n, size = 1, prob = params$pz1)
  z2  <- rbinom(n, size = 1, prob = params$pz2)
  a   <- rbinom(n, size = 1, prob = plogis(cbind(1, z1, z2) %*% params$a_params))
  
  # Regression model data simulation
  x_u <- make_x_u(z1, z2, a, params$u_ei)
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
#' @param params List of parameters
#' @param prior_type Variance type: unit variance ("unit"), informative priors on
#' non-identifiable parameters ("partial"), strong prior information 
#' ("strict"), or data-driven ("dd strict") based on MLE from a simulated data set
#' @param n Size of data set (only needed if data-driven)
#' @param n_frac Ratio of small n to analysis data n (only needed if data-driven)
#' @return Named list of prior means and variance-covariance matrices
#' @export
make_prior <- function(params, prior_type, n = NULL, n_frac = 0.1) {
  
  # Parse
  stopifnot(prior_type %in% c("unit", "partial", "strict", "dd"))
  if (prior_type %in% c("unit", "partial", "strict")) {
    P_y <- length(params$alpha)
    P_m <- length(params$beta)
    P_u <- length(params$gamma)
  } else if (prior_type == "dd") {
    
    if (is.null(n) | is.null(n_frac)) {
      stop("Need to specify sample size fraction for simulated external data")
    }
    
    s_dl <- simulate_data(n = floor(n * n_frac), params = params)
    fit_m <- glm(s_dl[["df"]][["m"]] ~ -1 + s_dl[["designs"]][["x_m"]] + s_dl[["df"]][["u"]], 
                 family = binomial(link = "logit"))
    fit_y <- glm(s_dl[["df"]][["y"]] ~ -1 + s_dl[["designs"]][["x_y"]] + s_dl[["df"]][["u"]], 
                 family = binomial(link = "logit"))
    fit_u <- glm(s_dl[["df"]][["u"]] ~ -1 + s_dl[["designs"]][["x_u"]], 
                 family = binomial(link = "logit"))
  }
  
  
  # Initialize container
  prior <- list()
  prior[["gamma"]] <- prior[["beta"]] <- prior[["alpha"]] <- list(mean = NA, vcov = NA)
  
  # Set prior means
  if (prior_type %in% c("unit", "partial", "strict")) {
    prior[["alpha"]][["mean"]] <- params$alpha
    prior[["beta"]][["mean"]]  <- params$beta
    prior[["gamma"]][["mean"]] <- params$gamma
  } else if (prior_type == "dd") {
    prior[["alpha"]][["mean"]] <- unname(coef(fit_y))
    prior[["beta"]][["mean"]]  <- unname(coef(fit_m))
    prior[["gamma"]][["mean"]] <- unname(coef(fit_u))
  }
  
  names(prior[["alpha"]][["mean"]]) <- names(params$alpha)
  names(prior[["beta"]][["mean"]])  <- names(params$beta)
  names(prior[["gamma"]][["mean"]]) <- names(params$gamma)
  
  # Set prior variances
  if (prior_type == "unit") {
    prior[["alpha"]][["vcov"]] <- diag(P_y)
    prior[["beta"]][["vcov"]]  <- diag(P_m)
    prior[["gamma"]][["vcov"]] <- diag(P_u)
  } else if (prior_type == "partial") {
    prior[["alpha"]][["vcov"]] <- diag(c(rep(1, P_y - 1), 0.001))
    prior[["beta"]][["vcov"]]  <- diag(c(rep(1, P_m - 1), 0.001))
    prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
  } else if (prior_type == "strict") {
    prior[["alpha"]][["vcov"]] <- diag(P_y) * 0.001
    prior[["beta"]][["vcov"]]  <- diag(P_m) * 0.001
    prior[["gamma"]][["vcov"]] <- diag(P_u) * 0.001
  } else if (prior_type == "dd") {
    prior[["alpha"]][["vcov"]] <- unname(vcov(fit_y))
    prior[["beta"]][["vcov"]]  <- unname(vcov(fit_m))
    prior[["gamma"]][["vcov"]] <- unname(vcov(fit_u))
  }
  
  # Return
  return(prior)
}


#' Run a single replicate of the binary-binary mediation sensitivity analysis
#' 
#' @param n Number of observations in data set
#' @param u_ei 0/1 flag for whether U should be exposure-induced
#' @param params List of data generating parameters (will be created by 
#' \code{\link{return_dgp_parameters}} if not specified)
#' @param prior_type See \code{\link{make_prior}} for details
#' @param ... Additional parameters to be passed to bin_bin_sens_stan
#' @return Stan model fit object
#' @export
run_sensitivity <- function(n, prior_type, u_ei, params = NULL, ...) {
  
  # Basic checks
  stopifnot(prior_type %in% c("unit", "partial", "strict", "dd"))
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei)
  }
  params$prior = prior_type
  
  # Simulate data
  dl <- simulate_data(n = n, params = params)
  
  # Prior 
  prior <- make_prior(params, prior_type, n, n_frac = 0.1)
  
  # Run Stan model and return fit
  sf <- bin_bin_sens_stan(dl$outcomes, dl$designs, prior, ...)
  return(list(stan_fit = sf, params = params, data_list = dl))
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


#' Calculate randomized interventional analog to natural direct effect
#'
#' @param params List of parameters (optional)
#' @param alpha Vector of Y regression coefficients (optional) 
#' @param beta Vector of M regression coefficients (optional) 
#' @param gamma Vector of U regression coefficients (optional)
#' @param u_ei Whether U is exposure-induced (optional)
#' @param z1 Vector containing levels of z1 covariate to evaluate conditional NDER
#' @param z2 Vector containing levels of z2 covariate to evaluate conditional NDER
#' @param mean Whether to average (defaults)
#' @param weights Vector of weights for z1 and z2 combinations (optional). Default is to 
#' give every combination equal weight
#' @return K x 1 Matrix of conditional NDERs, for each of K covariate patterns
#' @export
calculate_nder <- function(params = NULL, alpha = NULL, beta = NULL, gamma = NULL, 
                           u_ei = NULL, z1 = 0:1, z2 = 0:1, mean = FALSE, weights = NULL) {
  if (is.null(u_ei) & is.null(params)) {
    stop("Need to specify u_ei if not supplying params")
  }
    
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei)
    alpha <- params$alpha
    beta  <- params$beta
    gamma <- params$gamma
    u_ei  <- params$u_ei
  }
  
  # Make all possible combinations for provided z1, z2
  df <- expand.grid(z1 = z1, z2 = z2, u_for_y = 0:1, m = 0:1)
  bl_types <- make_bl_types()
  bl_types$bl_type <- 1:NROW(bl_types)
  if (mean == TRUE & is.null(weights)) {
    bl_types$weights <- 1/NROW(bl_types)
  } else {
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
  ey_cf1 <- plogis(make_x_y(df$z1, df$z2, a = 1, df$m, df$u_for_y, "full") 
                   %*% alpha)
  ey_cf2 <- plogis(make_x_y(df$z1, df$z2, a = 0, df$m, df$u_for_y, "full") 
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
#' Returns conditional on specific z1, z2 values
#' 
#' @param sens_res Result list with element stan_fit containing the stan fit
#' @param \dots Parameters to pass to \code{\link{calculate_nder}}
#' @return K x R matrix of conditional NDER values for the K covariate patterns requested
#' from R MCMC iterations (excludes warmup)
#' @export
calculate_nder_stan <- function(sens_res, ...) {
  as <- extract(sens_res$stan_fit, pars = "alpha")[["alpha"]]
  bs <- extract(sens_res$stan_fit, pars = "beta")[["beta"]]
  gs <- extract(sens_res$stan_fit, pars = "gamma")[["gamma"]]
  niter <- ncol(as)
  
  return(sapply(1:niter, function(r) {
    calculate_nder(params = NULL, 
                   alpha = as[r,], beta = bs[r,], gamma = gs[r,], 
                   ...)
  }))
}

