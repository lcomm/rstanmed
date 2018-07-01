#' Take a bootstrap sample
#' 
#' @param df Data frame of n rows
#' @return Data frame of n rows sampled with replacement from df
#' @export
boot_samp <- function(df) {
  n <- NROW(df)
  return(df[sample(1:n, size = n, replace = TRUE), ])
}



#' Return transportability parameters
#' 
#' Convert transport = T/F into strengths for Y-U and M-U relationships
#' 
#' @param transport Whether transportability should hold
#' @return Length-4 vector of (yu_strength, mu_strength, small_yu_strength, 
#' small_mu_strength)
#' @export
convert_transport_to_strengths <- function(transport) {
  if (transport) {
    strengths <- rep(1.5, 4)
  } else {
    strengths <- c(1.5, 1.5, 0, 0)
  }
  return(strengths)
}



#' Mean-center a design matrix, except for intercept
#' 
#' Centers all non-intercept variables to zero, and sets interaction
#' columns to be centered at product of component variable means
#' 
#' @param xmat Design matrix to mean-center
#' @return Matrix of same size as xmat
#' @export
mean_xmat <- function(xmat) {
  
  # Mean-center everything
  xmat_mc <- matrix(colMeans(xmat), nrow = 1, 
                    dimnames = list(NULL, colnames(xmat)))
  
  # Check for interactions
  # If any exist, fix to product of means, not mean of products
  if (grep(":", colnames(xmat))) {
    intxs <- grep(":", colnames(xmat))
    for (i in intxs) {
      vars <- strsplit(colnames(xmat)[i], split = ":")[[1]]
      xmat_mc[, i] <- xmat_mc[, vars[1]] * xmat_mc[, vars[2]]
    }
  }
  
  return(xmat_mc)
}



#' Function to calculate gamma correction factor (p. 84)
#' 
#' Ratio gamma factor for interaction-respecting correction
#' Not the same as gamma factor in CDE bias correction
#' 
#' @param formula_y Formula for outcome regression, including U
#' @param df Data frame for fitting regression model
#' @param covariate_df Data frame containing covariate distribution (need not have U)
#' @param mean Whether mean correction factor should be returned, if convergence
#' problems require the logistic regression model be used
#' @return Named list with gamma correction factor gam and 
#' convergence flag cprobflag
calculate_intx_gamma_factor <- function(formula_y, df, covariate_df = NULL, 
                                        mean = TRUE) {
  
  # Initial fit attempt - often has convergence issues
  fit_log <- try(glm(formula = formula_y, family = binomial(link = "log"), 
                     data = df), 
                 silent = TRUE)
  
  # Branch based on whether the model converged
  if (class(fit_log) != "try-error"){
    
    # Calculate gamma correction factor
    gamma_factor <- exp(coef(fit_log)["u"])
    
    # Set convergence problem flag to 0 for no problem at all
    cprobflag <- 0
    
  } else {
    
    # Try supplying starting values from a logistic regression
    fit_logit <- glm(formula = formula_y, family = binomial(link = "logit"), 
                     data = df)
    mustarts <- fitted(fit_logit)
    fit_log <- try(glm(formula = formula_y, family = binomial(link = "log"), 
                       data = df,
                       mustart = mustarts), silent = TRUE)
    
    # Branch based on whether the model converged after starting values
    if (class(fit_log) != "try-error") {
      
      # Calculate gamma correction factor
      gam <- exp(coef(fit_log)["u"])
      
      # Set convergence problem flag to 1, for intermediate problems
      cprobflag <- 1
      
    } else {
      
      # Out of options - use logistic regression approximation
      cov_df_u0 <- cov_df_u1 <- covariate_df
      cov_df_u0$u <- 0
      cov_df_u1$u <- 1
      
      # Calculate gamma correction factor(s)
      # Note: only an approximate - for average values of a and m
      num <- predict(fit_logit, newdata = cov_df_u1, type = "response")
      den <- predict(fit_logit, newdata = cov_df_u0, type = "response")
      gam <- as.numeric(exp(log(num) - log(den)))
      if (mean) {
        gam <- mean(gam)
      }
      # Set convergence problem flag to 2, the worst problem
      cprobflag <- 2
    }
  }
  
  # Return gamma and a flag for whether encountered convergence issue
  return(list(gam = gam, cprobflag = cprobflag))
  
}



#' Calculate delta(s) for gamma-delta correction
#' 
#' Calculates the difference in U prevalence by A given M, C (p. 77)
#' In correction factor for CDE because it assume no exp-med intx in outcome
#'
#' @param df Data containing n observations
#' @param covariate_df Data frame containing baseline confounder distribution, 
#' if different from from df
#' @param mean Whether to return mean across covariate pattern or 
#' @return Scalar (if mean = TRUE) or length-n vector of delta correction factors
#' @export
calculate_dg_delta_factor <- function(df, covariate_df = NULL, mean = TRUE) {
  
  # Default to covariate pattern in df if no covariate_df is provided
  if (is.null(covariate_df)) {
    covariate_df <- df
  }
  
  # Fit model with U as a function of M, despite time ordering
  fit_u_with_m <- glm(u ~ 1 + z1 + z2 + m + a, data = df,
                      family = binomial(link = "logit"))
  # m = 0 is reference level
  covariate_df$m <- 0
  cov_df_a0 <- cov_df_a1 <- covariate_df
  cov_df_a0$a <- 0
  cov_df_a1$a <- 1
  
  delta <- predict(fit_u_with_m, newdata = cov_df_a1, type = "response") -
           predict(fit_u_with_m, newdata = cov_df_a0, type = "response")
  
  if (mean) {
    delta <- mean(delta)
  }
  
  return(delta)
}



#' Calculate gamma(s) for gamma-delta correction
#' 
#' Calculates the difference in Y prevalence by U given A, M, Z (p. 77)
#' In correction factor for CDE because it assume no exp-med intx in outcome
#'
#' @param df Data containing n observations
#' @param covariate_df Data frame containing baseline confounder distribution, 
#' if different from from df
#' @param mean Whether to return mean across covariate pattern or 
#' @return Scalar (if mean = TRUE) or length-n vector of delta correction factors
#' @export
calculate_dg_gamma_factor <- function(df, covariate_df = NULL, 
                                      mean = TRUE) {
  
  # Default to covariate pattern in df if no covariate_df is provided
  if (is.null(covariate_df)) {
    covariate_df <- df
  }
  
  # Fit model of Y as function of U, M, and confounders
  fit_y <- glm(y ~ 1 + z1 + z2 + a + m + u, data = df,
                      family = binomial(link = "logit"))
  # U effect is not going to constant across levels of A, so evaluate at a = 0.5
  # m = 0 is reference level
  covariate_df$m <- 0
  covariate_df_a0p5 <- covariate_df
  covariate_df_a0p5$a <- 0.5
  cov_df_u0 <- cov_df_u1 <- covariate_df_a0p5
  cov_df_u0$u <- 0
  cov_df_u1$u <- 1
  
  gam <- predict(fit_y, newdata = cov_df_u1, type = "response") -
         predict(fit_y, newdata = cov_df_u0, type = "response")
  
  if (mean) {
    gam <- mean(gam)
  }
  
  return(gam)
}



#' Function to get counts of a (Z1, Z2) covariate pattern from a data frame
#' 
#' @param df Data frame for counting (Z1, Z2) patterns
#' @return Vector of counts for (0,0), (1,0), (0,1), and (1,1)
#' @export
get_z1_z2_counts <- function(df) {
  df$count <- 1
  cov_df <- expand.grid(z1 = 0:1, z2 = 0:1)
  cov_df$id <- 1:nrow(cov_df)
  cov_with_counts <- merge(x = cov_df, 
                           y = aggregate(count ~ z1 + z2, FUN = sum, 
                                         data = df),
                           by = c("z1", "z2"), 
                           all.x = TRUE)
  cov_with_counts <- cov_with_counts[order(cov_with_counts$id), ]
  cov_with_counts$count[is.na(cov_with_counts$count)] <- 0
  return(cov_with_counts$count)
}



#' Calculate the naive/uncorrected (r)NDE estimate
#' 
#' @param df Data frame containing y, z1, z2, a, m
#' @param am_intx Whether to include exposure-mediator interaction in
#' outcome model
#' @param mean Whether to report mean (r)NDE or covariate-specific
#' @return Average (r)NDE, marginalized over (z1, z2) pattern that appears
#' in df
#' @export
run_uncorrected <- function(df, am_intx, mean = TRUE) {
  
  if (am_intx) {
    formula_y <- formula(y ~ 1 + z1 + z2 + a + m + a:m)
  } else {
    formula_y <- formula(y ~ 1 + z1 + z2 + a + m)
  }
  
  # zeroes added on are for u coefficient
  alpha_naive <- c(coef(glm(formula_y, data = df,
                            family = binomial(link = "logit"))), 0)
  
  beta_naive  <- c(coef(glm(m ~ 1 + z1 + z2 + a, data = df,
                            family = binomial(link = "logit"))), 0)
  
  naive_nder <- calculate_nder(params = NULL, 
                               alpha = alpha_naive, 
                               beta = beta_naive, 
                               gamma = rep(0, 3), # does not matter
                               u_ei = 0, am_intx = am_intx,
                               yu_strength = 0, mu_strength = 0)
  if (mean) {
    naive_nder <- weighted.mean(naive_nder, w = get_z1_z2_counts(df))
  } else {
    #stop("Not implemented!")
  }
  return(naive_nder)
}



#' Get a calculated (r)NDE, corrected using the delta-gamma
#' correction derived from a small data set
#' 
#' Has to assume no exposure-mediator interaction
#' 
#' @param df Data frame containing y, a, m, z1, z2
#' @param small_df Data frame with secondary data (includes u)
#' @return Corrected (r)NDE estimate using delta-gamma method
#' @export
run_delta_gamma <- function(df, small_df) {
  
  # Calculate gamma and delta correction factors for each possible pattern
  # only need m = 0 because that is reference level for intervention
  cov_df <- expand.grid(z1 = 0:1, z2 = 0:1, m = 0)
  cov_df$delta   <- calculate_dg_delta_factor(df = small_df, 
                                              covariate_df = cov_df, 
                                              mean = FALSE)
  cov_df$gamma   <- calculate_dg_gamma_factor(df = small_df, 
                                              covariate_df = cov_df, 
                                              mean = FALSE)
  
  # Get big data frequency weights for each covariate pattern
  df$count <- 1
  cov_df$id <- 1:nrow(cov_df)
  cov_with_counts <- merge(x = cov_df, 
                           y = aggregate(count ~ z1 + z2, FUN = sum, 
                                         data = df),
                           by = c("z1", "z2"), 
                           all.x = TRUE)
  cov_with_counts <- cov_with_counts[order(cov_with_counts$id), ]
  cov_with_counts$count[is.na(cov_with_counts$count)] <- 0
  
  # Naive models for delta-gamma exclude A-M interaction
  # Zeros correspond to U
  alpha_naive <- c(coef(glm(y ~ 1 + z1 + z2 + a + m, data = df,
                            family = binomial(link = "logit"))), 0)
  beta_naive  <- c(coef(glm(m ~ 1 + z1 + z2 + a, data = df,
                            family = binomial(link = "logit"))), 0)
  naive_nder <- calculate_nder(params = NULL, 
                               alpha = alpha_naive, 
                               beta = beta_naive, 
                               gamma = rep(0, 3), # does not matter
                               u_ei = 0, am_intx = 0,
                               yu_strength = 0, mu_strength = 0)
  
  # Calculate corrected NDER values for each covariate pattern
  dg_nder <- naive_nder - (cov_with_counts$delta * cov_with_counts$gamma)
  # print((cov_with_counts$delta * cov_with_counts$gamma))
  
  # Take weighted mean by frequency of covariate pattern in big data set
  avg_dg_nder <- weighted.mean(dg_nder, w = cov_with_counts$count)
  return(avg_dg_nder)
  
}



#' Calculate the interaction-respecting corrected (r)NDE
#' 
#' See p. 486 for bias factor
#' 
#' @param df Big data frame containing y, a, m, z1, z2
#' @param small_df Small data containing y, a, m, z1, z2, u
#' @param am_intx Whether to include A-M interaction in outcome model
#' @return Estimated bias-corrected population (r)NDE
#' @export
run_intx_corr <- function(df, small_df, am_intx) {
  
  if (am_intx) {
    formula_y <- formula(y ~ 1 + z1 + z2 + a + m + a:m + u)
  } else {
    formula_y <- formula(y ~ 1 + z1 + z2 + a + m + u)
  }
  
  # Big data sample frequencies of covariate patterns
  cov_df <- expand.grid(z1 = 0:1, z2 = 0:1)
  cov_df$counts <- get_z1_z2_counts(df)

  # Reference level for bias correction is u'
  u_prime <- aggregate(u ~ z1 + z2, data = small_df[small_df$a == 0,], 
                       FUN = mean)[, "u"]
  a <- 1
  a_star <- 0
  
  # Fit models to calculate B_nde
  fit_y <- glm(formula_y, data = small_df, family = binomial(link = "logit"))
  fit_u <- glm(u ~ 1 + z1 + z2 + a + m, data = small_df, 
               family = binomial(link = "logit")) # weird but correct
  fit_m <- glm(m ~ 1 + z1 + z2 + a, data = small_df, 
               family = binomial(link = "logit"))
  
  # Calculate bias for every covariate pattern, summing across m and u
  cov_df$B_nde <- 0
  for (m in 0:1) {
    for (u in 0:1) {
      diff_y  <- predict(fit_y, newdata = cbind(cov_df, a = a, m = m, u = u), 
                         type = "response") - 
                 predict(fit_y, newdata = cbind(cov_df, a = a, m = m, u = u_prime), 
                         type = "response")
      diff_u1 <- predict(fit_u, newdata = cbind(cov_df, a = a, m = m),
                         type = "response") -
                 predict(fit_u, newdata = cbind(cov_df, a = a_star, m = m),
                         type = "response")
      diff_u0 <- (1 - predict(fit_u, newdata = cbind(cov_df, a = a, m = m),
                              type = "response")) -
                 (1 - predict(fit_u, newdata = cbind(cov_df, a = a_star, m = m),
                              type = "response"))
      p_m     <- m * predict(fit_m, newdata = cbind(cov_df, a = a_star), 
                             type = "response") +
                 (1 - m) * (1 - predict(fit_m, newdata = cbind(cov_df, a = a_star), 
                                        type = "response"))
      diff_u  <- (diff_u1 * u) + (diff_u0 * (1 - u))
      cov_df$B_nde <- cov_df$B_nde + diff_y * diff_u * p_m
    }
  }
  
  # Naive models can include A-M interaction
  # Zeros correspond to U
  alpha_naive <- c(coef(glm(update.formula(formula_y, . ~ . - u), data = df,
                            family = binomial(link = "logit"))), 0)
  beta_naive  <- c(coef(glm(m ~ 1 + z1 + z2 + a, data = df,
                            family = binomial(link = "logit"))), 0)
  gamma_naive <- coef(glm(u ~ 1 + z1 + z2, data = df,
                            family = binomial(link = "logit")))

  naive_nder <- calculate_nder(params = NULL, 
                               alpha = alpha_naive, 
                               beta = beta_naive, 
                               gamma = gamma_naive,
                               u_ei = 0, am_intx = am_intx,
                               yu_strength = 0, mu_strength = 0)
  
  # Subtract off bias and average by covariate pattern frequency
  intx_nder <- naive_nder - cov_df$B_nde
  avg_intx_nder <- weighted.mean(intx_nder, w = cov_df$counts)
  return(avg_intx_nder)
}



#' Run one full replicate with all frequentist comparator methods
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
#' @param result_type Whether to return full object ("raw") or only selected
#' information about NDER
#' @param B Number of bootstrap samples to take for CIs
#' @return Named list of results
#' @export
run_frequentist_replicate <- function(n, u_ei, am_intx, 
                              yu_strength, mu_strength, params = NULL,
                              small_yu_strength, small_mu_strength, 
                              small_params = NULL,
                              result_type = c("raw", "processed"),
                              B = 500,
                              ...) {
  
  # Basic checks
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei, am_intx, yu_strength, mu_strength)
  }

  if (is.null(small_params)) {
    small_params <- return_dgp_parameters(u_ei, am_intx, 
                                          small_yu_strength, small_mu_strength)
  }
  
  # Simulate data
  dl <- simulate_data(n = n, params = params)
  small_dl <- simulate_data(n = n, params = small_params)
  
  # Get true value
  truth_nder <- calculate_nder(params = NULL,
                               u_ei = params$u_ei,
                               am_intx = params$am_intx,
                               yu_strength = yu_strength,
                               mu_strength = mu_strength,
                               mean = TRUE)
  
  # Run different corrections
  nder_uc <- run_uncorrected(dl$df, am_intx = 0, mean = TRUE)
  nder_dg <- run_delta_gamma(dl$df, small_dl$df)
  nder_ix <- run_intx_corr(dl$df, small_dl$df, am_intx = am_intx)

  ci_uc <- quantile(replicate(B, run_uncorrected(boot_samp(dl$df), 
                                                 am_intx = am_intx, 
                                                 mean = TRUE)),
                    probs = c(0.025, 0.975))
  ci_dg <- quantile(replicate(B, run_delta_gamma(boot_samp(dl$df),
                                                 small_dl$df)),
                    probs = c(0.025, 0.975))
  ci_ix <- quantile(replicate(B, run_intx_corr(boot_samp(dl$df), 
                                               small_dl$df, 
                                               am_intx = am_intx)),
                    probs = c(0.025, 0.975))
  cov_uc  <- check_coverage(ci = ci_uc, truth = truth_nder)
  cov_dg  <- check_coverage(ci = ci_dg, truth = truth_nder)
  cov_ix  <- check_coverage(ci = ci_ix, truth = truth_nder)
  width_uc <- get_ci_width(ci_uc)
  width_dg <- get_ci_width(ci_dg)
  width_ix <- get_ci_width(ci_ix)
  
  if (result_type == "raw") {
    return(list(params = params, data_list = dl,
                truth_nder = truth_nder,
                estimate = c(uc = nder_uc, dg = nder_dg, ix = nder_ix),
                bias = c(uc = nder_uc - truth_nder,
                         dg = nder_dg - truth_nder,
                         ix = nder_ix - truth_nder),
                ci = c(uc = ci_uc, dg = ci_dg, ix = ci_ix),
                ci_cov = c(uc = cov_uc, dg = cov_dg, ix = cov_ix),
                ci_width = c(uc = width_uc, dg = width_dg, ix = width_ix)))
  } else {
    return(list(truth_nder = truth_nder,
                estimate = c(uc = nder_uc, dg = nder_dg, ix = nder_ix),
                bias = c(uc = nder_uc - truth_nder,
                         dg = nder_dg - truth_nder,
                         ix = nder_ix - truth_nder),
                ci = c(uc = ci_uc, dg = ci_dg, ix = ci_ix),
                ci_cov = c(uc = cov_uc, dg = cov_dg, ix = cov_ix),
                ci_width = c(uc = width_uc, dg = width_dg, ix = width_ix)))  
  }
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
#' @param result_type Whether to return full object ("raw") or only selected
#' information about NDER
#' @param ... Additional parameters to be passed to bin_bin_sens_stan
#' @return Stan model fit object
#' @export
run_bdf_replicate <- function(n, u_ei, am_intx, 
                              yu_strength, mu_strength, params = NULL,
                              small_yu_strength, small_mu_strength, 
                              small_params = NULL, dd_control = NULL,
                              prior_type = "dd",
                              result_type = c("raw", "processed"),
                              ...) {
  
  # Basic checks
  stopifnot(prior_type %in% c("unit", "partial", "strict", "dd"))
  if (is.null(params)) {
    params <- return_dgp_parameters(u_ei, am_intx, yu_strength, mu_strength)
  }
  params$prior <- prior_type
  
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
