#' Function to inflate variance for identifiable parameters
#' 
#' @param mat Variance-covariance matrix to inflate rows/columns
#' @param except_index Vector of row/column indices NOT to inflate
#' @param inf_fact Inflation factor (on standard deviation scale)
#' (identifiable parameters' variances are scaled by square of inf_fact)
#' @return Variance-covariance matrix with inflated variances for some parameters
make_da_noninformative <- function(mat, except_index, inf_fact){
  mat <- mat * inf_fact^2
  mat[except_index, ] <- mat[except_index, ] / inf_fact
  mat[, except_index] <- mat[, except_index] / inf_fact
  return(mat)
}



#' Function to make priors
#' 
#' @param df Data frame for fitting regressions to make priors
#' @param formulas List of regression formulas
#' @param inf_fact Inflation factor for non-informative parameters (default is 1000)
#' @return List of priors, made non-informative except for U
make_da_priors <- function(df, formulas, inf_fact = 1000) {
  
  # Extract priors from mediator regression
  if (length(unique(df$M)) > 2) {
    fitM <- VGAM::vglm(formulas$mediatorU, data = df, 
                       family = VGAM::multinomial(refLevel = 1))
  } else {
    fitM <- glm(formulas$mediatorU, family = binomial(link = "logit"), 
                data = df)
  }
  m_coef_M <- as.vector(coef(fitM))
  V_coef_M_raw <- unname(vcov(fitM))
  
  which_u_M <- grep("U:", names(coef(fitM)))
  V_coef_M <- make_da_noninformative(V_coef_M_raw, which_u_M, inf_fact)
  
  
  # Extract priors from outcome regression
  fitY <- glm(formulas$outcomeU, family = binomial(link = "logit"), data = df)
  m_coef_Y <- as.vector(coef(fitY))
  V_coef_Y_raw <- unname(vcov(fitY))
  
  which_u_Y <- grep("U", names(coef(fitY)))
  V_coef_Y <- make_da_noninformative(V_coef_Y_raw, which_u_Y, inf_fact)
  
  # Extract priors from unmeasured confounder regression
  fitU <- glm(formulas$unmeasured, family = binomial(link = "logit"), data = df)
  m_coef_U <- as.vector(coef(fitU))
  V_coef_U <- unname(vcov(fitU))
  
  # Make result
  priors <- list()
  priors$outcome[["mean"]] <- m_coef_Y
  priors$outcome[["vcov"]] <- V_coef_Y
  priors$mediator[["mean"]] <- m_coef_M
  priors$mediator[["vcov"]] <- V_coef_M
  priors$unmeasured[["mean"]] <- m_coef_U
  priors$unmeasured[["vcov"]] <- V_coef_U
  
  # Return
  return(priors)
  
}



#' Function to run the Stan medBSA function
#' 
#' @param big_df Main data set for fitting the models
#' @param small_df Small data set for constructing the priors
#' @param mediator Name of mediator variable in both data sets
#' @param outcome Name of outcome variable in both data sets
#' @param tx Exposure variable
#' @param unmeasured Name of mediator variable in small data set
#' @param formulas List of formulas
#' @param am_intx Exposure-mediator interaction indicator
#' @param inf_fact Inflation factor for making noninformative priors
#' @return Stan model fit
run_stan_BSA <- function(big_df, small_df, 
                         mediator, outcome, tx, unmeasured,
                         formulas, am_intx, inf_fact = 100, ...){
  
  # Alias the unmeasured confounder, mediator, and outcome
  small_df[["Y"]] <- small_df[[outcome]]
  small_df[["M"]] <- small_df[[mediator]]
  small_df[["U"]] <- small_df[[unmeasured]]
  big_df[["Y"]]   <- big_df[[outcome]]
  big_df[["M"]]   <- big_df[[mediator]]
  
  # Update formulas with shortcut outcome name and add U to other formulas
  formulas$u_noU      <- formulas$unmeasured
  formulas$unmeasured <- update.formula(formulas$unmeasured, U ~ .)
  formulas$mediator   <- update.formula(formulas$mediator, M ~ .)
  formulas$outcome    <- update.formula(formulas$outcome, Y ~ .)
  formulas$mediatorU  <- update.formula(formulas$mediator, . ~ . + U)
  formulas$outcomeU   <- update.formula(formulas$outcome, . ~ . + U)
  
  # Isolate baseline covariates from mediator regression model
  covs <- attr(terms(formulas$mediator), "term.labels")
  formulas$baseline <- as.formula(paste("~", paste(covs[covs != tx], 
                                                   collapse  = " + ")))
  
  # Build design matrices
  # Build aggregation and weighting objects
  N_tot <- NROW(big_df)
  cat(N_tot)
  big_df$w <- 1
  big_df <- aggregate(w ~ ., data = big_df, FUN = sum)
  big_df[["U"]] <- NA
  Xmat_M  <- model.matrix(formulas$mediator, big_df)
  Xmat_Y  <- model.matrix(formulas$outcome, big_df)
  Xmat_U  <- model.matrix(formulas$u_noU, big_df)
  Xmat_bl <- model.matrix(formulas$baseline, big_df)
  N       <- nrow(Xmat_Y)
  cat(" ")
  cat(N)
  
  # Construct priors
  priors <- make_da_priors(small_df, formulas, inf_fact = inf_fact)
  
  # Make stan data list
  stan_data <- list(K_m = length(unique(big_df$M)),
                    D_m = ncol(Xmat_M) + 1,
                    D_y = ncol(Xmat_Y) + 1,
                    D_u = ncol(Xmat_U),
                    D_bl = ncol(Xmat_bl),
                    N = N,
                    am_intx = am_intx,
                    x_mu1 = cbind(Xmat_M, U = rep(1, N)),
                    x_mu0 = cbind(Xmat_M, U = rep(0, N)),
                    x_yu1 = cbind(Xmat_Y, U = rep(1, N)),
                    x_yu0 = cbind(Xmat_Y, U = rep(0, N)),
                    x_u = Xmat_U,
                    x_bl = Xmat_bl,
                    m_coef_u = priors$unmeasured[["mean"]],
                    v_coef_u = priors$unmeasured[["vcov"]],
                    m_coef_m = priors$mediator[["mean"]],
                    v_coef_m = priors$mediator[["vcov"]],
                    m_coef_y = priors$outcome[["mean"]],
                    v_coef_y = priors$outcome[["vcov"]],
                    m = big_df[[mediator]],
                    y = big_df[[outcome]],
                    N_tot = N_tot,
                    w = big_df[["w"]])
  
  # Run
  samples <- rstan::sampling(stanmodels$mediation_model_ncp, 
                             data = stan_data, ...)
  
  # Return
  return(samples)
}



#' Run the data application
#' 
#' @param seer_file_path Path to SEER .txt file
#' @param cancors_file_path Path to CanCORS .csv file
#' @param samples_file_path Where to save Stan samples
#' @param am_intx Whether to have exposure-mediator interaction
#' @param inf_fact Inflation factor for priors
#' @param chains Number of chains to run
#' @param iter Number of MCMC iterations per chain
#' @param seed Random number seed for Stan
#' @param auto_write Whether to write the model to disk
#' @param mc.cores Number of cores for running Stan model
#' @return Stan fit (+ side effect of saving stan fit as RDS)
#' @export
#' @examples 
#' \dontrun{
#' run_data_app(seer_file_path = "Data/Raw/data_SEER_for_BSA_summer_project.txt", 
#'              cancors_file_path = "Data/Raw/selected_cancors_data_2016_3_14.csv",
#'              save_file_path = paste0("data_app_res_",
#'                                      format(Sys.Date(), "%Y%m%d"), ".rds"),
#'              am_intx = 1, inf_fact = 100, chains = 4, iter = 2000, seed = 42,
#'              auto_write = TRUE, mc.cores = 4)
#' }
run_data_app <- function(seer_file_path, 
                         cancors_file_path,
                         samples_file_path = paste0("data_app_samples_",
                                                    format(Sys.Date(), "%Y%m%d"),
                                                    ".csv"),
                         am_intx = 1,
                         inf_fact = 100,
                         chains = 4, iter = 2000, seed = 42,
                         auto_write = TRUE, mc.cores = 4, ...) {
  
  rstan::rstan_options(auto_write = auto_write)
  options(mc.cores = mc.cores)
  
  # Read in raw data
  cancors <- read.csv(cancors_file_path)
  seer    <- read.table(seer_file_path)
  
  # Process
  # Make middle age category the reference
  seer$age_cat <- relevel(as.factor(seer$age_cat), ref = 2)
  seer$female  <- ifelse(seer$female == 2, 1, 0)
  
  # Collapse regions to match CanCors
  seer$region[seer$region == 2] <- 1
  seer$region  <- factor(seer$region,
                         labels = c("Other", "South", "West"))
  
  # Make Cancors variable names match SEER
  cancors              <- na.omit(cancors)
  cancors$stage_cat    <- cancors$stage
  cancors$age_cat      <- relevel(as.factor(cancors$age_dx), ref = 2)
  cancors$female       <- ifelse(cancors$gender == 2, 1, 0)
  cancors$pov          <- ifelse(cancors$income == 1, 1, 0)
  cancors$black        <- cancors$race
  cancors$region       <- factor(cancors$region,
                                 labels = c("Other", "South", "West"))
  cancors$fiveyearsurv <- ifelse(cancors$surv < 365.25 * 5, 0, 1)
  
  # Make formulas to pass in
  # Leave out unmeasured confounder in mediator and outcome formulas!
  formulas <- list()
  formulas$mediator   <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  formulas$unmeasured <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  if (am_intx) {
    formulas$outcome <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                + black + as.factor(stage_cat) 
                                + black*as.factor(stage_cat))
  } else {
    formulas$outcome <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                + black + as.factor(stage_cat))
  }
  
  # Run
  bayes_res <- run_stan_BSA(big_df = seer, small_df = cancors,
                            mediator = "stage_cat",
                            outcome = "fiveyearsurv", tx = "black",
                            unmeasured = "pov",
                            formulas,
                            am_intx = am_intx,
                            inf_fact = inf_fact,
                            chains = chains, iter = iter, seed = seed,
                            sample_file = samples_file_path,
                            cores = mc.cores,
                            ...)
  
  return(bayes_res)
}



#' Calculate the total effect for the data application
#' 
#' The total effect here is the existing racial disparity
#' 
#' @param seer_file_path Path to SEER data
#' @return Scalar frequentist estimate of the total black-white disparity in 
#' 5-year survival
#' @export
calculate_total_effect <- function(seer_file_path) {
  
  # Read in SEER
  seer <- read.table(seer_file_path)
  
  # Make middle age category the reference
  seer$age_cat <- relevel(as.factor(seer$age_cat), ref = 2)
  seer$female  <- ifelse(seer$female == 2, 1, 0)
  
  # Collapse regions to match CanCors
  seer$region[seer$region == 2] <- 1
  seer$region  <- factor(seer$region,
                         labels = c("Other", "South", "West"))
  
  # Make formula for outcome regression and fit logistic regression model
  formula_outcome <- formula(fiveyearsurv ~ female + as.factor(age_cat) + 
                               as.factor(region) + black)
  fit_y <- glm(formula_outcome, data = seer, 
               family = binomial(link = "logit"))
  
  # Calculate individual-level predictions for marginalization over confounders
  seer_a1 <- seer_a0 <- seer
  seer_a0$black <- 0
  seer_a1$black <- 1
  te_indiv <- predict(fit_y, newdata = seer_a1, type = "response") - 
              predict(fit_y, newdata = seer_a0, type = "response")
  te <- mean(te_indiv)
  
  return(te)
}

