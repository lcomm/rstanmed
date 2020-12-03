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



#' Create a block matrix with diagonal repeats
#'
#' @param mat An RxC Matrix to repeat
#' @param times Number of times to repeat
#' @return A matrix with row and column counts Rxtimes rows and Cxtimes columns
make_repeated_block_matrix <- function(mat, times) {
  stopifnot(is.matrix(mat))
  stopifnot(is.numeric(times))
  return(kronecker(diag(times), mat))
}



#' Resize R matrix based on priors
#'
#' Make the R_mat matrix bigger with a block-identity, and/or through block
#' repeats. Adding a block-identity will happen if it is detected that the 
#' R matrix was missing a confounder (i.e., there is an unmeasured confounder),
#' and the block repeats will happen if the length of elements in the informative
#' prior suggests that the prior is actually on a coefficient matrix (i.e., the
#' model is for a multinomial logistic model outcome).
#' 
#' @param R_mat R matrix from QR decomposition of design matrix for that model,
#' with any rescaling already performed by \code{\link{resize_R}}
#' @param model_label Name of submodel to scale. Must be in \code{unmeasured}, 
#' \code{mediator}, or \code{outcome}.
#' @param prior_list Named list, as made by \code{\link{make_da_priors}}.
#' @return R matrix appropriate for pre- and most-multiplications
resize_Rmat <- function(R_mat, model_label, prior_list) {
  stopifnot(model_label %in% c("mediator", "unmeasured", "outcome"))
  stopifnot("mean" %in% names(prior_list[[model_label]]))
  stopifnot("vcov" %in% names(prior_list[[model_label]]))
  stopifnot(is.matrix(R_mat))
  stopifnot(dim(R_mat)[1] == dim(R_mat)[2])
  dimR <- dim(R_mat)[1]
  dimMean <- length(prior_list[[model_label]][["mean"]])
  # At present, code not set up for multiple missing confounders
  size_test <-c(dimMean %% dimR == 0, dimMean %% (dimR + 1) == 0)
  stopifnot(any(size_test))
  diffSize <- which(size_test) - 1
  stopifnot(diffSize %in% c(0, 1)) # extra sanity check
   
  # Augment with ones if prior includes an extra variable
  if (diffSize > 0) {
    newR <- diag(dimR + diffSize)
    newR[1:dimR, 1:dimR] <- R_mat
    R_mat <- newR
  }
  
  # Expand block-wise if necessary (i.e., for multicategory logit model)
  if (dimMean > dimR + 1) {
    n_block <- dimMean / (dimR + diffSize)
    stopifnot(dimMean %% (dimR + diffSize) == 0)
    R_mat <- make_repeated_block_matrix(mat = R_mat, times = n_block)
  }
  
  return(R_mat)
}



#' Basic submodel-level QR decomposition and scaling
#' 
#' @param Xmat Design matrix
#' @param scaling Scaling factor for Q (R multiplied by 1 / scaling)
#' @return Named list of Qstar and Rstar matrices
QR_to_scaled_list <- function(Xmat, scaling) {
  res <- list(qr_res = qr(Xmat))
  res$Qstar <- qr.Q(res$qr_res) * scaling
  res$Rstar <- qr.R(res$qr_res) / scaling
  return(res)
}



#' Apply QR decomposition to priors
#' 
#' Given an informative prior, adjust it based on the R (or, more like, R^*) 
#' matrix from a QR decomposition. This allows the Stan sampling to be performed
#' on the QR-scaled parameters because the informative prior is converted to the
#' QR-scaled version. 
#' 
#' @param prior_list Named list, as made by \code{\link{make_da_priors}}.
#' @param model_label Name of submodel to scale. Must be in \code{unmeasured}, 
#' \code{mediator}, or \code{outcome}.
#' @param R_mat R matrix from QR decomposition of design matrix for that model,
#' with any rescaling already performed by \code{\link{resize_R}}
#' @return \code{prior_list} list object with the \code{model_label} 
#' informative prior translated to the QR scale
QR_scale_submodel_prior <- function(prior_list, model_label, R_mat) {
  
  # Apply pre- and post-multiplications
  prior_list[[model_label]][["mean"]] <- 
    c(R_mat %*% prior_list[[model_label]][["mean"]])
  prior_list[[model_label]][["vcov"]] <- 
    c(R_mat %*% prior_list[[model_label]][["vcov"]] %*% t(R_mat))
  
  return(prior_list)
}



#' Do everything needed to apply a QR decomposition to the fit process
#' 
#' @param Xmat_M M matrix
#' @param Xmat_U U matrix
#' @param Xmat_Y Y matrix
#' @param scaling Scaling determined by \code{scale_factor} input to 
#' \code{\link{run_stan_BSA}}
#' @param prior_list List of priors as made by \code{\link{make_da_priors}}
#' @return Named list with Q, Rstar, Rstar_inv matrices and list of adjusted
#' priors
#' @export
calculate_QR <- function(Xmat_M, Xmat_U, Xmat_Y, scaling, prior_list) {
  
  # Basic QR decomposition output list
  res <- list(M = QR_to_scaled_list(Xmat = Xmat_M, scaling = scaling),
              U = QR_to_scaled_list(Xmat = Xmat_U, scaling = scaling),
              Y = QR_to_scaled_list(Xmat = Xmat_Y, scaling = scaling))
  
  # Resizing Rstar as needed
  res[["M"]][["Rstar"]] <- resize_Rmat(res[["M"]][["Rstar"]], 
                                       model_label = "mediator", 
                                       prior_list = prior_list)
  res[["U"]][["Rstar"]] <- resize_Rmat(res[["U"]][["Rstar"]], 
                                       model_label = "unmeasured", 
                                       prior_list = prior_list)
  res[["Y"]][["Rstar"]] <- resize_Rmat(res[["Y"]][["Rstar"]], 
                                       model_label = "outcome", 
                                       prior_list = prior_list)
  
  # Inverting the resized Rstars
  res[["M"]][["Rstar_inv"]] <- solve(res[["M"]][["Rstar"]])
  res[["U"]][["Rstar_inv"]] <- solve(res[["U"]][["Rstar"]])
  res[["Y"]][["Rstar_inv"]] <- solve(res[["Y"]][["Rstar"]])
  
  # Applying QR scaling to priors
  res$priors <- prior_list
  res$priors[["M"]] <- QR_scale_submodel_prior(prior_list = res$priors, 
                                               model_label = "mediator", 
                                               R_mat = res[["M"]][["Rstar"]])
  res$priors[["U"]] <- QR_scale_submodel_prior(prior_list = res$priors, 
                                               model_label = "unmeasured", 
                                               R_mat = res[["U"]][["Rstar"]])
  res$priors[["Y"]] <- QR_scale_submodel_prior(prior_list = res$priors, 
                                               model_label = "outcome", 
                                               R_mat = res[["Y"]][["Rstar"]])
  return(res)
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
#' @param mean_only Whether the final step should just take the mean of Y
#' @param scale_factor Scaling factor for QR decomposition. (Makes no difference
#' with respect to the assumed model, just for computational advantage.)
#' @return Stan model fit
run_stan_BSA <- function(big_df, small_df, 
                         mediator, outcome, tx, unmeasured,
                         formulas, am_intx, inf_fact = 1,
                         mean_only = 1, 
                         scale_factor = "sqrt", 
                         ...){
  
  scale_factor <- match.arg(scale_factor, choices = c("sqrt", "N", "1"))
  
  # Alias the unmeasured confounder, mediator, and outcome
  small_df[["Y"]] <- small_df[[outcome]]
  small_df[["M"]] <- small_df[[mediator]]
  small_df[["U"]] <- small_df[[unmeasured]]
  big_df[["Y"]]   <- big_df[[outcome]]
  big_df[["M"]]   <- big_df[[mediator]]
  K_m <- length(unique(big_df$M))
  
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
  big_df$w <- 1
  big_df <- aggregate(w ~ ., data = big_df, FUN = sum)
  big_df[["U"]] <- NA
  Xmat_M  <- model.matrix(formulas$mediator, big_df)
  Xmat_Y  <- model.matrix(formulas$outcome, big_df)
  Xmat_U  <- model.matrix(formulas$u_noU, big_df)
  N       <- nrow(Xmat_Y)
  
  # Construct priors
  priors <- make_da_priors(small_df, formulas, inf_fact = inf_fact)
  
  # Perform QR decomposition and modify priors accordingly priors
  scaling <- switch(scale_factor,
                    "sqrt" = sqrt(N_tot - 1),
                    "N" = N_tot,
                    "1" = 1)
  
  qr_res <- calculate_QR(Xmat_M = Xmat_M, Xmat_U = Xmat_U, Xmat_Y = Xmat_Y, 
                         scaling = scaling, prior_list = priors)
  
  # Make stan data list
  stan_data <- list(K_m = length(unique(big_df$M)),
                    D_m = ncol(Xmat_M) + 1,
                    D_y = ncol(Xmat_Y) + 1,
                    D_u = ncol(Xmat_U),
                    N = N,
                    am_intx = am_intx,
                    x_mu1 = cbind(qr_res[["M"]][["Qstar"]], U = rep(1, N)),
                    x_mu0 = cbind(qr_res[["M"]][["Qstar"]], U = rep(0, N)),
                    x_yu1 = cbind(qr_res[["Y"]][["Qstar"]], U = rep(1, N)),
                    x_yu0 = cbind(qr_res[["Y"]][["Qstar"]], U = rep(0, N)),
                    x_u = qr_res[["U"]][["Qstar"]],
                    m_coef_u = qr_res$priors$unmeasured[["mean"]],
                    v_coef_u = qr_res$priors$unmeasured[["vcov"]],
                    m_coef_m = qr_res$priors$mediator[["mean"]],
                    v_coef_m = qr_res$priors$mediator[["vcov"]],
                    m_coef_y = qr_res$priors$outcome[["mean"]],
                    v_coef_y = qr_res$priors$outcome[["vcov"]],
                    m = big_df[[mediator]],
                    y = big_df[[outcome]],
                    w = big_df[["w"]],
                    Rstar_inv_M = qr_res[["M"]][["Rstar_inv"]],
                    Rstar_inv_Y = qr_res[["Y"]][["Rstar_inv"]],
                    Rstar_inv_U = qr_res[["U"]][["Rstar_inv"]],
                    mean_only = mean_only)
  
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
#'              samples_file_path = paste0("data_app_res_",
#'                                  format(Sys.Date(), "%Y%m%d"), ".rds"),
#'              am_intx = 1, inf_fact = 100, chains = 4, iter = 2000, seed = 42,
#'              auto_write = TRUE, mc.cores = 4)
#' }
run_data_app <- function(seer_file_path, 
                         cancors_file_path,
                         samples_file_path = paste0("data_app_samples_",
                                                    format(Sys.Date(), "%Y%m%d"),
                                                    ".csv"),
                         am_intx = 1,
                         inf_fact = 1,
                         scale_factor = "sqrt",
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
                            scale_factor = scale_factor,
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

