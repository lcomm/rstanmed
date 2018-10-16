#' Take the softmax of a matrix, row-wise
#' 
#' @param x N x K Matrix
#' @return N x K matrix with row sums equal to 1
row_softmax <- function(x){
  score_exp <- exp(x)
  probs <- sweep(score_exp, 1, rowSums(score_exp), "/")
  return(probs)
}



#' Calculate the naive residual disparity for data application
#' 
#' @param seer_file_path Path to raw SEER txt file
#' @param am_intx Whether to include A * M interaction in Y model
#' @param B Number of bootstrap samples to take
#' @return Vector of B residual disparities
#' @export
calculate_data_app_naive <- function(seer_file_path, am_intx = 1, B = 2000) {
  
  # Read in seer from file path, and set 
  seer       <- read.table(path.expand(seer_file_path))
  mediator   <- "stage_cat"
  outcome    <- "fiveyearsurv" 
  tx         <- "black"
  unmeasured <- "pov"
  n          <- NROW(seer)
  K_m        <- length(unique(seer[[mediator]]))
  
  # Make middle age category the reference
  seer$age_cat <- relevel(as.factor(seer$age_cat), ref = 2)
  seer$female  <- ifelse(seer$female == 2, 1, 0)
  
  # Collapse regions to match CanCors
  seer$region[seer$region == 2] <- 1
  seer$region  <- factor(seer$region,
                         labels = c("Other", "South", "West"))
  
  formulas <- list()
  formulas$mediator   <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  formulas$unmeasured <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  if (am_intx) {
    formulas$outcome  <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black + as.factor(stage_cat)
                                 + black*as.factor(stage_cat))
  } else {
    formulas$outcome  <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black + as.factor(stage_cat))
  }
  
  # These actually get used now, except for U
  seer[["Y"]] <- seer[[outcome]]
  seer[["M"]] <- seer[[mediator]]
  seer[["U"]] <- sample(0:1, size = n, replace = TRUE) # gets overwritten
  seer[["A"]] <- seer[[tx]] #<- sample(0:1, size = n, replace = TRUE)
  
  # Construct formulas
  formulas$mediator   <- update.formula(formulas$mediator, M ~ .)
  formulas$outcome    <- update.formula(formulas$outcome, Y ~ .)
  
  nder_naive <- rep(NA, B)
  for (b in 1:B) {
    new_booted <- booted <- boot_samp(seer)
    fit_m <- VGAM::vglm(formulas$mediator, family = multinomial(refLevel = 1), 
                        data = booted)#, control = vglm.control(maxit = 500,
                                       #                       stepsize = 0.7))
    fit_y <- glm(formulas$outcome, family = binomial(link = "logit"), 
                 data = booted)
    
    new_booted$A <- new_booted$black <- 0
    pm_a0 <- VGAM::predict(fit_m, newdata = new_booted, type = "response")
    
    nder_indiv <- rep(0, NROW(booted))
    for (m in 1:4) {
      new_booted$M <- new_booted$stage_cat <- m
      new_booted$A <- new_booted$black <- 1
      py10 <- stats::predict(fit_y, newdata = new_booted, type = "response")
      new_booted$A <- new_booted$black <- 0
      py00 <- stats::predict(fit_y, newdata = new_booted, type = "response")
      nder_indiv <- nder_indiv + unname((py10 - py00) * pm_a0[, m])
    }
    nder_naive[b] <- mean(nder_indiv)
  }
  
  return(nder_naive)
}



#' Make arrays
#' 
#' @param Xmat_Y Design matrix for outcome regression
#' @param Xmat_M Design matrix for mediator regression
#' @param Xmat_U Design matrix for unmeasured confounder regression
#' @param am_intx 0/1 indicator of exposure-mediator interaction
#' @param n_patt Number of covariate patterns
#' @param K_m Number of levels for mediator
#' @return Named list containing arrays
#' @export
make_xy_xm_xu <- function(Xmat_Y, Xmat_M, Xmat_U, am_intx, n_patt, K_m) {
  
  # U
  D_u <- ncol(Xmat_U)
  xu <- array(Xmat_U, dim = c(n_patt, D_u, 2))
  wa_xu <- D_u
  
  # M
  D_m <- ncol(Xmat_M)
  xm <- array(Xmat_M, dim = c(n_patt, D_m, 2, 2))
  wa_xm <- D_u
  wu_xm <- D_m
  
  # Y
  D_y <- ncol(Xmat_Y)
  xy <- array(Xmat_Y, dim = c(n_patt, D_y, 2, 2, K_m))
  wa_xy <- D_u
  wm_xy <- (wa_xy + 1):(wa_xy + K_m - 1)
  if (am_intx) {
    wam_xy <- (max(wm_xy) + 1):(max(wm_xy) + K_m - 1)
  }
  wu_xy <- D_y
  
  # Fill in values for the design matrices
  for (a in 0:1) {
    xu[ , wa_xu, a + 1]    <- a
    xm[ , wa_xm, a + 1, ]  <- a
    xy[ , wa_xy, a + 1, , ] <- a
    
    for (u in 0:1) {
      xm[ , wu_xm,      , u + 1]   <- u
      xy[ , wu_xy, a + 1, u + 1, ] <- u
      
      for (m in 1:K_m) {
        # initialize all M dummies to 0
        xy[ , wm_xy, a + 1, u + 1, m] <- 0
        if (am_intx) {
          xy[ , wam_xy, a + 1, u + 1, m] <- 0
        }
        
        if (m != 1) {
          # make the one non-zero dummy
          xy[ , wm_xy[1] + (m - 2), a + 1, u + 1, m] <- 1
          if (am_intx && (a == 1)) {
            xy[ , wam_xy[1] + (m - 2), a + 1, u + 1, m] <- 1
          }
        }
        
      } # end m loop
    } # end u loop
  } # end a loop
  
  return(list(xy = xy, xm = xm, xu = xu))
}



#' Get distribution of M marginalized over U
#'
#' p(M = m | Z = z, A = a) = sum_u p(M = m | z, A = a, U = U(a)) P(U = u(a) | z)
#' 
#' @param pm Array with element (i,j,k,l) = P(M = j | Z_i, A = k - 1, U = l - 1)
#' @param pu1 N x 2 matrix with element (i,j) = P(U = 1 | Z_i, A = (j - 1))
#' @param a Value of A to evaluate 
#' @return N x 4 matrix with element (i,j) = P(M = j | Z = Z_i, A = a)
get_marginal_pm <- function(pm, pu1, a) {
  pmmarg <- pm[ , , a + 1, 1] * (1 - pu1[ , a + 1]) + 
            pm[ , , a + 1, 2] * (pu1[ , a + 1])  
  return(pmmarg)
}



#' Get distribution of Y marginalized over M and U
#' 
#' @param py1 Array with element (i,j,k,l) = 
#' P(Y = 1 | Z_i, A = (j - 1), U = (k - 1), M = l)
#' @param pmmarg N x 4 matrix with element (i,j) = P(M = j | Z = Z_i, A = a)
#' @param pu1 N x 2 matrix with element (i,j) = P(U = 1 | Z_i, A = (j - 1))
#' @param a Value of A to evaluate 
#' @return Length-N vector with element (i) = P(Y = 1 | A = a), where marginalization
#' over M and U has occurred with appropriate probabilities
get_marginal_ey <- function(py1, pmmarg, pu1, a) {
  py1margu <- py1[ , a + 1, 1, ] * (1 - pu1[ , a + 1]) + 
              py1[ , a + 1, 2, ] * (pu1[ , a + 1])
  return(rowSums(py1margu * pmmarg))
}



#' Make quartet of Y counterfactuals for (r)NDE and (r)NIE
#' 
#' @param py1 Array with element (i,j,k,l) = 
#' P(Y = 1 | Z_i, A = (j - 1), U = (k - 1), M = l)
#' @param pm Array with element (i,j,k,l) = P(M = j | Z_i, A = k - 1, U = l - 1)
#' @param pu1 N x 2 matrix with element (i,j) = P(U = 1 | Z_i, A = (j - 1))
#' @return N x 4 matrix
#' @export
make_quartet <- function(py1, pm, pu1) {
  
  pm_a1 <- get_marginal_pm(pm, pu1, a = 1)
  pm_a0 <- get_marginal_pm(pm, pu1, a = 0)
  
  ey00 <- get_marginal_ey(py1, pmmarg = pm_a0, pu1, a = 0)
  ey01 <- get_marginal_ey(py1, pmmarg = pm_a1, pu1, a = 0)
  ey10 <- get_marginal_ey(py1, pmmarg = pm_a0, pu1, a = 1)
  ey11 <- get_marginal_ey(py1, pmmarg = pm_a1, pu1, a = 1)
  return(cbind(ey00 = ey00, ey01 = ey01, ey10 = ey10, ey11 = ey11))
  
}



#' Function to calculate a Bayesian bootstrap weighted mean
#' 
#' @param statistic Length-n vector of a statistic at each covariate pattern
#' @param counts Length-n vector of counts for each covariate pattern
#' @return Length-B vector of population average statistics
#' @export
bayes_boot_weighted_mean <- function(statistic, counts) {
  stopifnot(length(statistic) == length(counts))
  w <- sapply(counts, FUN = function(x) rgamma(n = 1, shape = x, scale = 1))
  w <- w / sum(w)
  return(weighted.mean(statistic, w = w))
}



#' Process SEER data to make design matrices and counts for weighting
#' 
#' @param seer_file_path Path to raw SEER txt file
#' @param am_intx Whether to include A * M interaction in Y model
#' @return Named list of design matrices Xmat_U, ..., Xmat_Y with K unique rows, 
#' frequency counts for each of the K covariate pattern, number of unique 
#' patterns, and number of unique values for the mediator (K_m) 
process_seer_raw <- function(seer_file_path, am_intx = 1) {
  
  # Read in seer from file path, and set 
  seer       <- read.table(path.expand(seer_file_path))
  mediator   <- "stage_cat"
  outcome    <- "fiveyearsurv" 
  tx         <- "black"
  unmeasured <- "pov"
  n          <- NROW(seer)
  K_m        <- length(unique(seer[[mediator]]))
  
  # Make middle age category the reference
  seer$age_cat <- relevel(as.factor(seer$age_cat), ref = 2)
  seer$female  <- ifelse(seer$female == 2, 1, 0)
  
  # Collapse regions to match CanCors
  seer$region[seer$region == 2] <- 1
  seer$region  <- factor(seer$region,
                         labels = c("Other", "South", "West"))
  
  formulas <- list()
  formulas$mediator   <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  formulas$unmeasured <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black)
  if (am_intx) {
    formulas$outcome  <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black + as.factor(stage_cat)
                                 + black*as.factor(stage_cat))
  } else {
    formulas$outcome  <- formula(~ female + as.factor(age_cat) + as.factor(region) 
                                 + black + as.factor(stage_cat))
  }
  
  # These samples gets overwritten - just ensure the model matrices work
  seer[["Y"]] <- seer[[outcome]]
  seer[["M"]] <- seer[[mediator]] <- sample(1:K_m, size = n, replace = TRUE)
  seer[["U"]] <- sample(0:1, size = n, replace = TRUE) # gets overwritten
  seer[["A"]] <- seer[[tx]] <- sample(0:1, size = n, replace = TRUE)
  
  # Construct formulas
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
  
  # Simplify into unique patterns for Bayesian bootstrap
  full_bl <- cbind(count = 1, model.matrix(formulas$baseline, data = seer))
  all_counts <- aggregate(count ~ ., data = full_bl, FUN = sum)
  patts <- unique(full_bl[, -1])
  patts_with_counts <- merge(patts, all_counts, sort = FALSE)
  counts <- patts_with_counts$count
  n_patt <- nrow(patts)
  
  # Make shell design matrices for unique covariate patterns
  Xmat_Y <- matrix(NA,
                   ncol = ncol(model.matrix(formulas$outcomeU, seer)),
                   nrow = nrow(patts))
  Xmat_U <- matrix(NA,
                   ncol = ncol(model.matrix(formulas$unmeasured, seer)),
                   nrow = nrow(patts))
  Xmat_M <- matrix(NA,
                   ncol = ncol(model.matrix(formulas$mediatorU, seer)),
                   nrow = nrow(patts))
  Xmat_Y[1:nrow(patts), 1:ncol(patts)] <- patts
  Xmat_U[1:nrow(patts), 1:ncol(patts)] <- patts
  Xmat_M[1:nrow(patts), 1:ncol(patts)] <- patts
  
  return(list(Xmat_Y = Xmat_Y, Xmat_M = Xmat_M, Xmat_U = Xmat_U, 
              counts = counts, n_patt = n_patt, K_m = K_m))
}



#' Process data application stan fit to obtain posterior draws of mediation estimands
#' 
#' @param seer_file_path Path to SEER data
#' @param stan_fit Stan fit from data application 
#' @param am_intx Whether an exposure-mediator interaction was included
#' @return Named list of mediation quantities from closed form g-formula approach
#' @examples
#' \dontrun{
#' process_data_app_gf(seer_file_path = "../Data/Raw/data_SEER_for_BSA_summer_project.txt",
#'                     stan_fit = readRDS("../data_app_res_20180703.rds"),
#'                     am_intx = 1)
#' }
#' @export
process_data_app_gf <- function(seer_file_path, stan_fit, am_intx) {
  
  # Make design matrix arrays
  dfl <- process_seer_raw(seer_file_path = seer_file_path, am_intx = am_intx)
  xml <- make_xy_xm_xu(Xmat_Y = dfl$Xmat_Y, Xmat_M = dfl$Xmat_M, 
                       Xmat_U = dfl$Xmat_U, am_intx = am_intx, 
                       n_patt = dfl$n_patt, K_m = dfl$K_m)
  counts <- dfl$counts
  K_m <- dfl$K_m
  n_patt <- dfl$n_patt
  rm(dfl)
  xy <- xml$xy
  xm <- xml$xm
  xu <- xml$xu
  rm(xml)
  
  # Extract coefficients from model
  coefs <- as.array(rstan::extract(stan_fit, 
                                   pars = c("coef_y", "coef_m", "coef_u")))
  coef_u_arr <- t(coefs$coef_u)
  coef_m_arr <- aperm(coefs$coef_m, c(3, 2, 1)) # have last dimension be R
  coef_y_arr <- t(coefs$coef_y)
  R <- dim(coef_u_arr)[2]
  
  # Initialize containers
  # i = covariate patter, a = tx level, m = mediator level, u = unmeasured 
  # confounder level
  pu1 <- array(NA, dim = c(n_patt, 2)) #i, a
  pm  <- array(NA, dim = c(n_patt, K_m, 2, 2)) #i, m, a, u
  py1 <- array(NA, dim = c(n_patt, 2, 2, K_m)) #i, a, u, m
  nder_gfs <- nier_gfs <- ter_gfs <- rep(NA, R)
  
  # Loop over MCMC iterations
  for (r in 1:R) {
    coef_u <- coef_u_arr[, r]
    coef_m <- coef_m_arr[, , r]
    coef_y <- coef_y_arr[, r]
    
    for (a in 0:1) {
      pu1[ , a + 1] <- plogis(xu[ , , a + 1] %*% coef_u)
      for (u in 0:1) {
        pm[ , , a + 1, u + 1]  <- row_softmax(xm[ , , a + 1, u + 1] %*% coef_m)
        for (m in 1:4) {
          py1[ , a + 1, u + 1, m] <- plogis(xy[ , , a + 1, u + 1, m] %*% coef_y)
        }
      }
    }
    
    q <- make_quartet(py1, pm, pu1)
    nder_gfs[r] <- bayes_boot_weighted_mean(q[ , "ey10"] - q[ , "ey00"], 
                                            counts = counts)
    nier_gfs[r] <- bayes_boot_weighted_mean(q[ , "ey11"] - q[ , "ey10"],
                                            counts = counts)
    ter_gfs[r]  <- bayes_boot_weighted_mean(q[ , "ey11"] - q[ , "ey00"],
                                            counts = counts)
  }
  
  # Simulated versions
  nder_gcs <- c(unlist(extract(stan_fit, pars = "meffects[1]")))
  nier_gcs <- c(unlist(extract(stan_fit, pars = "meffects[2]")))
  ter_gcs  <- c(unlist(extract(stan_fit, pars = "meffects[3]")))
  nder_gc <- mean(nder_gcs)
  nier_gc <- mean(nier_gcs)
  ter_gc  <- mean(ter_gcs)
  ci_nder_gc <- make_ci(nder_gcs)
  ci_nier_gc <- make_ci(nier_gcs)
  ci_ter_gc  <- make_ci(ter_gcs)
  
  # Posterior summaries
  nder_gf <- mean(nder_gfs)
  nier_gf <- mean(nier_gfs)
  ter_gf  <- mean(ter_gfs)
  ci_nier_gf <- make_ci(nier_gfs)
  ci_nder_gf <- make_ci(nder_gfs)
  ci_ter_gf  <- make_ci(ter_gfs)
  
  return(list(nder_gcs = nder_gcs, nder_gc = nder_gc, ci_nder_gc = ci_nder_gc,
              nier_gcs = nier_gcs, nier_gc = nier_gc, ci_nier_gc = ci_nier_gc,
              ter_gcs  = ter_gcs,  ter_gc  = ter_gc,  ci_ter_gc  = ci_ter_gc,
              nder_gfs = nder_gfs, nder_gf = nder_gf, ci_nder_gf = ci_nder_gf,
              nier_gfs = nier_gfs, nier_gf = nier_gf, ci_nier_gf = ci_nier_gf,
              ter_gfs  = ter_gfs,  ter_gf  = ter_gf,  ci_ter_gf  = ci_ter_gf))
}


