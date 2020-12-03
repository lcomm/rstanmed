#' Make short variables for frequentist corrections
#' 
#' @param df Data application data frame
#' @return Data frame with additional aliased variables
#' @export
make_da_shortcut_vars <- function(df) {
  df$z1 <- df$female
  df$z2 <- as.factor(df$age_cat)
  df$z3 <- as.factor(df$region)
  df$a <- df$black
  df$m <- df$stage_cat
  df$y <- df$fiveyearsurv
  if ("pov" %in% colnames(df)) {
    df$u <- df$pov
  }
  return(df)
}



#' Process SEER for frequentist corrections
#' 
#' @param seer_file_path File path to raw SEER data
#' @return SEER data with aliased variables z1, z2, ...
#' @export
make_da_big_df <- function(seer_file_path) {
  seer <- read.table(seer_file_path)
  
  # Make middle age category the reference
  seer$age_cat <- relevel(as.factor(seer$age_cat), ref = 2)
  seer$female  <- ifelse(seer$female == 2, 1, 0)
  
  # Collapse regions to match CanCors
  seer$region[seer$region == 2] <- 1
  seer$region  <- factor(seer$region,
                         labels = c("Other", "South", "West"))
  
  seer <- make_da_shortcut_vars(seer)
  return(seer)
}



#' Process CanCORS for frequentist corrections
#' 
#' @param cancors_file_path File path to raw CanCORS data
#' @return CanCORS data with aliased variables z1, z2, ...
#' @export
make_da_small_df <- function(cancors_file_path) {
  # Read in raw data
  cancors <- read.csv(cancors_file_path)
  
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
  
  cancors <- make_da_shortcut_vars(cancors)
  return(cancors)
}



#' Count up baseline covariate patterns for predictions
#' 
#' @param df Data frame containing z1, z2, z3
#' @return Data frame containing unique combinations of z1, z2, z2, and frequency counts
#' @export
make_da_cov_df_with_counts <- function(df) {
  cov_df <- df[, c("z1", "z2", "z3")]
  cov_df <- dplyr::group_by(.data = cov_df, z1, z2, z3)
  cov_df <- as.data.frame(dplyr::tally(x = cov_df))
  colnames(cov_df)[colnames(cov_df) == "n"] <- "counts"
  return(cov_df)  
}


#' Data application interaction correction
#' 
#' @param df Big data frame containing y, a, m, z1, z2, z3
#' @param small_df Small data containing y, a, m, z1, z2, z3, u
#' @param nonparametric Whether to calculate everything nonparametrically
#' @return Estimated bias-corrected population (r)NDE
#' @export
run_da_intx_corr <- function(df, small_df, nonparametric = TRUE) {
  
  if (nonparametric) {
    formula_y <- formula(y ~ 1 + (z1 + z2 + z3 + a + m + u)^6)
  } else {
    formula_y <- formula(y ~ 1 + z1 + z2 + z3 + a + m + a:m + u)
  }
  
  # Big data sample frequencies of covariate patterns
  cov_df <- make_da_cov_df_with_counts(df)
  
  a <- 1
  a_star <- 0
  
  # Fit models to calculate B_nde
  fit_y <- glm(y ~ 1 + (z1 + z2 + z3 + a + m + u)^6, data = small_df, 
               family = binomial(link = "logit"))
  fit_u_cond <- glm(u ~ 1 + (z1 + z2 + z3 + a + m)^5, data = small_df, 
                    family = binomial(link = "logit")) # weird but correct
  fit_u_marg <- glm(u ~ 1 + (z1 + z2 + z3)^3, data = small_df, 
                    family = binomial(link = "logit"))
  
  # Calculate bias for every covariate pattern, summing across m and u
  cov_df$B_nde <- 0
  for (m in 1:4) {
    small_df$meqm <- (small_df$m == m) * 1
    fit_m_cond <- glm(meqm ~ 1 + (z1 + z2 + z3 + a + u)^5, data = small_df, 
                      family = binomial(link = "logit"))
    fit_m_marg <- glm(meqm ~ 1 + (z1 + z2 + z3 + a)^4, data = small_df, 
                      family = binomial(link = "logit"))
    
    for (u in 0:1) {
      ey1 <- predict(fit_y, newdata = cbind(cov_df, a = a, m = m, u = u), 
                     type = "response")
      ey2 <- predict(fit_y, newdata = cbind(cov_df, a = a_star, m = m, u = u), 
                     type = "response")
      pm3 <- predict(fit_m_marg, newdata = cbind(cov_df, a = a_star), 
                     type = "response")
      pm4 <- predict(fit_m_cond, newdata = cbind(cov_df, a = a_star, u = u), 
                     type = "response")
      qu5 <- predict(fit_u_cond, newdata = cbind(cov_df, a = a, m = m), 
                     type = "response")
      pu5 <- u * qu5 + (1 - u) * (1 - qu5)
      qu6 <- predict(fit_u_cond, newdata = cbind(cov_df, a = a_star, m = m), 
                     type = "response")
      pu6 <- u * qu6 + (1 - u) * (1 - qu6)
      qu7 <- predict(fit_u_marg, newdata = cbind(cov_df), 
                     type = "response")
      pu7 <- u * qu7 + (1 - u) * (1 - qu7)
      cov_df$B_nde <- cov_df$B_nde + 
        (ey1 * pu5 * pm3 
         - ey2 * pu6 * pm3
         - ey1 * pm4 * pu7 
         + ey2 * pm4 * pu7)
    }
  }
  
  naive_nder <- calculate_categorical_naive_nder(df = df, 
                                                 nonparametric = nonparametric)
  
  # Subtract off bias and average by covariate pattern frequency
  intx_nder <- naive_nder - cov_df$B_nde
  avg_intx_nder <- weighted.mean(intx_nder, w = cov_df$counts)
  return(avg_intx_nder)
}


#' Calculate the naive NDER for a categorical mediator
#' 
#' @param df Data frame for calculationg
#' @param nonparametric Whether to fit nearly-saturated models (i.e., do NP 
#' estimation). Note: this is not actually nonparametric because we need a little
#' bit of smoothing due to empty cells
#' @param avg Whether to take weighted average 
#' @return Vector (if avg = FALSE) of naive NDEs for each covariate pattern or scalar
#' NDE
#' @export
calculate_categorical_naive_nder <- function(df, nonparametric = FALSE, avg = FALSE) {
  a <- 1
  astar <- 0
  
  cov_df <- make_da_cov_df_with_counts(df)
  cov_df$NDE <- 0
  
  if (nonparametric) {
    fit_y <- lm(y ~ 1 + (z1 + z2 + z3 + a + m)^4, data = df)
  } else {
    fit_y <- glm(y ~ 1 + z1 + z2 + z3 + a * m, data = df,
                 family = binomial(link = "logit"))
  }
  for (m in 1:4) {
    df$meqm <- (df$m == m) * 1
    if (nonparametric) {
      fit_meqm <- glm(meqm ~ 1 + (z1 + z2 + z3 + a)^2, data = df,
                      family = binomial(link = "logit"))
    } else {
      fit_meqm <- glm(meqm ~ 1 + z1 + z2 + z3 + a, data = df,
                      family = binomial(link = "logit"))
    }
    # browser()
    pm <- predict(fit_meqm, newdata = cbind(cov_df, a = astar), type = "response")
    ey1 <- predict(fit_y, newdata = cbind(cov_df, m = m, a = a), type = "response")
    ey2 <- predict(fit_y, newdata = cbind(cov_df, m = m, a = astar), type = "response")
    cov_df$NDE <- cov_df$NDE + (ey1 - ey2) * pm
  }
  
  if (avg) {
    naive <- weighted.mean(x = cov_df$NDE, w = cov_df$counts)
  } else {
    naive <- cov_df$NDE
  }
  return(naive)
}

# seer <- make_da_big_df(seer_file_path = "~/Dropbox/mediation-BSA/Data/Raw/data_SEER_for_BSA_summer_project.txt")
# cancors <- make_da_small_df(cancors_file_path = "~/Dropbox/mediation-BSA/Data/Raw/selected_cancors_data_2016_3_14.csv")

# calculate_categorical_naive_nder(df = seer, nonparametric = TRUE, avg = TRUE)
# calculate_categorical_naive_nder(df = cancors, nonparametric = FALSE, avg = TRUE)
# calculate_categorical_naive_nder(df = cancors, nonparametric = TRUE, avg = TRUE)

# run_da_intx_corr(df = seer, small_df = cancors, nonparametric = FALSE)
# run_da_intx_corr(df = seer, small_df = cancors, nonparametric = TRUE)

