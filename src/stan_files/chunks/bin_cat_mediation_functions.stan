  /**
  * Simulate unmeasured confounder under counterfactual scenario A = a
  * (Assumes U model is bernoulli with logit link)
  *
  * @param gamma Length-P Vector of regression parameters from U model. If U is
  * exposure-induced then the A coefficient will be in the Pth position.
  * @param x_u NxP design matrix for U regression
  * @param a 0/1 level at which to set A
  * @param u_ei 0/1 flag for whether U is exposure-induced
  * return Length-N vector of simulated values for U_a
  **/
  vector ua_rng(vector gamma, matrix x_u, int a, int u_ei) {
    int N = rows(x_u);
    int P_u = cols(x_u);
    vector[N] ua;
    vector[N] lp;
    vector[P_u] gamma_zeroed = gamma;
    gamma_zeroed[P_u] = 0;
    
    if (u_ei == 0) {
      lp = x_u * gamma;
    } else if (u_ei == 1) {
      lp = x_u * gamma_zeroed + a * gamma[P_u];
    }
    
    for (n in 1:N) {
      ua[n] = bernoulli_logit_rng(lp[n]);
    }
    
    return ua;
  }
  
  /**
  * Simulate mediator under counterfactual scenario A = a
  * (Assumes M model is baseline category logit)
  * 
  * @param beta K_m x P matrix of regression parameters from M model. Column
  * (P - 1) corresponds to A and column P corresponds to U. Row 1 is reference level
  * and contains all zeros
  * @param x_m N x (P - 1) design matrix for M regression
  * @param u Length-N vector of U values (probably simulated)
  * @param a 0/1 level at which to set A
  * return Length-N vector of simulated values for M_a
  **/
  int[] ma_cat_rng(matrix beta, matrix x_m, vector u, int a) {
    //TODO(LCOMM): update to categorical (done but not checked)
    int N = rows(x_m);
    int P_m = cols(beta);
    int K_m = rows(beta);
    int ma[N];
    matrix[K_m, N] lp00;
    matrix[K_m, P_m - 1] beta_zeroed = beta[ , 1:(P_m - 1)];
    beta_zeroed[ , P_m - 1] = rep_vector(0, K_m);
    
    // this is lp if a = 0 and u = 0
    lp00 = beta_zeroed * x_m';
    
    for (n in 1:N) {
      ma[n] = categorical_logit_rng(lp00[ , n] + a * beta[ , P_m - 1] + 
                                    u[n] * beta[ , P_m]);
    }
    
    return ma;
  }
  
  /**
  * Simulate outcome under counterfactual scenario A = a
  * (Assumes Y model is bernoulli with logit link)
  * 
  * @param alpha Length-P Vector of regression parameters from Y model. Element
  * (P - 2) corresponds to A, element (P - 1) to M, and element P to U.
  * @param x_y NxP design matrix for Y regression
  * @param m Length-N vector of M values (probably simulated)
  * @param u Length-N vector of U values (probably simulated)
  * @param a 0/1 level at which to set A
  * @param am_intx 0/1 indicator for exposure-mediator interaction
  * @param K_m number of levels for M
  * @param mean_only 0/1 indicator for whether expected values of the potential
  * outcomes should be returned (1) or simulated values (0)
  * return Length-N vector of simulated values for Y_a
  **/
  vector ya_mcat_rng(vector alpha, matrix x_y, int[] m, vector u, int a, 
                     int am_intx, int K_m, int mean_only) {
    //TODO(LCOMM): change for categorical
    int N = rows(x_y);
    int P_y = cols(x_y) + 1;
    vector[N] ya;
    vector[N] lp;
    vector[P_y - 1] alpha_zeroed = alpha[1:(P_y - 1)];
    int a_index = am_intx ? (P_y - 2 * K_m + 1) : (P_y - K_m);
    real alpha_a = alpha[a_index];
    real alpha_u = alpha[P_y];
    alpha_zeroed[a_index:(P_y - 1)] = rep_vector(0, (P_y - a_index));
    
    lp = x_y * alpha_zeroed + a * alpha_a + u * alpha_u; // still need m and a * m
    
    for (n in 1:N) {
      // add components for m and a * m interaction only if necessary
      if (m[n] != 1) {
        lp[n] += alpha[a_index + m[n] - 1];
        if (am_intx == 1 && a == 1) {
          lp[n] += alpha[a_index + (K_m - 1) + m[n] - 1];
        }
      }
      // return mean for y or draw from y
      if (mean_only == 1) {
        ya[n] = inv_logit(lp[n]);  
      } else {
        ya[n] = 1.0 * bernoulli_logit_rng(lp[n]); 
      }
    } // end loop over n
    return ya;
  }
  
  /**
  * Simulate quartet of potential Y_(a, M_a*) for all (a, a*) combinations
  * 
  * @param alpha Regression coefficients from Y model. See \code{\link{ya_rng}}
  * for details on coefficient order.
  * @param beta Regression coefficients from M model. See \code{\link{ma_rng}}
  * for details on coefficient order.
  * @param gamma Regression coefficients from A model. See \code{\link{ua_rng}}
  * for details on coefficient order.
  * @param u_ei 0/1 indicator for whether U is exposure-induced. See 
  * \code{\link{ua_rng}} for details.
  * @param x_y Design matrix with N rows for Y regression model.
  * @param x_m Design matrix with N rows for Y regression model.
  * @param x_u Design matrix with N rows for U regression model.
  * @param am_intx 
  * @param K_m number of levels for M
  * @param mean_only 0/1 indicator for whether expected values of the potential
  * outcomes should be returned (1) or simulated values (0)
  * @return Nx4 matrix of means or simulated values
  * @export
  **/
  matrix quartet_mcat_rng(vector alpha, matrix beta, vector gamma, int u_ei,
                          matrix x_y, matrix x_m, matrix x_u, 
                          int am_intx, int K_m, int mean_only) {
    
    int N = rows(x_y);
    vector[N] u0;
    vector[N] u1;
    int m0[N];
    int m1[N];
    matrix[N, 4] quartet;
    
    // precursors
    u0 = ua_rng(gamma, x_u, 0, u_ei);
    if (u_ei == 1) {
      u1 = ua_rng(gamma, x_u, 1, u_ei);
    } else {
      u1 = u0;
    }
    m0 = ma_cat_rng(beta, x_m, u0, 0);
    m1 = ma_cat_rng(beta, x_m, u1, 1);
    
    // 0, 0
    quartet[,1] = ya_mcat_rng(alpha, x_y, m0, u0, 0, am_intx, K_m, mean_only);
    
    // 0, 1
    // recanting witness
    quartet[,2] = ya_mcat_rng(alpha, x_y, m1, u0, 0, am_intx, K_m, mean_only); 
    
    // 1, 0
    // recanting witness
    quartet[,3] = ya_mcat_rng(alpha, x_y, m0, u1, 1, am_intx, K_m, mean_only);
    
    // 1, 1
    quartet[,4] = ya_mcat_rng(alpha, x_y, m1, u1, 1, am_intx, K_m, mean_only);

    return quartet;
  }
  
  
  // Calculate simulation-based conditional NDER (Y_(1, g=M_0) - Y_(0, g=M_0)), 
  // NIER (Y_(1, g=M_1) - Y_(1, g=M_0)), TER (Y_(1, g=M_1) - Y_(0, g=M_0)
  // 
  // @param alpha Vector of regression coefficients from Y model.
  // See \code{\link{ya_rng}} for details on coefficient order. 
  // @param beta Matrix of regression coefficients from M model. 
  // See \code{\link{ma_rng}} for details on coefficient order.
  // @param gamma Vector of regression coefficients from A model.
  // See \code{\link{ua_rng}} for details on coefficient order.
  // @param u_ei 0/1 indicator for whether U is exposure-induced. See 
  // \code{\link{ua_rng}} for details.
  // @param x_y Design matrix with N rows for Y regression model 
  // (or fewer rows, if weighted).
  // @param x_m Design matrix with N rows for Y regression model
  // (or fewer rows, if weighted).
  // @param x_u Design matrix with N rows for U regression model
  // @param mean_only 0/1 indicator for whether expected values of the potential
  // outcomes should be returned (1) or simulated values (0)
  // @param N Number of unique observation types
  // @param N_tot Number of total observations 
  // @param w Vector of weights that sums to N_tot
  // @return N_tot x 3 matrix of means or differences in simulated values. 
  // Columns are NDER, NIER, and TER, rows are for bootstrapped observations.
  // @export
  matrix sim_ceffects_mcat_rng(vector alpha, matrix beta, vector gamma, int u_ei,
                               matrix x_y, matrix x_m, matrix x_u, 
                               int am_intx, int K_m, int mean_only, 
                               int N, int N_tot,
                               vector w) {
  vector[N] boot_probs = w ./ sum(w);
  int row_i;
  matrix[N_tot, 3] ceffects;
  matrix[N_tot, 4] quartet;
  matrix[N_tot, cols(x_y)] x_ynew;
  matrix[N_tot, cols(x_m)] x_mnew;
  matrix[N_tot, cols(x_u)] x_unew;
  
  // make vector of indices for the bootstrap sample
  // make "new design matrices"
  for (n in 1:N_tot) {
    row_i = categorical_rng(boot_probs);
    x_ynew[n,] = x_y[row_i,];
    x_mnew[n,] = x_m[row_i,];
    x_unew[n,] = x_u[row_i,];
  }
  
  // get counterfactuals (or mean of counterfactuals)  
  quartet = quartet_mcat_rng(alpha, beta, gamma, u_ei, x_ynew, x_mnew, x_unew, 
                             am_intx, K_m, mean_only);
  ceffects[,1] = quartet[,3] - quartet[,1]; // direct effect
  ceffects[,2] = quartet[,4] - quartet[,3]; // indirect effect
  ceffects[,3] = quartet[,4] - quartet[,1]; // total effect
  
  return ceffects;
  }

