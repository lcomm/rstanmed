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
  * (Assumes M model is bernoulli with logit link)
  * 
  * @param beta Length-P Vector of regression parameters from M model. Element
  * (P - 1) corresponds to A and element P corresponds to U.
  * @param x_m Nx(P - 1) design matrix for M regression
  * @param u Length-N vector of U values (probably simulated)
  * @param a 0/1 level at which to set A
  * return Length-N vector of simulated values for M_a
  **/
  vector ma_rng(vector beta, matrix x_m, vector u, int a) {
    int N = rows(x_m);
    int P_m = cols(x_m) + 1;
    vector[N] ma;
    vector[N] lp;
    vector[P_m - 1] beta_zeroed = beta[1:(P_m - 1)];
    beta_zeroed[P_m - 1] = 0;
    
    lp = x_m * beta_zeroed + a * beta[P_m - 1] + u * beta[P_m];
    
    for (n in 1:N) {
      ma[n] = bernoulli_logit_rng(lp[n]);
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
  * @param mean_only 0/1 indicator for whether expected values of the potential
  * outcomes should be returned (1) or simulated values (0)
  * return Length-N vector of simulated values for Y_a
  **/
  vector ya_rng(vector alpha, matrix x_y, vector m, vector u, int a, 
                int mean_only) {
    int N = rows(x_y);
    int P_y = cols(x_y) + 1;
    vector[N] ya;
    vector[N] lp;
    vector[P_y - 1] alpha_zeroed = alpha[1:(P_y - 1)];
    alpha_zeroed[(P_y - 2):(P_y - 1)] = rep_vector(0, 2);
    
    lp = x_y * alpha_zeroed + a * alpha[P_y - 2] + m * alpha[P_y - 1] + 
         u * alpha[P_y];
    
    if (mean_only == 1) {
      ya = inv_logit(lp);
    } else if (mean_only == 0) {
      for (n in 1:N) {
        ya[n] = bernoulli_logit_rng(lp[n]);
      }  
    }
    
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
  * @param mean_only 0/1 indicator for whether expected values of the potential
  * outcomes should be returned (1) or simulated values (0)
  * @return Nx4 matrix of means or simulated values
  * @export
  **/
  matrix quartet_rng(vector alpha, vector beta, vector gamma, int u_ei,
                     matrix x_y, matrix x_m, matrix x_u, int mean_only) {
    
    int N = rows(x_y);
    vector[N] u0;
    vector[N] u1;
    vector[N] m0;
    vector[N] m1;
    matrix[N, 4] quartet;
    
    // precursors
    u0 = ua_rng(gamma, x_u, 0, u_ei);
    u1 = ua_rng(gamma, x_u, 1, u_ei);
    m0 = ma_rng(beta, x_m, u0, 0);
    m1 = ma_rng(beta, x_m, u1, 1);
    
    // 0, 0
    quartet[,1] = ya_rng(alpha, x_y, m0, u0, 0, mean_only);
    
    // 0, 1
    quartet[,2] = ya_rng(alpha, x_y, m1, u0, 0, mean_only); // recanting witness
    
    // 1, 0
    quartet[,3] = ya_rng(alpha, x_y, m0, u1, 1, mean_only); // recanting witness
    
    // 1, 1
    quartet[,4] = ya_rng(alpha, x_y, m1, u1, 1, mean_only);

    return quartet;
  }
  
  /**
  * Calculate simulation-based NDER for Y_(1, M_0) - Y_(0, M_0) from matrix block
  * of B parameter draws
  * 
  * @param alpha Matrix with B rows containing regression coefficients from Y model.
  * See \code{\link{ya_rng}} for details on coefficient order. 
  * @param beta Matrix with B rows containing regression coefficients from M model. 
  * See \code{\link{ma_rng}} for details on coefficient order.
  * @param gamma Matrix with B rows containing regression coefficients from A model.
  * See \code{\link{ua_rng}} for details on coefficient order.
  * @param u_ei 0/1 indicator for whether U is exposure-induced. See 
  * \code{\link{ua_rng}} for details.
  * @param x_y Design matrix with N rows for Y regression model.
  * @param x_m Design matrix with N rows for Y regression model.
  * @param x_u Design matrix with N rows for U regression model.
  * @param mean_only 0/1 indicator for whether expected values of the potential
  * outcomes should be returned (1) or simulated values (0)
  * @return B x N matrix of means or simulated values
  * @export
  **/
  matrix sim_nder_rng(matrix alpha, matrix beta, matrix gamma, int u_ei,
                      matrix x_y, matrix x_m, matrix x_u, int mean_only) {
  // TODO: Fix for mean_only = FALSE case
  // Need better bootstrap sampling approach first
  int B = rows(alpha);
  int N = rows(x_y);
  matrix[cols(alpha), B] alphat = alpha';
  matrix[cols(beta), B] betat = beta';
  matrix[cols(gamma), B] gammat = gamma';
  vector[N] boot_probs = rep_vector(1.0/N, N);
  int row_i;
  matrix[N, B] nder;
  matrix[N, 4] quartet;
  matrix[N, cols(x_y)] x_ynew;
  matrix[N, cols(x_m)] x_mnew;
  matrix[N, cols(x_u)] x_unew;
  
  // loop over MCMC iterations
  for (b in 1:B) {
    // make vector of indices for the bootstrap sample
    for (n in 1:N) {
      row_i = categorical_rng(boot_probs);
      x_ynew[n,] = x_y[row_i,];
      x_mnew[n,] = x_m[row_i,];
      x_unew[n,] = x_u[row_i,];
    }
    
    quartet = quartet_rng(alphat[,b], betat[,b], gammat[,b], u_ei,
                          x_ynew, x_mnew, x_unew, mean_only);
    nder[,b] = quartet[,3] - quartet[,1];
  }
  
  return nder';
  }
