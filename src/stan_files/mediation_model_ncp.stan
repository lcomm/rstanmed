functions {
#include /chunks/basic_functions.stan
#include /chunks/bin_cat_mediation_functions.stan
}
data {
    // CONSTANTS
    int<lower=2> K_m;                   // number of categories for M
    int<lower=2> D_m;                   // design matrix dimension, M model
    int<lower=2> D_y;                   // design matrix dimension, Y model
    int<lower=2> D_u;                   // design matrix dimension, U model
    int<lower=1> N;                     // number of unique observation types
    int<lower=0,upper=1> am_intx;       // 0/1 exposure-mediator interaction
    // DESIGN MATRICES
    matrix[N, D_m] x_mu1;            // design matrix, M model if U = 1
    matrix[N, D_m] x_mu0;            // design matrix, M model if U = 0
    matrix[N, D_y] x_yu1;            // design matrix, Y model if U = 1
    matrix[N, D_y] x_yu0;            // design matrix, Y model if U = 0
    matrix[N, D_u] x_u;              // design matrix, U model
    // R INVERSE MATRICES;
    matrix[D_m*(K_m-1), D_m*(K_m-1)] Rstar_inv_M;
    matrix[D_y, D_y] Rstar_inv_Y;
    matrix[D_u, D_u] Rstar_inv_U;
    // REGRESSION PARAMETERS
    vector[D_u] m_coef_u;               // mean, coef_u prior
    cov_matrix[D_u] v_coef_u;           // vcov, coef_u prior
    vector[D_m*(K_m-1)] m_coef_m;       // mean, coef_m prior
    cov_matrix[D_m*(K_m-1)] v_coef_m;   // vcov, coef_m prior
    vector[D_y] m_coef_y;               // mean, coef_y prior
    cov_matrix[D_y] v_coef_y;           // vcov, coef_y prior
    // OUTCOME DATA
    int<lower=1, upper=K_m> m[N];       // observations of mediator M
    int<lower=0, upper=1> y[N];         // observations of outcome Y
    // WEIGHTS
    int w[N];                        // number of times data pattern is observed
    // SIMULATION
    int mean_only;
    
}
transformed data {
    row_vector[D_m] zs = rep_row_vector(0, D_m);
    matrix[D_m*(K_m-1), D_m*(K_m-1)] beta_L  = cholesky_decompose(v_coef_m);
    matrix[D_y, D_y] alpha_L = cholesky_decompose(v_coef_y);
    matrix[D_u, D_u] gamma_L = cholesky_decompose(v_coef_u);
    int u_ei = 1;
    matrix[N, D_m - 1] x_m = x_mu1[, 1:(D_m - 1)];
    matrix[N, D_y - 1] x_y = x_yu1[, 1:(D_y - 1)];
    int N_tot = sum(w);
    vector[N] w_vec = to_vector(w);
    // uncollapsed (i.e., non-count) versions of design matrices
    matrix[N_tot, cols(x_y)] big_x_y = uncollapse_matrix(x_y, w);
    matrix[N_tot, cols(x_m)] big_x_m = uncollapse_matrix(x_m, w);
    matrix[N_tot, cols(x_u)] big_x_u = uncollapse_matrix(x_u, w);
}
parameters {
    vector[D_y] coef_y_unscaled;
    vector[D_u] coef_u_unscaled;
    vector[(K_m - 1)*D_m] coef_m_raw_unscaled; // non-zero parameters
}
transformed parameters {
    // coef_X is on QR-decomposed scale
    // greek parameter versions are on original (interpretable) scale
    vector[D_y] coef_y = m_coef_y + alpha_L * coef_y_unscaled;
    vector[D_u] coef_u = m_coef_u + gamma_L * coef_u_unscaled;
    vector[(K_m - 1)*D_m] coef_m_raw_scaled = m_coef_m + 
                                              beta_L * coef_m_raw_unscaled;
    matrix[K_m, D_m] coef_m = append_row(zs, 
                              to_matrix(coef_m_raw_scaled, (K_m - 1), D_m, 1));
    matrix[K_m, D_m] beta = append_row(zs, 
                                       to_matrix(Rstar_inv_M * coef_m_raw_scaled, 
                                                 (K_m - 1), D_m, 1));
    vector[D_y] alpha = Rstar_inv_Y * coef_y;
    vector[D_u] gamma = Rstar_inv_U * coef_u;
}
model {
    /** Declarations **/
    // Linear predictors
    // P(U = 1 | design matrix for U regression)
    vector[N] eta_u; 
    
    // P(M = m | U = 1, design matrix for M regression)
    matrix[K_m, N] eta_mu1;
    
    // P(M = m | U = 1, design matrix for M regression)
    matrix[K_m, N] eta_mu0;

    // P(Y = 1 | U = 1, design matrix for Y regression)
    vector[N] eta_yu1;
    
    // P(Y = 1 | U = 0, design matrix for Y regression)
    vector[N] eta_yu0;
    
    // Log-likelihood contributions
    vector[N] ll_marg;
    
    /** Calculations **/
    eta_u = x_u * coef_u;
    eta_mu1 = coef_m * x_mu1';
    eta_mu0 = coef_m * x_mu0';
    eta_yu1 = x_yu1 * coef_y;
    eta_yu0 = x_yu0 * coef_y;

    for (n in 1:N) {
        real x0;
        real x1;

        if (y[n] == 0){
            // u = 0 portion
            x0 = log1m_inv_logit(eta_yu0[n]) + 
                 log_softmax(col(eta_mu0, n))[m[n]] + 
                 log1m_inv_logit(eta_u[n]);
                 
            // u = 1 portion
            x1 = log1m_inv_logit(eta_yu1[n]) + 
                 log_softmax(col(eta_mu1, n))[m[n]] + 
                 log_inv_logit(eta_u[n]);

        } else {

            // u = 0 portion
            x0 = log_inv_logit(eta_yu0[n]) + 
                 log_softmax(col(eta_mu0, n))[m[n]] + 
                 log1m_inv_logit(eta_u[n]);
            
            // u = 1 portion
            x1 = log_inv_logit(eta_yu1[n]) + 
                 log_softmax(col(eta_mu1, n))[m[n]] + 
                 log_inv_logit(eta_u[n]);
    
        }

        // joint (M,Y) log-probabilities, marginal over U
        ll_marg[n] = log_sum_exp(x0, x1);
    }
    
    // Add likelihood to target
    target += dot_product(w_vec, ll_marg);
    
    /** Priors **/
    coef_y_unscaled ~ normal(0, 1);
    coef_u_unscaled ~ normal(0, 1);
    coef_m_raw_unscaled ~ normal(0, 1);
    
}
generated quantities {

  vector[3] meffects;
  meffects = bbootColMeans_rng(sim_ceffects_mcat_rng(alpha, beta, gamma, u_ei, 
                                                     big_x_y, big_x_m, big_x_u, 
                                                     am_intx, K_m, mean_only));
}
