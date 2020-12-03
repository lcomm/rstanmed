functions {
#include /chunks/basic_functions.stan
#include /chunks/bin_bin_mediation_functions.stan
}
data {
  int<lower=0> N; // # of unique observation patterns in big data set
  int<lower=0> P_m; // # of covariates in mediator model
  int<lower=0> P_y; // # of covariates in outcome model
  int<lower=0> P_u; // # of covariates in outcome model
  int<lower=0,upper=1> y[N]; // collapsed observed outcome
  int<lower=0,upper=1> m[N]; // collapsed observed mediator
  matrix[N, P_m - 1] x_m; // collapsed design matrix for mediator (excludes u)
  matrix[N, P_y - 1] x_y; // collapsed design matrix for outcome (excludes u)
  matrix[N, P_u] x_u; // collapsed design matrix for u
  vector[P_m] location_beta; // NOT used for prior, just scaling
  vector[P_y] location_alpha; // NOT used for prior, just scaling
  vector[P_u] location_gamma; // NOT used for prior, just scaling
  cov_matrix[P_m] scale_V_beta; // NOT used for prior, just scaling
  cov_matrix[P_y] scale_V_alpha; // NOT used for prior, just scaling
  cov_matrix[P_u] scale_V_gamma; // NOT used for prior, just scaling
  int<lower=0> N_s; // # of unique observation patterns in small data set
  int<lower=0,upper=1> small_y[N_s]; // collapsed small observed outcome
  int<lower=0,upper=1> small_m[N_s]; // collapsed small observed mediator
  int<lower=0,upper=1> small_u[N_s]; // collapsed small observed mediator
  matrix[N_s, P_m] small_x_m; // collapsed design small matrix for mediator (includes u)
  matrix[N_s, P_y] small_x_y; // collapsed design small matrix for outcome (includes u)
  matrix[N_s, P_u] small_x_u; // collapsed design small matrix for u
  int u_ei; // whether or not u is exposure-induced
  int am_intx; // whether there is exposure-mediator interaction
  int mean_only; // whether to do means of final counterfactual
  int w[N]; //count pattern for frequency of observations
  int small_w[N_s]; //count pattern for frequency of observations
  real<lower=0,upper=1> pp_alpha; // downscaling factor for power prior (1 = simple data pooling)
} 
transformed data {
  matrix[P_m, P_m] beta_L;
  matrix[P_y, P_y] alpha_L;
  matrix[P_u, P_u] gamma_L;
  int N_tot = sum(w);
  vector[N] w_vec = to_vector(w);
  vector[N_s] small_w_vec = to_vector(small_w);
  // uncollapsed (i.e., non-count) versions of design matrices
  // used only in generated quantities block
  matrix[N_tot, cols(x_y)] big_x_y = uncollapse_matrix(x_y, w);
  matrix[N_tot, cols(x_m)] big_x_m = uncollapse_matrix(x_m, w);
  matrix[N_tot, cols(x_u)] big_x_u = uncollapse_matrix(x_u, w);
  beta_L  = cholesky_decompose(scale_V_beta);
  alpha_L = cholesky_decompose(scale_V_alpha);
  gamma_L = cholesky_decompose(scale_V_gamma);
}
parameters {
  vector[P_m] beta_unscaled;
  vector[P_y] alpha_unscaled;
  vector[P_u] gamma_unscaled;
}
transformed parameters {
  vector[P_m] beta = location_beta + beta_L * beta_unscaled;
  vector[P_y] alpha = location_alpha + alpha_L * alpha_unscaled;
  vector[P_u] gamma = location_gamma + gamma_L * gamma_unscaled;
  vector[P_y - 1] alpha_no_u = head(alpha, P_y - 1);
  real alpha_u = alpha[P_y];
  vector[P_m - 1] beta_no_u = head(beta, P_m - 1);
  real beta_u = beta[P_m];
}
model {
  vector[N] ll_marg;
  real pu1; 
  real ll_if_u1;
  real ll_if_u0;
  // vector[N_s] small_ll;
  
  // likelihood for small data set (i.e., data pooling / power prior)
  for (n in 1:N_s) {
    // small_ll[n] = bernoulli_logit_lpmf(small_u[n] | small_x_u[n, ] * gamma) + 
    //               bernoulli_logit_lpmf(small_y[n] | small_x_y[n, ] * alpha) + 
    //               bernoulli_logit_lpmf(small_m[n] | small_x_m[n, ] * beta);
    target += pp_alpha * small_w_vec[n] * (
                  bernoulli_logit_lpmf(small_u[n] | small_x_u[n, ] * gamma) + 
                  bernoulli_logit_lpmf(small_y[n] | small_x_y[n, ] * alpha) + 
                  bernoulli_logit_lpmf(small_m[n] | small_x_m[n, ] * beta));
  }
  // target += pp_alpha * dot_product(small_w_vec, small_ll);
  
  // likelihood for big data set
  for (n in 1:N) {
    pu1 = inv_logit(x_u[n, ] * gamma);
    ll_if_u0 = bernoulli_logit_lpmf(y[n] | x_y[n, ] * alpha_no_u) + 
               bernoulli_logit_lpmf(m[n] | x_m[n, ] * beta_no_u);
    ll_if_u1 = bernoulli_logit_lpmf(y[n] | x_y[n, ] * alpha_no_u + alpha_u) + 
               bernoulli_logit_lpmf(m[n] | x_m[n, ] * beta_no_u + beta_u);
    ll_marg[n] = log_mix(pu1, ll_if_u1, ll_if_u0);
  }
  target += dot_product(w_vec, ll_marg);
  
}
generated quantities {

  vector[3] meffects;
  meffects = bbootColMeans_rng(sim_ceffects_rng(alpha, beta, gamma, u_ei, 
                                                big_x_y, big_x_m, big_x_u, 
                                                am_intx, mean_only));
}
