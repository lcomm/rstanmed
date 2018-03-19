functions {
#include /chunks/basic_functions.stan
#include /chunks/bin_bin_mediation_functions.stan
}
data {
  int<lower=0> N; // # obs
  int<lower=0> P_m; // # of covariates in mediator model
  int<lower=0> P_y; // # of covariates in outcome model
  int<lower=0> P_u; // # of covariates in outcome model
  int<lower=0,upper=1> y[N]; // observed outcome
  int<lower=0,upper=1> m[N]; // observed mediator
  matrix[N, P_m - 1] x_m; // design matrix for mediator (excludes u)
  matrix[N, P_y - 1] x_y; // design matrix for outcome (excludes u)
  matrix[N, P_u] x_u; // design matrix for u
  vector[P_m] prior_mean_beta;
  vector[P_y] prior_mean_alpha;
  vector[P_u] prior_mean_gamma;
  cov_matrix[P_m] prior_vcov_beta;
  cov_matrix[P_y] prior_vcov_alpha;
  cov_matrix[P_u] prior_vcov_gamma;
  int u_ei; // whether or not u is exposure-induced
  int am_intx; // whether there is exposure-mediator interaction
  int mean_only; // whether to do means of final counterfactual
} 
parameters {
  vector[P_m] beta;
  vector[P_y] alpha;
  vector[P_u] gamma;
}
transformed parameters {
  vector[P_y - 1] alpha_no_u = head(alpha, P_y - 1);
  real alpha_u = alpha[P_y];
  vector[P_m - 1] beta_no_u = head(beta, P_m - 1);
  real beta_u = beta[P_m];
}
model {
  vector[N] lp;
  real pu1; 
  real ll_if_u1;
  real ll_if_u0;
  
  // priors
  target += multi_normal_lpdf(alpha | prior_mean_alpha, prior_vcov_alpha);
  target += multi_normal_lpdf(beta | prior_mean_beta, prior_vcov_beta);
  target += multi_normal_lpdf(gamma | prior_mean_gamma, prior_vcov_gamma);
  
  // likelihood
  for (n in 1:N) {
    pu1 = inv_logit(x_u[n, ] * gamma);
    ll_if_u0 = bernoulli_logit_lpmf(y[n] | x_y[n, ] * alpha_no_u) + 
               bernoulli_logit_lpmf(m[n] | x_m[n, ] * beta_no_u);
    ll_if_u1 = bernoulli_logit_lpmf(y[n] | x_y[n, ] * alpha_no_u + alpha_u) + 
               bernoulli_logit_lpmf(m[n] | x_m[n, ] * beta_no_u + beta_u);
    lp[n] = log_mix(pu1, ll_if_u1, ll_if_u0);
    target += lp[n];
  }
}
generated quantities {
  matrix[N, 3] ceffects;
  vector[3] meffects;
  ceffects = sim_ceffects_rng(alpha, beta, gamma, u_ei, x_y, x_m, x_u, 
             am_intx, mean_only);
  meffects = colMeans(ceffects);
}
