data {
  int<lower=0> N; // # obs
  int<lower=0> P_m; // # of covariates in mediator model
  int<lower=0> P_y; // # of covariates in outcome model
  int<lower=0> P_u; // # of covariates in outcome model
  int<lower=0,upper=1> y[N]; // observed outcome
  int<lower=0,upper=1> m[N]; // observed mediator
  int<lower=0,upper=1> u[N]; // observed confounder
  matrix[N, P_m] x_m; // design matrix for mediator
  matrix[N, P_y] x_y; // design matrix for outcome
  matrix[N, P_u] x_u; // design matrix for u
} 
parameters {
  vector[P_m] beta;
  vector[P_y] alpha;
  vector[P_u] gamma;
}
transformed parameters {
  real alpha_u = alpha[P_y];
  real beta_u = beta[P_m];
}
model {
  vector[N] lp;
  for (n in 1:N) {
    lp[n] = bernoulli_logit_lpmf(y[n] | x_y[n, ] * alpha) + 
            bernoulli_logit_lpmf(m[n] | x_m[n, ] * beta) +
            bernoulli_logit_lpmf(u[n] | x_u[n, ] * gamma);
    target += lp[n];
  }
}

