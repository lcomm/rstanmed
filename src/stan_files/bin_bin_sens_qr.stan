// functions {
// #include /chunks/basic_functions.stan
// #include /chunks/bin_bin_mediation_functions.stan
//   row_vector get_A(matrix S) {
//     int K = cols(S);
//     row_vector[K - 1] S21 = sub_col(S, 1, K, K - 1)';
//     return S21 * inverse(block(S, 1, 1, K - 1, K - 1));
//   }
//   real get_last_cmean_A(vector y, vector mu, row_vector A) {
//     int K = num_elements(mu);
//     real ans;
//     ans = mu[K] + A * (y - head(mu, K - 1));
//     return ans;
//   } 
//   real get_last_cmean_S(vector y, vector mu, matrix S) {
//     int K = cols(S);
//     real ans;
//     row_vector[K - 1] S21 = sub_col(S, 1, K, K - 1)';
//     ans = mu[K] + S21 * inverse(block(S, 1, 1, K - 1, K - 1)) * 
//         (y - head(mu, K - 1));
//     return ans;
//   } 
//   real get_csd(matrix S) {
//     int K = cols(S);
//     vector[K - 1] S12 = sub_col(S, 1, K, K - 1);
//     return sqrt(S[K, K] - S12' * inverse(block(S, 1, 1, K - 1, K - 1)) * S12);
//   }
// }
// data {
//   int<lower=0> N; // # of unique observation patterns
//   int<lower=0> P_m; // # of covariates in mediator model
//   int<lower=0> P_y; // # of covariates in outcome model
//   int<lower=0> P_u; // # of covariates in outcome model
//   int<lower=0,upper=1> y[N]; // collapsed observed outcome
//   int<lower=0,upper=1> m[N]; // collapsed observed mediator
//   matrix[N, P_m - 1] x_m; // collapsed design matrix for mediator (excludes u)
//   matrix[N, P_y - 1] x_y; // collapsed design matrix for outcome (excludes u)
//   matrix[N, P_u] x_u; // collapsed design matrix for u
//   vector[P_m] prior_mean_beta;
//   vector[P_y] prior_mean_alpha;
//   vector[P_u] prior_mean_gamma;
//   cov_matrix[P_m] prior_vcov_beta;
//   cov_matrix[P_y] prior_vcov_alpha;
//   cov_matrix[P_u] prior_vcov_gamma;
//   int u_ei; // whether or not u is exposure-induced
//   int am_intx; // whether there is exposure-mediator interaction
//   int mean_only; // whether to do means of final counterfactual
//   int w[N]; //count pattern for frequency of observations
// } 
// transformed data {
//   // qr stuff
//   // u 
//   matrix[N, P_u] XQu = qr_thin_Q(x_u) * sqrt(N - 1);
//   matrix[P_u, P_u] XRu = qr_thin_R(x_u) / sqrt(N - 1);
//   matrix[P_u, P_u] XRu_inv = inverse(XRu);
//   matrix[P_u, P_u] prior_gt_vcov = XRu * prior_vcov_gamma * XRu';
//   vector[P_u] prior_gt_mean = XRu * prior_mean_gamma;
//   // m
//   matrix[N, P_m - 1] XQm = qr_thin_Q(x_m) * sqrt(N - 1);
//   matrix[P_m - 1, P_m - 1] XRm = qr_thin_R(x_m) / sqrt(N - 1);
//   matrix[P_m, P_m] XRm_ss = ;
//   matrix[P_m - 1, P_m - 1] XRm_inv = inverse(XRm);
//   matrix[P_m - 1, P_m - 1] prior_btnou_vcov = XRm * block(prior_vcov_beta, 1, 1, P_m - 1, P_m - 1) * XRm';
//   vector[P_m - 1] prior_btnou_mean = XRm * head(prior_btnou_mean, P_m - 1);
//   row_vector[P_m  - 1] Am = get_A(XRm * prior_vcov_beta * XRm');
//   // y 
//   matrix[N, P_y - 1] XQy = qr_thin_Q(x_y) * sqrt(N - 1);
//   matrix[P_y - 1, P_y - 1] XRy = qr_thin_R(x_y) / sqrt(N - 1);
//   matrix[P_y - 1, P_y - 1] XRy_inv = inverse(XRy);
//   matrix[P_y - 1, P_y - 1] prior_atnou_vcov = XRy * block(prior_vcov_alpha, 1, 1, P_y - 1, P_y - 1) * XRy';
//   vector[P_y - 1] prior_atnou_mean = XRy * head(prior_atnou_mean, P_y - 1);
//   row_vector[P_y  - 1] Ay = get_A(XRy * prior_vcov_alpha * XRy');
//   //
//   //matrix[P_m, P_m] beta_L;
//   //matrix[P_y, P_y] alpha_L;
//   //matrix[P_u, P_u] gamma_L;
//   int N_tot = sum(w);
//   vector[N] w_vec = to_vector(w);
//   // uncollapsed (i.e., non-count) versions of design matrices
//   // used only in generated quantities block
//   matrix[N_tot, cols(x_y)] big_x_y = uncollapse_matrix(x_y, w);
//   matrix[N_tot, cols(x_m)] big_x_m = uncollapse_matrix(x_m, w);
//   matrix[N_tot, cols(x_u)] big_x_u = uncollapse_matrix(x_u, w);
//   //beta_L  = cholesky_decompose(prior_vcov_beta);
//   //alpha_L = cholesky_decompose(prior_vcov_alpha);
//   //gamma_L = cholesky_decompose(prior_vcov_gamma);
// }
// parameters {
//   vector[P_u] gamma_tilde;
//   vector[P_m - 1] beta_tilde_nou;
//   real beta_u;
//   vector[P_y - 1] alpha_tilde_nou;
//   real alpha_u;
// }
// transformed parameters {
//   //vector[P_m] beta;
//   //vector[P_y] alpha;
//   
//   //vector[P_y - 1] alpha_no_u = head(alpha, P_y - 1);
//   //real alpha_u = alpha[P_y];
//   //vector[P_m - 1] beta_no_u = head(beta, P_m - 1);
//   //real beta_u = beta[P_m];
// }
// model {
//   vector[N] ll_marg;
//   real pu1; 
//   real ll_if_u1;
//   real ll_if_u0;
//   
//   // priors
//   target += multi_normal_lpdf(beta_tilde_nou | prior_btnou_mean, prior_btnou_vcov);
//   target += normal_lpdf(beta_u | , );
//   target += multi_normal_lpdf(alpha_tilde_nou | prior_atnou_mean, prior_atnou_vcov);
//   target += normal_lpdf(alpha_u | , );
//   
//   // likelihood
//   for (n in 1:N) {
//     pu1 = inv_logit(XQu[n, ] * gamma_tilde);
//     ll_if_u0 = bernoulli_logit_lpmf(y[n] | XQy[n, ] * alpha_tilde_nou) + 
//                bernoulli_logit_lpmf(m[n] | XQm[n, ] * beta_tilde_nou);
//     ll_if_u1 = bernoulli_logit_lpmf(y[n] | XQy[n, ] * alpha_tilde_nou + alpha_u) + 
//                bernoulli_logit_lpmf(m[n] | XQm[n, ] * beta_tilde_nou + beta_u);
//     ll_marg[n] = log_mix(pu1, ll_if_u1, ll_if_u0);
//   }
//   target += dot_product(w_vec, ll_marg);
//   
// }
// generated quantities {
// 
//   vector[P_y] alpha = [ XRy_inv * alpha_tilde_nou, alpha_u ];
//   vector[P_m] beta = [ XRm_inv * beta_tilde_nou, beta_u ];
//   vector[P_u] gamma = XRu_inv * gamma_tilde;
//   
//   vector[3] meffects;
//   meffects = bbootColMeans_rng(sim_ceffects_rng(alpha, beta, gamma, u_ei, 
//                                                 big_x_y, big_x_m, big_x_u, 
//                                                 am_intx, mean_only));
// }
