data {
	int<lower=1> N;      // # of subjects
	int N_cond; // # of conditions
	//int N_time; // # of test_index N_time=2
	int<lower=0> Nu_max; // Max(across subjects) number of upper boundary(correct) responses
	int<lower=0> Nl_max; // Max(across subjects) number of lower boundary(incorrect) responses
	int<lower=0> Nu[N]; // # of upper boundary(correct) response for each subject
	int<lower=0> Nl[N]; // # of lower boundary(incorrect) response for each subject

	real RTu[N, N_cond, Nu_max]; // Upper boundary(correct) RTs for each subject, condition, and trial
  real RTl[N, N_cond, Nl_max]; // Lower boundary(incorrect) RTs for each subject, condition, and trial
  real minRT[N, N_cond]; // Minimum RTs for each subject, condition
  real RTbound; // lower bound for RT across all subjects
}

parameters {
  // Group-level correlation matrix (cholesky factor for faster computation)
  // cholesky_factor_corr[2] L_R_mu;
  // cholesky_factor_corr[2] L_R_sigma;

  // Group-level parameter means
  real alpha_mean_base;
  real alpha_mean_delta;
  real delta_mean_base;
  real delta_mean_delta;
  real tau_mean_base;
  real tau_mean_delta;

  // Group-level parameter SDs
  real<lower=0> alpha_sd_base;
  real<lower=0> alpha_sd_delta;
  real<lower=0> delta_sd_base;
  real<lower=0> delta_sd_delta;
  real<lower=0> tau_sd_base;
  real<lower=0> tau_sd_delta;

  // Individual-level raw parameters (before being transformed)
  matrix[N,N_cond] alpha_pr;
  matrix[N,N_cond] delta_pr;
  matrix[N,N_cond] tau_pr;

}
transformed parameters {
  // Individual-level parameters
  matrix<lower=0>[N, N_cond] alpha; //boundary separation
  real beta = 0.5; //bias
  matrix<lower=0>[N, N_cond] delta; //drift rate
  matrix[N, N_cond] tau; //non-decision time

  // Construct inidividual offsets (for non-centered parameterization)
  // mu_i_delta_tilde = diag_pre_multiply(mu_sd_delta, L_R_mu) * mu_i_delta_pr;
  // sigma_i_delta_tilde = diag_pre_multiply(sigma_sd_delta, L_R_sigma) * sigma_i_delta_pr;

  // Compute individual-level parameters from non-centered parameterization
  for (i in 1:N) {
    // Congruent at time 1
    alpha[i,1] = alpha_mean_base + alpha_sd_base * alpha_pr[i,1];
    delta[i,1] = delta_mean_base + delta_sd_base * delta_pr[i,1];
    tau[i,1] = Phi_approx(tau_mean_base + tau_sd_base * tau_pr[i,1])*minRT[i,1];
    // Incongruent at time 1
    alpha[i,2] = alpha_mean_delta + alpha_sd_delta * alpha_pr[i,2];
    delta[i,2] = delta_mean_delta + delta_sd_delta * delta_pr[i,2];
    tau[i,2] = Phi_approx(tau_mean_delta + tau_sd_delta * tau_pr[i,2])*minRT[i,2];
  }
}
model {
  // Prior on cholesky factor of correlation matrix
  // L_R_mu    ~ lkj_corr_cholesky(1);
  // L_R_sigma ~ lkj_corr_cholesky(1);

  // Priors on group-level means
  alpha_mean_base ~ normal(0, 1);
  alpha_mean_delta ~ normal(0, 1);
  delta_mean_base ~ normal(0, 1);
  delta_mean_delta ~ normal(0, 1);
  tau_mean_base ~ normal(0, 1);
  tau_mean_delta ~ normal(0, 1);
  // mu_mean_base     ~ normal(0, 1);
  // mu_mean_delta    ~ normal(0, 1);
  // sigma_mean_base  ~ normal(0, 1);
  // sigma_mean_delta ~ normal(0, 1);
  // ndt_mean         ~ normal(0, 1);

  // Priors on group-level SDs
  alpha_sd_base ~ cauchy(0, 5);
  alpha_sd_delta ~ cauchy(0, 5);
  delta_sd_base ~ cauchy(0, 5);
  delta_sd_delta ~ cauchy(0, 5);
  tau_sd_base ~ cauchy(0, 5);
  tau_sd_delta ~ cauchy(0, 5);
  // mu_sd_base     ~ normal(0, 1);
  // mu_sd_delta    ~ normal(0, 1);
  // sigma_sd_base  ~ normal(0, 1);
  // sigma_sd_delta ~ normal(0, 1);
  // ndt_sd         ~ normal(0, 1);

  // Priors on individual-level parameters
  for (i in 1:N) {
    for (j in 1:N_cond) {
      alpha_pr[i, j] ~ normal(0, 1);
      delta_pr[i, j] ~ normal(0, 1);
      tau_pr[i, j] ~ normal(0, 1);
    }
  }


  // For each subject
  for (i in 1:N) {
    // Congruent at time 1
    RTu[i,1,:Nu[i]] ~ wiener(alpha[i,1],tau[i,1],beta,delta[i,1]);
    RTl[i,1,:Nl[i]] ~ wiener(alpha[i,1],tau[i,1],1-beta,-delta[i,1]);
    // Incongruent at time 1
    RTu[i,2,:Nu[i]] ~ wiener(alpha[i,2],tau[i,2],beta,delta[i,2]);
    RTl[i,2,:Nl[i]] ~ wiener(alpha[i,2],tau[i,2],1-beta,-delta[i,2]);
    }

  }

generated quantities {
  // test-retest correlations
  // corr_matrix[2] R_mu;
  // corr_matrix[2] R_sigma;

  // posterior predictions and log-likelihood
  // real post_pred_c1_t1[N, T_max];
  // real post_pred_c1_t2[N, T_max];
  // real post_pred_c2_t1[N, T_max];
  // real post_pred_c2_t2[N, T_max];
  real log_lik[N, N_cond];

	// Reconstruct correlation matrix from cholesky factor
  // R_mu = L_R_mu * L_R_mu';
  // R_sigma = L_R_sigma * L_R_sigma';

  // set LL and post_pred arrays to -1
  // post_pred_c1_t1 = rep_array(-1.0, N, T_max);
  // post_pred_c1_t2 = rep_array(-1.0, N, T_max);
  // post_pred_c2_t1 = rep_array(-1.0, N, T_max);
  // post_pred_c2_t2 = rep_array(-1.0, N, T_max);
  for (i in 1:N) {
    log_lik[i,] = rep_array(-1.0, N_cond);
  }

  // Generate posterior predictions
  for (i in 1:N) {
    // Congruent at time 1
      // post_pred_c1_t1[i,t] = shiftlnorm_rng(ndt_i[i,1],
      //                                       mu_i_base[i,1],
      //                                       exp(sigma_i_base[i,1]));
      log_lik[i,1] = wiener_lpdf(RTu[i,1,:Nu[i]] |alpha[i,1],tau[i,1],beta,delta[i,1]);
      log_lik[i,1] += wiener_lpdf(RTl[i,1,:Nl[i]] |alpha[i,1],tau[i,1],1-beta,-delta[i,1]);

    // Incongruent at time 1
      log_lik[i,2] = wiener_lpdf(RTu[i,2,:Nu[i]] |alpha[i,2],tau[i,2],beta,delta[i,2]);
      log_lik[i,2] += wiener_lpdf(RTl[i,2,:Nl[i]] |alpha[i,2],tau[i,2],1-beta,-delta[i,2]);
    }
}

