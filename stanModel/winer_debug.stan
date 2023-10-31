data {
	int<lower=1> N;      // # of subjects
	int N_cond; // # of conditions
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
  real alpha_mean_base;
  real alpha_mean_delta;
  real delta_mean_base;
  real delta_mean_delta;
  real tau_mean_base;
  real tau_mean_delta;
  real<lower=0> alpha_sd_base;
  real<lower=0> alpha_sd_delta;
  real<lower=0> delta_sd_base;
  real<lower=0> delta_sd_delta;
  real<lower=0> tau_sd_base;
  real<lower=0> tau_sd_delta;
  matrix[N,N_cond] alpha_pr;
  matrix[N,N_cond] delta_pr;
  matrix[N,N_cond] tau_pr;
}
transformed parameters {
  matrix<lower=0>[N, N_cond] alpha; //boundary separation
  real beta = 0.5; //bias
  matrix<lower=0>[N, N_cond] delta; //drift rate
  matrix[N, N_cond] tau; //non-decision time
  for (i in 1:N) {
    alpha[i,1] = alpha_mean_base + alpha_sd_base * alpha_pr[i,1];
    delta[i,1] = delta_mean_base + delta_sd_base * delta_pr[i,1];
    tau[i,1] = Phi_approx(tau_mean_base + tau_sd_base * tau_pr[i,1])*minRT[i,1];
    alpha[i,2] = alpha_mean_delta + alpha_sd_delta * alpha_pr[i,2];
    delta[i,2] = delta_mean_delta + delta_sd_delta * delta_pr[i,2];
    tau[i,2] = Phi_approx(tau_mean_delta + tau_sd_delta * tau_pr[i,2])*minRT[i,2];
  }
}
model {
  alpha_mean_base ~ normal(0, 1);
  alpha_mean_delta ~ normal(0, 1);
  delta_mean_base ~ normal(0, 1);
  delta_mean_delta ~ normal(0, 1);
  tau_mean_base ~ normal(0, 1);
  tau_mean_delta ~ normal(0, 1);
  alpha_sd_base ~ cauchy(0, 5);
  alpha_sd_delta ~ cauchy(0, 5);
  delta_sd_base ~ cauchy(0, 5);
  delta_sd_delta ~ cauchy(0, 5);
  tau_sd_base ~ cauchy(0, 5);
  tau_sd_delta ~ cauchy(0, 5);
  alpha_pr ~ normal(0, 1);
  delta_pr ~ normal(0, 1);
  tau_pr ~ normal(0, 1);
  for (i in 1:N) {
    RTu[i,1,:Nu[i]] ~ wiener(alpha[i,1],tau[i,1],beta,delta[i,1]);
    RTl[i,1,:Nl[i]] ~ wiener(alpha[i,1],tau[i,1],1-beta,-delta[i,1]);
    RTu[i,2,:Nu[i]] ~ wiener(alpha[i,2],tau[i,2],beta,delta[i,2]);
    RTl[i,2,:Nl[i]] ~ wiener(alpha[i,2],tau[i,2],1-beta,-delta[i,2]);
    }
}

generated quantities {
  real log_lik[N, N_cond];
  for (i in 1:N) {
    log_lik[i,] = rep_array(-1.0, N_cond);
  }
  for (i in 1:N) {
      log_lik[i,1] = wiener_lpdf(RTu[i,1,:Nu[i]] |alpha[i,1],tau[i,1],beta,delta[i,1]);
      log_lik[i,1] += wiener_lpdf(RTl[i,1,:Nl[i]] |alpha[i,1],tau[i,1],1-beta,-delta[i,1]);
      log_lik[i,2] = wiener_lpdf(RTu[i,2,:Nu[i]] |alpha[i,2],tau[i,2],beta,delta[i,2]);
      log_lik[i,2] += wiener_lpdf(RTl[i,2,:Nl[i]] |alpha[i,2],tau[i,2],1-beta,-delta[i,2]);
    }
}

