// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
// from hBayesDM, modified by Alex Pike on 15/02/2022

data {
  int<lower=1> N;               // Number of subjects
  int<lower=1> N_time;          // Number of timepoints
  int<lower=0> Nu_max;          // Max (across subjects) number of upper boundary responses
  int<lower=0> Nl_max;          // Max (across subjects) number of lower boundary responses
  int<lower=0> Nu[N, N_time];   // Number of upper boundary responses for each subj
  int<lower=0> Nl[N, N_time];   // Number of lower boundary responses for each subj
  real RTu[N, N_time, Nu_max];  // upper boundary response times
  real RTl[N, N_time, Nl_max];  // lower boundary response times
  matrix [N, N_time] minRT;     // minimum RT for each subject of the observed data
  real RTbound;                 // lower bound or RT across all subjects (e.g., 0.1 second)
}

parameters {
  // parameters of the DDM (parameter names in Ratcliffs DDM), from https://github.com/gbiele/stan_wiener_test/blob/master/stan_wiener_test.R
  // also see: https://groups.google.com/forum///!searchin/stan-users/wiener%7Csort:relevance/stan-users/-6wJfA-t2cQ/Q8HS-DXgBgAJ
  // alpha (a): Boundary separation or Speed-accuracy trade-off (high alpha means high accuracy). alpha > 0
  // beta (b): Initial bias Bias for either response (beta > 0.5 means bias towards "upper" response 'A'). 0 < beta < 1
  // delta (v): Drift rate Quality of the stimulus (delta close to 0 means ambiguous stimulus or weak ability). 0 < delta
  // tau (ter): Nondecision time + Motor response time + encoding time (high means slow encoding, execution). 0 < ter (in seconds)
  ///* upper boundary of tau must be smaller than minimum RT
  //to avoid zero likelihood for fast responses.
  //tau can for physiological reasons not be faster than 0.1 s.*/
  
  // Group-level correlation matrix (cholesky factor for faster computation)
  cholesky_factor_corr[N_time] L_R_alpha;
  cholesky_factor_corr[N_time] L_R_beta;
  cholesky_factor_corr[N_time] L_R_delta;
  cholesky_factor_corr[N_time] L_R_tau;

  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[N_time] alpha_mean;
  vector[N_time] beta_mean;
  vector[N_time] delta_mean;
  vector[N_time] tau_mean;
  
  vector<lower=0>[N_time] alpha_sd;
  vector<lower=0>[N_time] beta_sd;
  vector<lower=0>[N_time] delta_sd;
  vector<lower=0>[N_time] tau_sd;

  // Individual-level raw parameters (before being transformed)
  matrix[N_time,N] alpha_pr;
  matrix[N_time,N] beta_pr;
  matrix[N_time,N] delta_pr;
  matrix[N_time,N] tau_pr;
}

transformed parameters {
  //Individual-level parameter off-sets (for non-centred parameterization)
  matrix[N_time,N] alpha_tilde;
  matrix[N_time,N] beta_tilde;
  matrix[N_time,N] delta_tilde;
  matrix[N_time,N] tau_tilde;
  
  //Individual_level parameters
  matrix<lower=0>                         [N,N_time] alpha; // boundary separation
  matrix<lower=0, upper=1>                [N,N_time] beta;  // initial bias
  matrix                                  [N,N_time] delta; // drift rate
  matrix<lower=RTbound, upper=max(minRT)> [N,N_time] tau; // nondecision time
  
  //Construct individual offsets (for non-centred parameterization)
  alpha_tilde = diag_pre_multiply(alpha_sd, L_R_alpha)   * alpha_pr;
  beta_tilde  = diag_pre_multiply(beta_sd, L_R_beta)     * beta_pr;
  delta_tilde = diag_pre_multiply(delta_sd, L_R_delta)   * delta_pr;
  tau_tilde   = diag_pre_multiply(tau_sd, L_R_tau)       * tau_pr;

  for (time in 1:N_time){
    for (i in 1:N) {
      alpha[i, time] =        exp(alpha_mean[time] + alpha_tilde[time,i]);
      beta [i, time] = Phi_approx(beta_mean[time]  + beta_tilde[time,i]);
      delta[i, time] =            delta_mean[time] + delta_tilde[time,i];
      tau  [i, time] = Phi_approx(tau_mean[time]   + tau_tilde[time,i]) * (minRT[i, time] - RTbound) + RTbound;
    } // end of subj loop
  }// end of time loop
}

model {
  //Prior on cholesky factors of correlation matrix
  L_R_alpha    ~ lkj_corr_cholesky(1);
  L_R_beta     ~ lkj_corr_cholesky(1);
  L_R_delta    ~ lkj_corr_cholesky(1);
  L_R_tau      ~ lkj_corr_cholesky(1);
  
  // Priors on group-level means
  alpha_mean ~ normal(0, 1);
  beta_mean  ~ normal(0, 1);
  delta_mean ~ normal(0, 1);
  tau_mean   ~ normal(0, 1);
  
  //Priors on group level SDs
  alpha_sd ~ normal(0, 0.2);
  beta_sd  ~ normal(0, 0.2);
  delta_sd ~ normal(0, 0.2);
  tau_sd   ~ normal(0, 0.2);

  // Individual parameters for non-centered parameterization
  to_vector(alpha_pr) ~ normal(0, 1);
  to_vector(beta_pr)  ~ normal(0, 1);
  to_vector(delta_pr) ~ normal(0, 1);
  to_vector(tau_pr)   ~ normal(0, 1);

  // Begin time loop
  for (time in 1:N_time){
  // Begin subject loop
    for (i in 1:N) {
      // Response time distributed along wiener first passage time distribution
      RTu[i, time, :Nu[i, time]] ~ wiener(alpha[i, time], tau[i, time], beta[i, time], delta[i, time]);
      RTl[i, time, :Nl[i, time]] ~ wiener(alpha[i, time], tau[i, time], 1-beta[i, time], -delta[i, time]);
  
    } // end of subject loop
  } // end of time loop
}

generated quantities {
  // test-retest correlations
  corr_matrix[N_time] R_alpha;
  corr_matrix[N_time] R_beta;
  corr_matrix[N_time] R_delta;
  corr_matrix[N_time] R_tau;
  
  // For log likelihood calculation
  real log_lik[N, N_time];
  
  // Reconstruct correlation matrices from cholesky factor
  R_alpha  = L_R_alpha  * L_R_alpha';
  R_beta   = L_R_beta   * L_R_beta';
  R_delta  = L_R_delta  * L_R_delta';
  R_tau    = L_R_tau    * L_R_tau';

  { // local section, this saves time and space
    // Begin time loop
    for (time in 1:N_time){
      // Begin subject loop
      for (i in 1:N) {
        log_lik[i, time] = wiener_lpdf(RTu[i, time, :Nu[i, time]] | alpha[i, time], tau[i, time], beta[i, time], delta[i, time]);
        log_lik[i, time] += wiener_lpdf(RTl[i,time, :Nl[i, time]] | alpha[i, time], tau[i, time], 1-beta[i, time], -delta[i, time]);
      } // end subject loop
    } // end time loop
  } // end local section
} 

