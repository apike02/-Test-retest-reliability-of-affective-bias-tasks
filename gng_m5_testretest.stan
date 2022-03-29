// gng_m5 test-retest inspired by https://github.com/Nathaniel-Haines/Reliability_2020/blob/master/Code/Stan/joint_RT_normal.stan

data {
    int<lower=1> N;      // # of subjects
    int<lower=1> T;  // max # of trials across subjects
    int<lower=1> N_time; // number of timepoints
    int<lower=1,upper=T> T_subj[N, N_time];  // # of trials within subjects, timepoints
    int<lower=1, upper=4> cue[N, N_time, T]; //cue displayed
    int<lower=-1, upper=1> pressed[N, N_time, T]; //their response (go or nogo)
    real outcome[N, N_time, T]; //the outcome they received
}

transformed data {
  vector[4] initV;
  initV = rep_vector(0.0, 4);
}

parameters {
  // Group-level correlation matrix (cholesky factor for faster computation)
  cholesky_factor_corr[N_time] L_R_xi; 
  cholesky_factor_corr[N_time] L_R_ep;
  cholesky_factor_corr[N_time] L_R_b; 
  cholesky_factor_corr[N_time] L_R_app;
  cholesky_factor_corr[N_time] L_R_av; 
  cholesky_factor_corr[N_time] L_R_rhoRew;
  cholesky_factor_corr[N_time] L_R_rhoPun; 

  // Group-level parameter means
  vector[N_time] xi_mean;
  vector[N_time] ep_mean;
  vector[N_time] b_mean;
  vector[N_time] app_mean;
  vector[N_time] av_mean;
  vector[N_time] rhoRew_mean;
  vector[N_time] rhoPun_mean;

  // Group-level parameter SDs
  vector<lower=0>[N_time] xi_sd;
  vector<lower=0>[N_time] ep_sd;
  vector<lower=0>[N_time] b_sd;
  vector<lower=0>[N_time] app_sd;
  vector<lower=0>[N_time] av_sd;
  vector<lower=0>[N_time] rhoRew_sd;
  vector<lower=0>[N_time] rhoPun_sd;

  // Individual-level parameters (before being transformed)
  matrix[N_time,N] xi_pr; 
  matrix[N_time,N] ep_pr; 
  matrix[N_time,N] b_pr; 
  matrix[N_time,N] app_pr; 
  matrix[N_time,N] av_pr; 
  matrix[N_time,N] rhoRew_pr; 
  matrix[N_time,N] rhoPun_pr; 
}
transformed parameters {
  // Individual-level parameter off-sets (for non-centered parameterization)
  matrix[N_time,N] xi_tilde;
  matrix[N_time,N] ep_tilde;
  matrix[N_time,N] b_tilde;
  matrix[N_time,N] app_tilde;
  matrix[N_time,N] av_tilde;
  matrix[N_time,N] rhoRew_tilde;
  matrix[N_time,N] rhoPun_tilde;

  // Individual-level parameters 
  matrix[N,N_time] xi;
  matrix[N,N_time] ep;
  matrix[N,N_time] b;
  matrix[N,N_time] app;
  matrix[N,N_time] av;
  matrix[N,N_time] rhoRew;
  matrix[N,N_time] rhoPun;
  
  // Construct individual offsets (for non-centered parameterization)
  xi_tilde     = diag_pre_multiply(xi_sd, L_R_xi)         * xi_pr;
  ep_tilde     = diag_pre_multiply(ep_sd, L_R_ep)         * ep_pr; 
  b_tilde      = diag_pre_multiply(b_sd, L_R_b)           * b_pr;
  app_tilde    = diag_pre_multiply(app_sd, L_R_app)       * app_pr; 
  av_tilde     = diag_pre_multiply(av_sd, L_R_av)         * av_pr;
  rhoRew_tilde = diag_pre_multiply(rhoRew_sd, L_R_rhoRew) * rhoRew_pr; 
  rhoPun_tilde = diag_pre_multiply(rhoPun_sd, L_R_rhoPun) * rhoPun_pr;

  // Compute individual-level parameters from non-centered parameterization
  for (time in 1:N_time){
    for (i in 1:N) {
      xi[i,time]     = Phi_approx(xi_mean[time] + xi_tilde[time,i]);
      ep[i,time]     = Phi_approx(ep_mean[time] + ep_tilde[time,i]);
      b[i,time]      = b_mean[time]             + b_tilde[time,i];
      app[i,time]    = app_mean[time]           + app_tilde[time,i];
      av[i,time]     = av_mean[time]            + av_tilde[time,i];
      rhoRew[i,time] = exp(rhoRew_mean[time]    + rhoRew_tilde[time,i]);
      rhoPun[i,time] = exp(rhoPun_mean[time]    + rhoPun_tilde[time,i]);
    }
  }
}
model {
  // Prior on cholesky factor of correlation matrix
  L_R_xi     ~ lkj_corr_cholesky(1);
  L_R_ep     ~ lkj_corr_cholesky(1); 
  L_R_b      ~ lkj_corr_cholesky(1);
  L_R_app    ~ lkj_corr_cholesky(1); 
  L_R_av     ~ lkj_corr_cholesky(1);
  L_R_rhoRew ~ lkj_corr_cholesky(1); 
  L_R_rhoPun ~ lkj_corr_cholesky(1);

  
  // Priors on group-level means 
  xi_mean      ~ normal(0, 1);
  ep_mean      ~ normal(0, 1);
  b_mean       ~ normal(0, 1);
  app_mean     ~ normal(0, 1);
  av_mean      ~ normal(0, 1);
  rhoRew_mean  ~ normal(0, 1);
  rhoPun_mean  ~ normal(0, 1);

  // Priors on group-level SDs
  xi_sd      ~ normal(0, 1);
  ep_sd      ~ normal(0, 1);
  b_sd       ~ normal(0, 1);
  app_sd     ~ normal(0, 1);
  av_sd      ~ normal(0, 1);
  rhoRew_sd  ~ normal(0, 1);
  rhoPun_sd  ~ normal(0, 1);


  // Priors on individual-level parameters
  to_vector(xi_pr)     ~ normal(0, 1);
  to_vector(ep_pr)     ~ normal(0, 1); 
  to_vector(b_pr)      ~ normal(0, 1);
  to_vector(app_pr)    ~ normal(0, 1); 
  to_vector(av_pr)     ~ normal(0, 1);
  to_vector(rhoRew_pr) ~ normal(0, 1); 
  to_vector(rhoPun_pr) ~ normal(0, 1);

  // For each subject
  for (time in 1:N_time){
    for (i in 1:N) {
      vector[4] wv_g;  // action weight for go
      vector[4] wv_ng; // action weight for nogo
      vector[4] qv_g;  // Q value for go
      vector[4] qv_ng; // Q value for nogo
      vector[4] sv;    // stimulus value
      vector[4] pGo;   // prob of go (press)
  
      wv_g  = initV;
      wv_ng = initV;
      qv_g  = initV;
      qv_ng = initV;
      sv    = initV;
      
      for (t in 1:T_subj[i,time]){
        // Choices at time (time) for participant (i) for trial (t)
        real pavlovian;
        if (sv[cue[i, time, t]]>0){
          pavlovian=app[i, time];
        } else {
          pavlovian=av[i, time];
        }
        wv_g[cue[i, time, t]]  = qv_g[cue[i, time, t]] + b[i, time] + pavlovian * sv[cue[i, time, t]];
        wv_ng[cue[i, time, t]] = qv_ng[cue[i, time, t]];  // qv_ng is always equal to wv_ng (regardless of action)
        pGo[cue[i, time, t]]   = inv_logit(wv_g[cue[i, time, t]] - wv_ng[cue[i, time, t]]);
        {  // noise
          pGo[cue[i, time, t]]   *= (1 - xi[i, time]);
          pGo[cue[i, time, t]]   += xi[i, time]/2;
        }
        pressed[i, time, t] ~ bernoulli(pGo[cue[i, time, t]]);
  
        // after receiving feedback, update sv[t + 1]
        if (outcome[i, time, t] >= 0) {
          sv[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - sv[cue[i, time, t]]);
        } else {
          sv[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - sv[cue[i, time, t]]);
        }
  
        // update action values
        if (pressed[i, time, t]) { // update go value
          if (outcome[i, time, t] >=0) {
            qv_g[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - qv_g[cue[i, time, t]]);
          } else {
            qv_g[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - qv_g[cue[i, time, t]]);
          }
        } else { // update no-go value
          if (outcome[i, time, t] >=0) {
            qv_ng[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - qv_ng[cue[i, time, t]]);
          } else {
            qv_ng[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - qv_ng[cue[i, time, t]]);
          }
        }
      }//end of trial loop
    }//end of participant loop
  }//end of time loop
}
generated quantities { 
  // test-retest correlations
  corr_matrix[N_time] R_xi;
  corr_matrix[N_time] R_ep;
  corr_matrix[N_time] R_b;
  corr_matrix[N_time] R_app;
  corr_matrix[N_time] R_av;
  corr_matrix[N_time] R_rhoRew;
  corr_matrix[N_time] R_rhoPun;

  // posterior predictions and log-likelihood
  int post_pred[N, N_time, T];
  real log_lik[N, N_time, T];
  
  // Reconstruct correlation matrices from cholesky factor
  R_xi     = L_R_xi     * L_R_xi';
  R_ep     = L_R_ep     * L_R_ep';
  R_b      = L_R_b      * L_R_b';
  R_app    = L_R_app    * L_R_app';
  R_av     = L_R_av     * L_R_av';
  R_rhoRew = L_R_rhoRew * L_R_rhoRew';
  R_rhoPun = L_R_rhoPun * L_R_rhoPun';

  // initialize LL and post_pred arrays to -1
  for (i in 1:N) {
    post_pred[i,,] = rep_array(-1, N_time, T);
    log_lik[i,,] = rep_array(-1.0, N_time, T);
  }
  
  { // local section, this saves time and space
    for (time in 1:N_time){
      for (i in 1:N) {
        vector[4] wv_g;  // action weight for go
        vector[4] wv_ng; // action weight for nogo
        vector[4] qv_g;  // Q value for go
        vector[4] qv_ng; // Q value for nogo
        vector[4] sv;    // stimulus value
        vector[4] pGo;   // prob of go (press)
  
        wv_g  = initV;
        wv_ng = initV;
        qv_g  = initV;
        qv_ng = initV;
        sv    = initV;
  
        for (t in 1:T_subj[i,time]) {
          real pavlovian;
          if (sv[cue[i, time, t]]>0){
            pavlovian=app[i, time];
          } else {
            pavlovian=av[i, time];
          }
          wv_g[cue[i, time, t]]  = qv_g[cue[i, time, t]] + b[i, time] + pavlovian * sv[cue[i, time, t]];
          wv_ng[cue[i, time, t]] = qv_ng[cue[i, time, t]];  // qv_ng is always equal to wv_ng (regardless of action)
          pGo[cue[i, time, t]]   = inv_logit(wv_g[cue[i, time, t]] - wv_ng[cue[i, time, t]]);
          {  // noise
            pGo[cue[i, time, t]]   *= (1 - xi[i, time]);
            pGo[cue[i, time, t]]   += xi[i, time]/2;
          }
          log_lik[i, time, t] = bernoulli_lpmf(pressed[i, time, t] | pGo[cue[i, time, t]]);
  
          // generate posterior prediction for current trial
          post_pred[i, time, t] = bernoulli_rng(pGo[cue[i, time, t]]);
  
          // after receiving feedback, update sv[t + 1]
          if (outcome[i, time, t] >= 0) {
            sv[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - sv[cue[i, time, t]]);
          } else {
            sv[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - sv[cue[i, time, t]]);
          }
  
          // update action values
          if (pressed[i, time, t]) { // update go value
            if (outcome[i, time, t] >=0) {
              qv_g[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - qv_g[cue[i, time, t]]);
            } else {
              qv_g[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - qv_g[cue[i, time, t]]);
            }
          } else { // update no-go value
            if (outcome[i, time, t] >=0) {
              qv_ng[cue[i, time, t]] += ep[i, time] * (rhoRew[i, time] * outcome[i, time, t] - qv_ng[cue[i, time, t]]);
            } else {
              qv_ng[cue[i, time, t]] += ep[i, time] * (rhoPun[i, time] * outcome[i, time, t] - qv_ng[cue[i, time, t]]);
            }
          }
        } // end of t loop
      } // end of i loop
    } // end of time loop
  } // end of local section
} 
