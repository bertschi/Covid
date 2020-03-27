// Simple hierarchical growth model
functions {
  real growth(real t, vector theta) {
    // Simple logistic function
    return theta[1] * 1 / (1 + exp(- (t - theta[2]) / theta[3]));
  }
}
data {
  int<lower=0> T; // Number of observations
  int<lower=1> C; // Number of countries
  int<lower=0> cases[C, T];
  int<lower=0> deaths[C, T];
}
parameters {
  vector<lower=0>[3] theta[C];
  real<lower=0> phi[C, 2];
  
  // Hierarchical specification on all other parameters
  real<lower=0> alpha_obs;
  real<lower=0> beta_obs;
  real<lower=0, upper=1> p_obs[C];

  real tau_obs_mu;
  real<lower=0> tau_obs_sigma;
  real tau_obs_raw[C];

  real<lower=0> alpha_death;
  real<lower=0> beta_death;  
  real<lower=0, upper=1> p_death[C];

  real tau_die_mu;
  real<lower=0> tau_die_sigma;  
  real tau_die_raw[C];
}
transformed parameters {
  real tau_obs[C];
  real tau_die[C];

  for (c in 1:C) {
    tau_obs[c] = tau_obs_mu + tau_obs_sigma * tau_obs_raw[c];
    tau_die[c] = tau_die_mu + tau_die_sigma * tau_die_raw[c];
  }
}
model {
  // Priors
  alpha_obs ~ normal(0, 10);
  beta_obs ~ normal(0, 10);
  p_obs ~ beta(alpha_obs, beta_obs);

  alpha_death ~ normal(0, 10);
  beta_death ~ normal(0, 10);
  p_death ~ beta(alpha_death, beta_death);

  tau_obs_mu ~ normal(10, 5);
  tau_obs_sigma ~ student_t(3, 0, 1);
  tau_obs_raw ~ normal(0, 1);

  tau_die_mu ~ normal(10, 5);
  tau_die_sigma ~ student_t(3, 0, 1);
  tau_die_raw ~ normal(0, 1);
  
  for (c in 1:C) {
    theta[c] ~ student_t(3, 0, 1);
    phi[c] ~ normal(0, 10);
  }

  // Likelihood
  for (c in 1:C) {
    for (t in 2:T) {
      cases[c, t] - cases[c, t-1] ~ neg_binomial_2(p_obs[c] * (growth(t - tau_obs[c], theta[c]) - growth(t - tau_obs[c] - 1, theta[c]) + 1e-8), phi[c, 1]);
      deaths[c, t] - deaths[c, t-1] ~ neg_binomial_2(p_death[c] * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
}
generated quantities {
  real log_lik[C, T-1];
  int cases_pred[C, T];
  int deaths_pred[C, T];
  real hidden_pred[C, T];
  
  for (c in 1:C) {
    for (t in 2:T) {
      log_lik[c, t-1] = neg_binomial_2_lpmf(cases[c, t] - cases[c, t-1]
					 | p_obs[c] * (growth(t - tau_obs[c], theta[c]) - growth(t - tau_obs[c] - 1, theta[c]) + 1e-8), phi[c, 1])
                    + neg_binomial_2_lpmf(deaths[c, t] - deaths[c, t-1]
					 | p_death[c] * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
  // Generate predictions
  for (c in 1:C) {
    for (t in 1:T)
      hidden_pred[c, t] = growth(t, theta[c]);
  }
  
  // Note: First case of predicted counts is an approximation!
  for (c in 1:C) {
    cases_pred[c, 1] = neg_binomial_2_rng(p_obs[c] * (growth(1 - tau_obs[c], theta[c]) + 1e-8), phi[c, 1]);
    deaths_pred[c, 1] = neg_binomial_2_rng(p_death[c] * (growth(1 - tau_obs[c] - tau_die[c], theta[c]) + 1e-8), phi[c, 2]);
    for (t in 2:T) {
      cases_pred[c, t] = cases_pred[c, t-1] + neg_binomial_2_rng(p_obs[c] * (growth(t - tau_obs[c], theta[c]) - growth(t - tau_obs[c] - 1, theta[c]) + 1e-8), phi[c, 1]);
      deaths_pred[c, t] = deaths_pred[c, t-1] + neg_binomial_2_rng(p_death[c] * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
}
