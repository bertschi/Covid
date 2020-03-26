// Simple hierarchical growth model
functions {
  real growth(real t, vector theta) {
    // Simple logistic function
    return theta[1] * 1 / (1 + exp(- (t - theta[2]) / theta[3])) + theta[4];
  }
}
data {
  int<lower=0> T; // Number of observations
  int<lower=1> C; // Number of countries
  int<lower=0> cases[C, T];
  int<lower=0> deaths[C, T];
}
parameters {
  vector<lower=0>[4] theta[C];     // Country specific growth parameters
  real<lower=0, upper=1> p_obs[C]; // Country specific observation prevalence
  real<lower=0> tau_obs[C];        // Delay until infection can be identified
  real<lower=0, upper=1> p_death;  // Global death probability
  real<lower=0> tau_die[C];        // Additional delay until death is reported
  real<lower=0> phi[C, 2];         // Observation noise overdispersion
}
model {
  // Priors
  for (c in 1:C) {
    theta[c] ~ student_t(3, 0, 1);
    // p_obs uniform
    tau_obs[c] ~ normal(0, 10);
    {
      real mu = 0.01;
      real nu = 10;
      p_death ~ beta(mu * nu, (1 - mu) * nu);
    }
    tau_die[c] ~ normal(0, 10);
    phi[c] ~ normal(0, 10);
  }
  // Likelihood
  for (c in 1:C) {
    for (t in 2:T) {
      cases[c, t] - cases[c, t-1] ~ neg_binomial_2(p_obs[c] * (growth(t - tau_obs[c], theta[c]) - growth(t - tau_obs[c] - 1, theta[c]) + 1e-8), phi[c, 1]);
      deaths[c, t] - deaths[c, t-1] ~ neg_binomial_2(p_death * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
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
					 | p_death * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
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
    deaths_pred[c, 1] = neg_binomial_2_rng(p_death * (growth(1 - tau_obs[c] - tau_die[c], theta[c]) + 1e-8), phi[c, 2]);
    for (t in 2:T) {
      cases_pred[c, t] = cases_pred[c, t-1] + neg_binomial_2_rng(p_obs[c] * (growth(t - tau_obs[c], theta[c]) - growth(t - tau_obs[c] - 1, theta[c]) + 1e-8), phi[c, 1]);
      deaths_pred[c, t] = deaths_pred[c, t-1] + neg_binomial_2_rng(p_death * (growth(t - tau_obs[c] - tau_die[c], theta[c]) - growth(t - tau_obs[c] - tau_die[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
}
