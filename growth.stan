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
  vector<lower=0>[3] theta[C];     // Country specific growth parameters
  real<lower=0, upper=1> p_obs[C]; // Country specific observation prevalence
  real<lower=0, upper=1> p_death;  // Global death probability
  real<lower=0> tau_death[C];      // Delay in reporting deaths
  real<lower=0> phi[C, 2];         // Observation noise overdispersion
}
model {
  // Priors
  for (c in 1:C) {
    theta[c] ~ student_t(3, 0, 1);
    // p_obs uniform
    {
      real mu = 0.01;
      real nu = 10;
      p_death ~ beta(mu * nu, (1 - mu) * nu);
    }
    tau_death[c] ~ normal(0, 10);
    phi[c] ~ normal(0, 10);
  }
  // Likelihood
  for (c in 1:C) {
    for (t in 2:T) {
      cases[c, t] - cases[c, t-1] ~ neg_binomial_2(p_obs[c] * (growth(t, theta[c]) - growth(t - 1, theta[c]) + 1e-8), phi[c, 1]);
      deaths[c, t] - deaths[c, t-1] ~ neg_binomial_2(p_death * (growth(t - tau_death[c], theta[c]) - growth(t - tau_death[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
}
generated quantities {
  real log_lik[C, T-1];

  for (c in 1:C) {
    for (t in 2:T) {
      log_lik[c, t-1] = neg_binomial_2_lpmf(cases[c, t] - cases[c, t-1]
					 | p_obs[c] * (growth(t, theta[c]) - growth(t - 1, theta[c]) + 1e-8), phi[c, 1])
                    + neg_binomial_2_lpmf(deaths[c, t] - deaths[c, t-1]
					 | p_death * (growth(t - tau_death[c], theta[c]) - growth(t - tau_death[c] - 1, theta[c]) + 1e-8), phi[c, 2]);
    }
  }
}
