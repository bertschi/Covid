// Stochastic model with delayed dynamics
functions {
  // Weighting function for infection window
  real win(real tau, real mu, real sigma) {
    return exp(normal_lpdf(tau | mu, sigma));
  }
  // Reparameterized lognormal
  real log_normal_z(real n, real p, real delta, real z) {
    // Match mean and variance of overdispersed binomial
    real muN  = n * p;
    real varN = n * p * (1 - p) * delta;
    // Compute corresponding log normal parameters (see ProbOnto) and transform standard normal variate z
    return exp(log(muN / sqrt(1 + varN / muN^2)) + sqrt(log(1 + varN / muN^2)) * z);
  }

  real beta_change(real t, real beta_1, real beta_2, real t_change, real temp_change) {
    return beta_1 + (beta_2 - beta_1) * inv_logit((t - t_change) / temp_change);
  }
}
data {
  int<lower=0> T;           // number/days of observations
  int<lower=0> infected[T]; // daily new infection counts
  int<lower=0> dead[T];     // daily new casualties
  int<lower=1> N;           // population size
  int<lower=0> T_pred;      // number of desired predictions
}
transformed data {
  int tau_max = 21;
  real pobs = 0.5; // Fraction of observed infections
}
parameters {
  // Spreading parameters
  real<lower=0> beta_1;
  real<lower=0> beta_2;
  real<lower=0, upper=T> t_change;
  real<lower=0> temp_change;
  // Window parameters
  real<lower=0, upper=tau_max> mu_win;
  real<lower=0> sigma_win;
  // Observation parameters
  real<lower=0> phi_inf;
  real<lower=0, upper=1> pdie;
  real<lower=0> phi_die;
  // Latent infections
  real<lower=0> delta;
  real<lower=0> inf_0;
  real inf_raw[T];
}
transformed parameters {
  // Model dynamics
  real inf[T + 1];

  inf[1] = inf_0;
  for (t in 1:T) {
    real S = N - sum(inf[1:t]);

    real I = 0;
    for (tau in 0:tau_max) {
      if (t - tau > 0)
	I += win(tau, mu_win, sigma_win) * inf[t - tau];
      else
	I += win(tau, mu_win, sigma_win) * inf[1];
    }
    // Compute next state infection probability
    {
      real beta = beta_change(t, beta_1, beta_2, t_change, temp_change);
      real p = 1 - exp(- beta * I / N);
      inf[t + 1] = log_normal_z(S, p, delta, inf_raw[t]);
    }
  }
  // print(inf);
}
model {
  // Priors
  beta_1 ~ normal(0, 5);
  beta_2 ~ normal(0, 5);
  temp_change ~ student_t(3, 0, 10); // t_change implicitly uniform

  sigma_win ~ normal(0, tau_max / 2); // mu_win implictly uniform
  phi_inf ~ normal(0, 10);
  phi_die ~ normal(0, 10); // pdie implicitly uniform
  
  delta ~ normal(0, 10);
  inf_0 ~ student_t(3, 0, 10);
  inf_raw ~ normal(0, 1); // transformed to log normal above
  // Likelihood
  for (t in 1:T) {
    infected[t] ~ neg_binomial_2(pobs * inf[t + 1], phi_inf);
    dead[t] ~ neg_binomial_2(pdie * inf[t + 1], phi_die);
  }
}
generated quantities {
  real log_lik[T];
  real inf_pred[T + T_pred + 1];
  int infected_pred[T + T_pred];
  int dead_pred[T + T_pred];
  
  for (t in 1:T)
    log_lik[t] = neg_binomial_2_lpmf(infected[t] | pobs * inf[t + 1], phi_inf)
      + neg_binomial_2_lpmf(dead[t] | pdie * inf[t + 1], phi_inf);

  // Generate predictions
  inf_pred[1] = inf[1];
  for (t in 1:(T + T_pred)) {
    real S = N - sum(inf_pred[1:t]);

    real I = 0;
    for (tau in 0:tau_max) {
      if (t - tau > 0)
	I += win(tau, mu_win, sigma_win) * inf_pred[t - tau];
      else
	I += win(tau, mu_win, sigma_win) * inf_pred[1];
    }
    // Compute next state infection probability
    {
      real beta = beta_change(t, beta_1, beta_2, t_change, temp_change);
      real p = 1 - exp(- beta * I / N);
      real tmp_z;
      if (t > T)
	tmp_z = normal_rng(0, 1);
      else
	tmp_z = inf_raw[t];
      
      inf_pred[t + 1] = log_normal_z(S, p, delta, tmp_z);
    }
  }

  for (t in 1:(T + T_pred)) {
    infected_pred[t] = neg_binomial_2_rng(pobs * inf_pred[t + 1], phi_inf);
    dead_pred[t] = neg_binomial_2_rng(pdie * inf_pred[t + 1], phi_die);
  }
}
