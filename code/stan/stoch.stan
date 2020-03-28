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
  
}
data {
  int<lower=0> T;           // number/days of observations
  int<lower=0> infected[T]; // daily new infection counts
  int<lower=1> N;           // population size
  int<lower=0> T_pred;      // number of desired predictions
}
transformed data {
  int tau_max = 21;
  real pobs = 0.5; // Fraction of observed infections
}
parameters {
  real<lower=0> beta;
  // Window parameters
  real<lower=0, upper=tau_max> mu_win;
  real<lower=0> sigma_win;
  // Observation parameters
  real<lower=0> phi;
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
      real p = 1 - exp(- beta * I / N);
      inf[t + 1] = log_normal_z(S, p, delta, inf_raw[t]);
    }
  }
  // print(inf);
}
model {
  // Priors
  beta ~ normal(0, 5);
  sigma_win ~ normal(0, tau_max / 2); // mu_win implictly uniform
  phi ~ normal(0, 10);
  delta ~ normal(0, 10);
  inf_0 ~ student_t(3, 0, 10);
  inf_raw ~ normal(0, 1); // transformed to log normal above
  // Likelihood
  for (t in 1:T)
    infected[t] ~ neg_binomial_2((1 - exp(- pobs)) * inf[t + 1], phi);
}
generated quantities {
  real log_lik[T];
  real inf_pred[T + T_pred + 1];
  int infected_pred[T + T_pred];
  
  for (t in 1:T)
    log_lik[t] = neg_binomial_2_lpmf(infected[t] | (1 - exp(- pobs)) * inf[t + 1], phi);

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
      real p = 1 - exp(- beta * I / N);
      real tmp_z;
      if (t > T)
	tmp_z = normal_rng(0, 1);
      else
	tmp_z = inf_raw[t];
      
      inf_pred[t + 1] = log_normal_z(S, p, delta, tmp_z);
    }
  }

  for (t in 1:(T + T_pred))
        infected_pred[t] = neg_binomial_2_rng((1 - exp(- pobs)) * inf_pred[t + 1], phi);
}
