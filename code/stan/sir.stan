// Stan model assuming SIR dynamics
functions {
  real[] sir(real t, real[] x, real[] theta, real[] x_r, int[] x_i) {
    // x = (S, I, R)
    real dxdt[3];

    // Parameters of logistic beta dynamics
    real level_1 = theta[1];
    real level_2 = theta[2];
    real tau     = theta[3];
    real temp    = theta[4];
    // SIR parameters
    real gamma = theta[5];
    real beta = level_1 + (level_2 - level_1) * inv_logit((t - tau) / temp);
    
    int N = x_i[1];

    dxdt[1] = - beta * x[3] / N * x[1];
    dxdt[2] = beta * x[3] / N * x[1] - gamma * x[2];
    dxdt[3] = gamma * x[2];
    return dxdt;
  }
}
data {
  int<lower=0> T;         // number/days of observations
  int<lower=0> cases[T];      // total case counts
  int<lower=0> deaths[T];     // total death counts
  int<lower=1> N;             // population size
  int<lower=0> T_pred;        // number of desired predictions
  real<lower=0, upper=1> cfr; // case fatality rate
}
transformed data {
  // Parameters passed to ODE solver
  real x_r[0];
  int  x_i[1];
  // observation times
  real times[T];
  
  x_i[1] = N;
  for (t in 1:T)
    times[t] = t;
}
parameters {
  // SIR parameters
  real<lower=0> level_1;
  real<lower=0> level_2;
  real tau;
  real<lower=0> temp;
  real<lower=0> gamma;
  // Fraction of observed cases
  real<lower=0, upper=1> p_obs;
  // Initial condition
  real<lower=0> xini[3];
  // Observation parameters;
  real<lower=0> phiC;
  real<lower=0> phiD;
}
transformed parameters {
  real theta[5];
  real x0[3] = xini;
  real x[T, 3];

  theta[1] = level_1; theta[2] = level_2;
  theta[3] = tau; theta[4] = temp;
  theta[5] = gamma;
  // Integrate ODE system
  x0[1] = N - xini[1];
  x = integrate_ode_rk45(sir, x0, 0.0, times, theta, x_r, x_i);
}
model {
  // Priors
  level_1 ~ normal(0, 5);
  level_2 ~ normal(0, 5);
  tau ~ normal(0, 10);
  temp ~ normal(0, 10);
  gamma ~ normal(0, 5);

  p_obs ~ uniform(0, 1);
  xini ~ student_t(3, 0, 5);

  phiC ~ normal(0, 10);
  phiD ~ normal(0, 10);

  // Likelihood
  for (t in 2:T) {
    real dx[3] = sir(t - 1, x[t - 1, :], theta, x_r, x_i);
    cases[t] - cases[t - 1] ~ neg_binomial_2(p_obs * (dx[2] + dx[3]), phiC);
    deaths[t] - deaths[t - 1] ~ neg_binomial_2(cfr * dx[3], phiD);
  }
}
generated quantities {
  real log_lik[T - 1];
  real times_pred[T + T_pred];
  real x_pred[T + T_pred, 3];
  int cases_pred[T + T_pred];
  int deaths_pred[T + T_pred];

  for (t in 2:T) {
    real dx[3] = sir(t - 1, x[t - 1, :], theta, x_r, x_i);
    log_lik[t - 1]
      = neg_binomial_2_lpmf(cases[t] - cases[t - 1] | p_obs * (dx[2] + dx[3]), phiC)
      + neg_binomial_2_lpmf(deaths[t] - deaths[t - 1] | cfr * dx[3], phiD);
  }

  // Generate predictions
  for (t in 1:(T + T_pred))
    times_pred[t] = t;
  x_pred = integrate_ode_rk45(sir, x0, 0.0, times_pred, theta, x_r, x_i);
  cases_pred[1] = cases[1];
  deaths_pred[1] = deaths[1];
  for (t in 2:(T + T_pred)) {
    real dx[3] = sir(t - 1, x_pred[t - 1, :], theta, x_r, x_i);
    cases_pred[t] = cases_pred[t - 1] + neg_binomial_2_rng(p_obs * (dx[2] + dx[3]), phiC);
    deaths_pred[t] = deaths_pred[t - 1] + neg_binomial_2_rng(cfr * dx[3], phiD);
  }
}
