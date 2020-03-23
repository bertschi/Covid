// Stan model assuming SEIR dynamics
functions {
  real[] seir(real t, real[] x, real[] theta, real[] x_r, int[] x_i) {
    // x = (S, E, I, R)
    real dxdt[4];

    real beta = theta[1]; real a = theta[2]; real gamma = theta[3];
    int N = x_i[1];
    dxdt[1] = - beta * x[3] / N * x[1];
    dxdt[2] = beta * x[3] / N * x[1] - a * x[2];
    dxdt[3] = a * x[2] - gamma * x[3];
    dxdt[4] = gamma * x[3];
    return dxdt;
  }
}
data {
  int<lower=0> T;           // number/days of observations
  int<lower=0> infected[T]; // daily new infection counts
  int<lower=1> N;           // population size
  int<lower=0> T_pred;      // number of desired predictions
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
  // SEIR parameters
  real<lower=0> beta;
  real<lower=0> a;
  real<lower=0> gamma;
  // Initial condition
  real<lower=0> xini[4];
  // Observation parameters;
  real<lower=0> phi;
}
transformed parameters {
  real theta[3];
  real x0[4] = xini;
  real x[T,4];

  theta[1] = beta; theta[2] = a; theta[3] = gamma;
  // Integrate ODE system
  x0[1] = N - xini[1];
  x = integrate_ode_rk45(seir, x0, 0.0, times, theta, x_r, x_i);
}
model {
  // Priors
  xini ~ student_t(3, 0, 5);
  beta ~ normal(0, 5);
  a ~ normal(0, 5);
  gamma ~ normal(0, 5);
  phi ~ normal(0, 10);
  // Likelihood
  infected ~ neg_binomial_2(x[:, 3], phi);
}
generated quantities {
  real log_lik[T];
  real times_pred[T + T_pred];
  real x_pred[T + T_pred, 4];
  int infected_pred[T + T_pred];

  for (t in 1:T)
    log_lik[t] = neg_binomial_2_lpmf(infected[t] | x[t, 3], phi);

  // Generate predictions
  for (t in 1:(T + T_pred))
    times_pred[t] = t;
  
  x_pred = integrate_ode_rk45(seir, x0, 0.0, times_pred, theta, x_r, x_i);
  for (t in 1:(T + T_pred))
    infected_pred[t] = neg_binomial_2_rng(x[t, 3], phi);
}
