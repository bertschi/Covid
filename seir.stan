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
  int<lower=0> T;
  int<lower=0>  infected[T];
  real<lower=0> times[T]; // observation times
  int<lower=1> N; // population size
}
transformed data {
  real x_r[0];
  int  x_i[1];
  
  x_i[1] = N;
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
  x = integrate_ode_rk45(seir, x0, -0.1, times, theta, x_r, x_i);
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

  for (t in 1:T)
    log_lik[t] = neg_binomial_2_lpmf(infected[t] | x[t, 3], phi);
}
