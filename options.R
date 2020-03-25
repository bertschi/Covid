## Look at financial markets via options

library(tidyverse)
library(ggthemes)
library(quantmod)
library(rstan)

sp500_opts <- getOptionChain("^GSPC", NULL)

df_opts <-
    bind_rows(
        map_df(names(sp500_opts),
               ~ sp500_opts[[.x]]$call %>%
                   mutate(Expiration = .x,
                          Type = "Call")),
        map_df(names(sp500_opts),
               ~ sp500_opts[[.x]]$put %>%
                   mutate(Expiration = .x,
                          Type = "Put"))) %>%
    as_tibble() %>%
    mutate(Expiration =
               readr::parse_date(Expiration,
                                 format = "%b.%d.%Y",
                                 locale = locale(date_names(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"),
                                                            c("Jan", "Feb", "Mär", "Apr", "Mai", "Jun", "Jul", "Aug", "Sep", "Okt", "Nov", "Dez"),
                                                            c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")))))

df_opts %>%
    gather(Quote, Price,
           Ask, Bid) %>%
    filter(Vol > 100) %>%
    ggplot(aes(Strike, Price,
               group = interaction(Quote, Type, Expiration),
               color = factor(Expiration),
               linetype = Quote)) +
    geom_line() +
    theme_tufte() +
    scale_color_colorblind()

## Stan for fitting mixture of log-normal to option prices

code <- "
functions {
  // Black-Scholes formulas using forward prices
  real dplus(real F, real K, real r, real tau, real sigma) {
    return 1 / (sigma * sqrt(tau)) * (log(F) - log(K) + sigma^2 / 2 * tau);
  }
  real dminus(real F, real K, real r, real tau, real sigma) {
    return dplus(F, K, r, tau, sigma) - sigma * sqrt(tau);
  }
  real call(real F, real K, real r, real tau, real sigma) {
    return exp(- r *tau) * (F * Phi(dplus(F, K, r, tau, sigma))
                          - K * Phi(dminus(F, K, r, tau, sigma)));
  }
  real put(real F, real K, real r, real tau, real sigma) {
    return K * exp(- r * tau) - F + call(F, K, r, tau, sigma);
  }
}
data {
  int<lower=0> N;                   // Number of prices
  real<lower=0> prices[N];
  real<lower=0> strikes[N];
  int<lower=0, upper=1> callput[N]; // Call = 1, Put = 0
  int<lower=1> K;                   // Number of mixture components
  // Basic parameters
  real<lower=0> tau; // Time to maturity
  real r;            // Riskless interest rate
}
parameters {
  // Log-normal parameters
  positive_ordered[K] F;   // Means, i.e. forward prices
  real<lower=0> sigma[K];  // Volatilities
  // Other parameters
  simplex[K] theta;
  real<lower=0> sigma_obs; // Observation noise
  real<lower=0> inv_nu;
}
transformed parameters {
  real price_hat[N]; // Theoretical prices

  for (n in 1:N) {
    price_hat[n] = 0;
    for (k in 1:K) {
      if (callput[n] == 1)
        price_hat[n] += theta[k] * call(F[k], strikes[n], r, tau, sigma[k]);
      else
        price_hat[n] += theta[k] * put(F[k], strikes[n], r, tau, sigma[k]);
    }
  }
}
model {
  sigma ~ normal(0, 1);
  theta ~ dirichlet(rep_vector(1.0, K));
  sigma_obs ~ student_t(3, 0, 5);
  inv_nu ~ normal(0, 1);
  F ~ normal(mean(strikes), max(strikes) - min(strikes)); // Keep it at reasonable values

  prices ~ student_t(1 / inv_nu, price_hat, sigma_obs);
}
generated quantities {
  real log_lik[N];

  for (n in 1:N)
    log_lik[n] = student_t_lpdf(prices[n] | 1 / inv_nu, price_hat[n], sigma_obs);
}
"

model <- stan_model(model_code = code)

data <- df_opts %>%
    mutate(Price = (Ask + Bid) / 2,
           Maturity = Expiration - today()) %>%
    filter(Vol > 100) %>%
    filter(Expiration == parse_date("Dec.17.2021", "%b.%d.%Y"))

fit <- sampling(model,
                data = list(N = nrow(data),
                            prices = data$Price,
                            strikes = data$Strike,
                            callput = (data$Type == "Call"),
                            K = 2,
                            tau = as.numeric(unique(data$Maturity)),
                            r = 0.0),
                chains = 2)
