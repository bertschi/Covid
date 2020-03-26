## Look at financial markets via options

library(tidyverse)
library(ggthemes)
library(quantmod)
library(lubridate)
library(tidybayes)
library(rstan)
library(loo)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

sp500_opts <- getOptionChain("^GSPC", NULL)
sp500 <- getQuote("^GSPC") %>%
    as_tibble()

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

spot <- sp500 %>%
    mutate(Date = date(`Trade Time`),
           Spot = Last) %>%
    select(Date, Spot)

delta <- 0.025
data <- df_opts %>%
    mutate(Date = spot$Date,
           Spot = spot$Spot,
           Price = (Ask + Bid) / 2,
           Maturity = Expiration - Date) %>%
    ## select out of the money options that are traded
    filter(Bid > 0) %>%
    filter(ifelse(Type == "Call",
                  (1 - delta) * Spot <= Strike,
                  (1 + delta) * Spot >= Strike))

r <- 0.0

fits <- data %>%
    crossing(K = 1:3) %>%
    group_by(Expiration, Maturity, K, Date, Spot) %>%
    nest() %>%
    mutate(fits = pmap(list(data, Maturity, K),
                      function (data, tau, K)
                          sampling(model,
                                   data = list(N = nrow(data),
                                               prices = data$Price,
                                               strikes = data$Strike,
                                               callput = (data$Type == "Call"),
                                               K = K,
                                               tau = as.numeric(unique(tau)),
                                               r = r),
                                   chains = 3))) %>%
    ungroup()

## Visual model comparison based on LOO
fits %>%
    mutate(log_lik = map(fits,
                         ~ loo(.x)$pointwise %>% as_tibble())) %>%
    unnest(log_lik) %>%
    ggplot(aes(elpd_loo,
               fill = factor(K))) +
    geom_density(alpha = 0.4) +
    facet_grid(Expiration ~ .) +
    scale_fill_colorblind() +
    theme_tufte() +
    theme(text = element_text(size = 16),
          legend.position = "top") +
    labs(x = "LOO log likelihood",
         y = NULL,
         fill = "K",
         title = "S&P 500",
         subtitle = paste("Option chain from", today(),
                          sep = " "))

ggsave("reports/model_loo.pdf")

## Compute spot price from put-call parity
## spots <- data %>%
##     select(Expiration, Strike, Type, Price, Maturity) %>%
##     group_by(Expiration) %>%
##     spread(Type, Price) %>%
##     drop_na() %>%
##     mutate(Spot = Call - Put + exp(- r * as.numeric(Maturity)) * Strike) %>%
##     ungroup()

pdfs <- fits %>%
    mutate(pdfs = pmap(list(fits, Spot, Maturity),
                      function (fit, spot, tau) {
                          fit %>%
                              spread_draws(F[k], sigma[k], theta[k]) %>%
                              group_by(.chain, .iteration, .draw, k) %>%
                              mutate(pdf = pmap(list(theta, F, sigma),
                                                function (theta, F, sigma)
                                                    tibble(S = seq(spot / 3, spot * 2,
                                                                   length.out = 101)) %>%
                                                    mutate(pdf = theta * dlnorm(S, log(F), sqrt(as.numeric(tau)) * sigma)))) %>%
                              ungroup() %>%
                              unnest(pdf)
                      })) %>%
    unnest(pdfs)

pdfs %>%
    filter(K == 2) %>%
    filter(.iteration < 100) %>%
    ## group_by(.chain, .iteration, S) %>%
    ## summarize(pdf = sum(pdf)) %>%
    ggplot(aes(S, pdf,
               ## group = .iteration)) +
               color = factor(k),
               group = interaction(.iteration, .chain, k))) +
    geom_line(alpha = 0.1) +
    geom_vline(aes(xintercept = Spot)) +
    scale_color_colorblind() +
    facet_grid(Expiration ~ .) +
    ## theme_economist_white() +
    theme_tufte() +
    theme(text = element_text(size = 12),
          legend.position = "top") +
    labs(x = "Stock price",
         y = "Implied risk-neutral density",
         color = "Component",
         title = "S&P 500",
         subtitle = paste("Option chain from", today(),
                          sep = " "))

ggsave("reports/risk_neutral.pdf")

fits %>%
    mutate(theta = map(fits,
                       ~ .x %>%
                           spread_draws(theta[k]))) %>%
    unnest(theta) %>%
    filter(K == 2 & k == 1) %>%
    ggplot(aes(factor(Expiration), theta)) +
    geom_tufteboxplot() +
    theme_tufte() +
    theme(text = element_text(size = 12)) +
    labs(x = "Expiration date",
         y = "Crash probability",
         title = "S&P 500",
         subtitle = paste("Option chain from", today(),
                          sep = " ")) +
    coord_cartesian(ylim = c(0, 0.5))

ggsave("reports/crash_prob.pdf")
