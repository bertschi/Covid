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

## Show VIX
getSymbols("^VIX")
vix <- VIX %>%
    as_tibble() %>%
    mutate(Date = index(VIX)) %>%
    gather(OHLC, Value,
           matches("VIX"))

vix %>%
    filter(OHLC == "VIX.Close") %>%
    ggplot(aes(Date, Value)) +
    geom_line() +
    theme_tufte() +
    theme(text = element_text(size = 14),
          panel.grid.major.y = element_line(color = "grey")) +
    labs(x = NULL, y = NULL)

ggsave("../../reports/figs/VIX.pdf")

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

model <- stan_model("../stan/black_scholes_mix.stan")

spot <- sp500 %>%
    mutate(Date = date(`Trade Time`),
           Spot = Last) %>%
    select(Date, Spot)

delta <- 0.025 ## Allow 2.5% in the money
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

data %>%
    ggplot(aes(Strike, Price,
               shape = factor(Expiration),
               color = Type)) +
    geom_point() +
    geom_vline(aes(xintercept = Spot),
               color = "grey") +
    scale_color_colorblind() +
    theme_tufte() +
    theme(text = element_text(size = 12),
          legend.position = "top") +
    labs(x = "Strike", y = "Mid price",
         color = NULL, shape = "Expiration date",
         title = "S&P 500",
         subtitle = paste("Option chain from", today()))

ggsave("../../reports/figs/option_chain.pdf")

r <- 0.0

fits <- data %>%
    crossing(K = 1:2) %>%
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
                                   chains = 2))) %>%
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

ggsave("../../reports/figs/model_loo.pdf")

## Compute spot price from put-call parity
## spots <- data %>%
##     select(Expiration, Strike, Type, Price, Maturity) %>%
##     group_by(Expiration) %>%
##     spread(Type, Price) %>%
##     drop_na() %>%
##     mutate(Spot = Call - Put + exp(- r * as.numeric(Maturity)) * Strike) %>%
##     ungroup()

## Compute pricing error of each fit
theo_price <- function(data, fit) {
    data %>%
        mutate(i = 1:n()) %>%
        left_join(fit %>%
                  spread_draws(price_hat[i]) %>%
                  rename(TheoPrice = price_hat))
}

theo_df <- fits %>%
    mutate(pred = map2(data, fits, theo_price)) %>%
    unnest(pred)

theo_df %>%
    gather(LowHigh, Val,
           Bid, Ask) %>%
    ggplot(aes(Strike, TheoPrice,
               color = Type,
               group = interaction(Strike, Type))) +
    geom_boxplot(outlier.size = 0) +
    geom_line(aes(y = Val,
                  group = interaction(Type, LowHigh))) +
    facet_grid(K ~ Expiration) +
    theme_tufte() +
    theme(text = element_text(size = 12),
          legend.position = "top") +
    labs(x = "Strike",
         y = "Price",
         color = NULL)

ggsave("../../reports/figs/price_pred.pdf")

## Show implied risk neutral densities
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

ggsave("../../reports/figs/risk_neutral.pdf")

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

ggsave("../../reports/figs/crash_prob.pdf")
