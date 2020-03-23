## Run analysis on some data

library(tidyverse)
library(ggthemes)
library(bayesplot)
library(tidybayes)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Download data ... should be done each day to retrieve new cases
df_cases <- read_csv("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_19-covid-Confirmed.csv&filename=time_series_2019-ncov-Confirmed.csv")
df_deaths <- read_csv("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_19-covid-Deaths.csv&filename=time_series_2019-ncov-Deaths.csv")
df_recovered <- read_csv("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_19-covid-Recovered.csv&filename=time_series_2019-ncov-Recovered.csv")

## Extract data for Germany
data <- df %>%
    filter(`Country/Region` == "Germany") %>%
    arrange(Date) %>%
    mutate(Infected = Cases - lag(Cases)) %>%
    select(Date, Infected) %>%
    drop_na()

## Fit model
model <- stan_model("stoch.stan")
fit <- sampling(model,
                data = list(T = nrow(data),
                            infected = data$Infected,
                            N = 80000000,
                            ## Predict one month ahead
                            T_pred = 30))

## Plot results ...
gather_draws(fit,
             infected_pred[day]) %>%
    ggplot(aes(factor(day), .value)) +
    geom_tufteboxplot() +
    theme_tufte() +
    scale_y_log10() +
    geom_point(aes(factor(day), Infected),
               data = data %>%
                   mutate(day = 1:n()),
               color = "red",
               size = 2)
