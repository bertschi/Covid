## Impressive data collapse by rescaling time via growth rate

library(tidyverse)
library(ggthemes)
library(lubridate)
library(latex2exp)

## Copied from visual.R

df_raw <- read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv") %>%
    mutate(date = make_date(year, month, day),
           country = countryterritoryCode)

## Long form data including absolute and relative cumulative totals
df <- df_raw %>%
    gather(series, daily,
           cases, deaths) %>%
    ## Compute cumulative total cases
    group_by(country, series) %>%
    arrange(date) %>%
    mutate(absolute = cumsum(daily)) %>%
    ungroup() %>%
    ## Compute relative values
    mutate(relative = absolute / popData2018) %>%
    gather(type, value,
           absolute, relative)

## country_codes <- c("CHN", "ITA", "DEU", "USA", "FRA", "ESP", "GBR", "KOR", "AUT", "SWE", "JPN", "TWN")
country_codes <- c("USA", "CHN", "ESP", "FRA", "ITA", "DEU", "GBR", "AUT", "AUS", "KOR")

country_codes <- df %>%
    filter(popData2018 > 1e6) %>%
    select(date, country, type, series, value) %>%
    group_by(type, series, country) %>%
    summarize(value = max(value)) %>%
    group_by(type, series) %>%
    arrange(desc(value)) %>%
    slice(1:8) %>%
    ungroup() %>%
    distinct(country) %>%
    .$country

country_codes <- c("KOR", country_codes)

day_from_thresh <- function(date, val, thresh) {
    thresh_date <- date[which(val > thresh)[1]]
    date - thresh_date
}

df_thresh <- tribble(~ type, ~ threshold,
                     ## "absolute", 10,
                     ## "absolute", 100,
                     ## "absolute", 1000,
                     ## "relative", 1e-7,
                     ## "relative", 1e-5,
                     "relative", 1e-6,
                     "relative", 2e-6,
                     "relative", 4e-6,
                     "relative", 8e-6,
                     ## "relative", 2e-6,
                     ## "relative", 4e-6,
                     ## "relative", 8e-6,
                     ## "relative", 16e-6
                     ) %>%
    full_join(df %>%
              select(date, country, type, series, value, popData2018)) %>%
    ## Note: spread needs unique index
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(-uid) %>%
    ## Align at first day exceeding threshold for deaths
    group_by(country, type, threshold) %>%
    arrange(date) %>%
    mutate(relday = day_from_thresh(date, deaths, threshold)) %>%
    gather(series, value,
           cases, deaths) %>%
    ungroup()

df_cfr <- df_thresh %>%
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(- uid) %>%
    crossing(delay = seq(0, 30)) %>%
    group_by(type, threshold, country, delay) %>%
    arrange(relday) %>%
    mutate(cfr = deaths / lag(cases, unique(delay))) %>%
    ungroup()

model <- function (data) {
    fit <- coef(lm(cfr ~ relday,
                   data = data))
    tibble(## cfr_est = last(arrange(data, relday)$cfr),
           cfr_est = median(data$cfr),
           lin_loss = fit["relday"],
           num_data = nrow(data))
}

df_est <- df_cfr %>%
    filter(country %in% country_codes) %>%
    ## Only estimate from reliable values
    filter(relday >= 0 & is.finite(cfr) & delay > 0) %>%
    group_by(type, threshold, country, delay) %>%
    ## remove groups with few data points
    filter(n() >= 10) %>%
    nest() %>%
    mutate(fit = map(data, model)) %>%
    unnest(fit) %>%
    group_by(type, threshold, country) %>%
    arrange(abs(lin_loss)) %>%
    summarize(delay_est = first(delay),
              cfr_est = first(cfr_est),
              num_data = first(num_data)) %>%
    ungroup() %>%
    ## lower bound on unobserved cases
    group_by(type, threshold) %>%
    ## mutate(case_fct = cfr_est / min(cfr_est))
    mutate(case_fct = cfr_est / 0.01)

################################################################################
## New figures and rescaling
################################################################################

tau <- 5

df_thresh %>%
    filter(country %in% country_codes &
           series == "deaths") %>%
    group_by(country, threshold) %>%
    arrange(relday) %>%
    mutate(chg = value - lag(value, tau)) %>%
    drop_na() %>%
    mutate(scaling = cumsum(chg)) %>%
    ungroup() %>%
    filter(relday > -10) %>%
    ggplot(aes(scaling, value,
               group = country,
               color = as.numeric(relday))) +
    geom_path() +
    scale_color_viridis_c(direction = -1) +
    scale_y_log10() +
    scale_x_log10() +
    theme_tufte() +
    facet_wrap(~ threshold)

df_thresh %>%
    filter(country %in% country_codes &
           series == "deaths") %>%
    group_by(country, threshold) %>%
    arrange(relday) %>%
    mutate(chg = value - lag(value, tau)) %>%
    drop_na() %>%
    mutate(scaling = cumsum(chg)) %>%
    ungroup() %>%
    filter(relday > -10) %>%
    ggplot(aes(scaling, value,
               color = country)) +
    geom_path(size = 0.8) +
    scale_color_viridis_d() +
    scale_y_log10() +
    scale_x_log10() +
    facet_wrap(~ threshold) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "grey")) +
    labs(x = "Rescaled time",
         y = "Fraction",
         color = NULL)

## Add scaling by deaths globally

df_scaling <- df_thresh %>%
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(- uid) %>%
    group_by(type, threshold, country) %>%
    arrange(relday) %>%
    mutate(chg = deaths - lag(deaths, tau)) %>%
    drop_na() %>%
    mutate(scaling = cumsum(chg)) %>%
    ungroup() %>%
    gather(series, value,
           deaths, cases)

## and replicate above figure including cases

df_scaling %>%
    filter(relday > -10 & country %in% country_codes) %>%
    ggplot(aes(scaling, value,
               color = country)) +
    geom_path(size = 0.8) +
    scale_color_viridis_d() +
    scale_y_log10() +
    scale_x_log10() +
    theme_tufte() +
    facet_grid(threshold ~ series) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "grey")) +
    labs(x = "Rescaled time",
         y = "Fraction",
         color = NULL)

ggsave("../../reports/figs/ecdc_scaling_all.pdf")

df_scaling %>%
    filter(country %in% country_codes) %>%
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(- uid) %>%
    ## filter(threshold %in% c(1000, 1e-6, 1e-5)) %>%
    left_join(df_est) %>%
    drop_na() %>%
    group_by(type, threshold, country) %>%
    mutate(act_cases = lag(case_fct * cases, unique(delay_est)),
           ## Note: due to delay we don't know scaling for most recent data yet!
           act_day = relday + delay_est) %>%
    ungroup() %>%
    filter(relday > -2) %>%
    ggplot(aes(scaling, act_cases,
               color = country)) +
    geom_line(size = 0.8) +
    scale_y_log10() +
    scale_x_log10() +
    scale_color_viridis_d() +
    facet_wrap(~ type + threshold) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "grey")) +
    labs(x = "Rescaled time",
         y = "Fraction",
         color = NULL) +
    coord_cartesian(ylim = c(1e-4, 1e-1))

ggsave("../../reports/figs/ecdc_scaling_estdelay_cases.pdf")
