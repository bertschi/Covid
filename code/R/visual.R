## Visual exploration of Covid data ... follows notes, but based on
## ECDC data

library(tidyverse)
library(ggthemes)
library(lubridate)

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

df_chg <- df_raw %>%
    gather(series, daily,
           cases, deaths) %>%
    ## Compute relative values
    mutate(absolute = daily,
           relative = absolute / popData2018) %>%
    gather(type, change,
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

df %>%
    ggplot(aes(date, value,
               color = country)) +
    geom_line() +
    scale_color_viridis_d() +
    facet_grid(type ~ series,
               scales = "free") +
    scale_y_log10() +
    theme_tufte() +
    theme(legend.position = "none")

day_from_thresh <- function(date, val, thresh) {
    thresh_date <- date[which(val > thresh)[1]]
    date - thresh_date
}

df_thresh <- tribble(~ type, ~ threshold,
                     "absolute", 10,
                     "absolute", 100,
                     "absolute", 1000,
                     "relative", 1e-7,
                     "relative", 1e-6,
                     "relative", 1e-5) %>%
    full_join(df %>%
              select(date, country, type, series, value)) %>%
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

## Replicate figure from notes.R
df_thresh %>%
    filter(country %in% country_codes) %>%
    ggplot(aes(as.numeric(relday), value,
               color = country)) +
    geom_line(size = 1.2) +
    scale_y_log10() +
    ## scale_color_brewer(palette = "Paired") +
    scale_color_viridis_d() +
    facet_grid(type + threshold ~ series,
               scales = "free") +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12)) +
    labs(x = "Relative days",
         y = "Count/fraction")

plot_nytimes <- function (df) {
    df %>%
        ggplot(aes(as.numeric(relday), value)) +
        geom_line(aes(group = country2),
                  data = df %>%
                      rename(country2 = country),
                  color = "grey") +
        geom_line(size = 1.2) +
        scale_y_log10() +
        facet_wrap(~ country) +
        theme_tufte() +
        theme(text = element_text(size = 12)) +
        labs(x = "Relative days",
             y = "Count/fraction")    
}

## Now select nicest alignment and compute case fatality rates for
## different delays between cases and deaths
df_cfr <- df_thresh %>%
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(- uid) %>%
    crossing(delay = seq(0, 21)) %>%
    group_by(type, threshold, country, delay) %>%
    arrange(relday) %>%
    mutate(cfr = deaths / lag(cases, unique(delay))) %>%
    ungroup()

df_cfr %>%
    filter(country %in% country_codes) %>%
    filter(type == "relative" & threshold == 1e-6) %>%
    filter(relday >= 0) %>%
    filter(delay %in% seq(0, 21, by = 4)) %>%
    ggplot(aes(as.numeric(relday), cfr,
               color = delay, group = delay)) +
    geom_line() +
    geom_smooth(method = "lm",
                se = FALSE) +
    ## facet_grid(country ~ type + threshold,
    ##            scales = "free") +
    facet_wrap(~ country) +
    theme_tufte() +
    scale_y_log10() +
    scale_color_viridis_c()

plot_cfr <- function (df_cfr, country) {
    df_cfr %>%
        filter(country == !!country) %>%
        filter(type == "relative" & threshold == 1e-6) %>%
        filter(relday >= 0) %>%
        filter(delay < 12) %>%
        ggplot(aes(as.numeric(relday), cfr)) +
        geom_smooth(method = "lm",
                    se = FALSE,
                    color = "grey") +
        geom_line() +
        facet_wrap(~ delay) +
        theme_tufte() +
        theme(panel.grid.major.y = element_line(color = "grey")) +
        coord_cartesian(ylim = c(0, 0.25))
}

## This is it ... purely visual identification of crucial unknowns!!!

model <- function (data) {
    fit <- coef(lm(cfr ~ relday,
                   data = data))
    tibble(cfr_est = mean(data$cfr), ## fit["(Intercept)"],
           lin_loss = fit["relday"])
}

## model <- function (data) {
##     tibble(cfr_est = mean(data$cfr),
##            lin_loss = var(data$cfr))
## }

df_est <- df_cfr %>%
    filter(country %in% country_codes) %>%
    ## Only estimate from reliable values
    filter(relday >= 0 & is.finite(cfr)) %>%
    group_by(type, threshold, country, delay) %>%
    ## remove groups with few data points
    filter(n() > 10) %>%
    nest() %>%
    mutate(fit = map(data, model)) %>%
    unnest(fit) %>%
    group_by(type, threshold, country) %>%
    arrange(abs(lin_loss)) %>%
    summarize(delay_est = first(delay),
              cfr_est = first(cfr_est)) %>%
    ungroup() %>%
    ## lower bound on unobserved cases
    group_by(type, threshold) %>%
    mutate(case_fct = cfr_est / min(cfr_est))

df_est %>%
    gather(est, value,
           delay_est, cfr_est,
           case_fct) %>%
    ggplot(aes(country, value,
               color = interaction(type, threshold))) +
    geom_point(size = 3,
               alpha = 0.4) +
    facet_grid(est ~ .,
               scales = "free") +
    scale_color_colorblind() +
    theme_tufte() +
    theme(panel.grid.major.y = element_line(color = "grey"))

## Show case data shifted and scaled according to estimates

df_thresh %>%
    group_by_at(vars(-value)) %>%
    mutate(uid = 1:n()) %>%
    ungroup() %>%
    spread(series, value) %>%
    select(- uid) %>%
    filter(country %in% country_codes) %>%
    ## filter(threshold %in% c(1000, 1e-6, 1e-5)) %>%
    left_join(df_est) %>%
    drop_na() %>%
    group_by(type, threshold, country) %>%
    mutate(act_cases = case_fct * cases) %>%
    ggplot(aes(relday + delay_est, act_cases,
               color = country)) +
    geom_line(size = 1.2) +
    scale_y_log10() +
    scale_color_viridis_d() +
    facet_wrap(~ type + threshold,
               scales = "free") +
    coord_cartesian(xlim = c(-10, 40)) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "grey"))

plot_nytimes(df_thresh %>%
             group_by_at(vars(-value)) %>%
             mutate(uid = 1:n()) %>%
             ungroup() %>%
             spread(series, value) %>%
             select(- uid) %>%
             filter(country %in% country_codes) %>%
             filter(type == "relative" & threshold == 1e-6) %>%
             left_join(df_est) %>%
             drop_na() %>%
             group_by(type, threshold, country) %>%
             mutate(relday = relday + delay_est,
                    value = case_fct * cases)) +
    coord_cartesian(xlim = c(-10, 40),
                    ylim = c(1e-7, 1e-1)) +
    theme(panel.grid.major.y = element_line(color = "grey"),
          panel.grid.minor.y = element_line(color = "lightgrey"))
