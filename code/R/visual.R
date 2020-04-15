## Visual exploration of Covid data ... follows notes, but based on
## ECDC data

library(tidyverse)
library(ggthemes)
library(lubridate)
library(latex2exp)

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

## Show raw data

df %>%
    filter(country %in% country_codes) %>%
    filter(type == "absolute") %>%
    ggplot(aes(date, value,
               color = country)) +
    geom_line(size = 0.8) +
    scale_color_viridis_d() +
    facet_grid(series ~ .) +
    scale_y_log10(breaks  = c(1, 10, 100, 1000, 10000, 100000),
                  labels = function (x) format(x, scientific = FALSE)) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 16)) +
    labs(x = NULL,
         y = "Total count",
         color = NULL)

ggsave("../../reports/figs/ecdc_raw_absolute.pdf")

df %>%
    filter(country %in% country_codes) %>%
    filter(type == "relative") %>%
    ggplot(aes(date, value,
               color = country)) +
    geom_line(size = 0.8) +
    scale_color_viridis_d() +
    facet_grid(series ~ .) +
    ## scale_y_log10(breaks  = c(0.01, 0.1, 1, 10, 100, 1000),
    ##               labels = function (x) format(x, scientific = FALSE)) +
    scale_y_log10(breaks = 10 ^ (-7:-3)) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12)) +
    labs(x = NULL,
         y = "Relative count [per inhabitants]",
         color = NULL)

ggsave("../../reports/figs/ecdc_raw_relative.pdf")

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

## Replicate figure from notes.R
df_thresh %>%
    filter(country %in% country_codes) %>%
    filter(type == "relative") %>%
    ggplot(aes(as.numeric(relday), value,
               color = country)) +
    geom_line(size = 0.8) +
    ## scale_color_brewer(palette = "Paired") +
    scale_color_viridis_d() +
    facet_grid(threshold ~ series) +
    theme_tufte() +
    coord_cartesian(xlim = c(-14, 40)) +
    scale_y_log10(breaks = 10 ^ (-7:-3)) +
    theme(legend.position = "top",
          text = element_text(size = 12)) +
    labs(x = "Relative days",
         y = "Relative count",
         color = NULL)

ggsave("../../reports/figs/ecdc_align_all.pdf")

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
             y = "Relative count")    
}

df_thresh %>%
    filter(threshold == 2e-6 &
           series == "deaths" &
           relday > -10 &
           country %in% country_codes) %>%
    plot_nytimes()

ggsave("../../reports/figs/ecdc_aligned_twopermill_nyt.pdf")

## Now select nicest alignment and compute case fatality rates for
## different delays between cases and deaths
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
        theme(panel.grid.major.y = element_line(color = "grey"),
              text = element_text(size = 16)) +
        labs(x = "Relative days since one death per mill.",
             y = "Estimated CFR")
}

plot_cfr(df_cfr, "ITA") +
    coord_cartesian(ylim = c(0, 0.4)) +
    labs(title = "ITA",
         subtitle = TeX("Estimated $cfr_{\\tau}$ for varying delays $\\tau$"))
ggsave("../../reports/figs/ecdc_cfr_delay_ITA.pdf")

plot_cfr(df_cfr, "DEU") +
    coord_cartesian(ylim = c(0, 0.1)) +
    labs(title = "DEU",
         subtitle = TeX("Estimated $cfr_{\\tau}$ for varying delays $\\tau$"))
ggsave("../../reports/figs/ecdc_cfr_delay_DEU.pdf")

## This is it ... purely visual identification of crucial unknowns!!!

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

df_est %>%
    mutate(frac_obs = 1 / case_fct) %>%
    gather(est, value,
           delay_est, cfr_est,
           frac_obs) %>%
    mutate(est = map_chr(est,
                         ~ list(delay_est = "Delay $\\tau$",
                                cfr_est = "CFR",
                                frac_obs = "Fraction of observed cases")[[.x]])) %>%
    ggplot(aes(country, value,
               shape = factor(threshold))) +
    geom_point(## aes(size = num_data),
               size = 1,
               fill = "black",
               position = position_dodge(width = 0.8)) +
    facet_grid(est ~ .,
               scales = "free",
               labeller = function (labels, multi_line = TRUE) {
                   label_parsed(map(labels, ~ TeX(.x)),
                                multi_line)
               }) +
    scale_shape_manual(values = c(1, 5, 21, 23)) +
    theme_tufte() +
    theme(text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "lightgrey"),
          ## panel.border = element_rect(fill = NA, color = "grey")
          ) +
    labs(x = NULL, y = NULL,
         shape = "Threshold",
         size = "# data")

ggsave("../../reports/figs/ecdc_estimates.pdf")

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
    facet_wrap(~ threshold) +
    coord_cartesian(xlim = c(-10, 40)) +
    theme_tufte() +
    theme(legend.position = "top",
          text = element_text(size = 12),
          panel.grid.major.y = element_line(color = "grey")) +
    labs(x = TeX("Relative days, shifted by estimated delay $\\tau$"),
         y = TeX("Relative cases, including estimated unobserved cases"),
         color = NULL)

ggsave("../../reports/figs/ecdc_aligned_case_est.pdf")

plot_nytimes(df_thresh %>%
             group_by_at(vars(-value)) %>%
             mutate(uid = 1:n()) %>%
             ungroup() %>%
             spread(series, value) %>%
             select(- uid) %>%
             filter(country %in% country_codes) %>%
             filter(type == "relative" & threshold == 2e-6) %>%
             left_join(df_est) %>%
             group_by(type, threshold, country) %>%
             mutate(relday = relday + delay_est,
                    value = case_fct * cases)) +
    coord_cartesian(xlim = c(-10, 40),
                    ylim = c(1e-7, 1e-1)) +
    theme(panel.grid.major.y = element_line(color = "grey"),
          panel.grid.minor.y = element_line(color = "lightgrey")) +
    labs(x = "Relative days since two death per mill.",
         y = "Relative count")

ggsave("../../reports/figs/ecdc_nyt_case_est_twopermill.pdf")

## SIR fitting

library(rstan)
library(bayesplot)
library(tidybayes)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

model <- stan_model("../stan/sir.stan")

foo <- df_thresh %>%
    filter(country == "DEU" &
           threshold == 2e-6 &
           relday > -10) %>%
    arrange(relday) %>%
    mutate(value = value * popData2018) %>%
    spread(series, value)

fit <- sampling(model,
                data = list(T = nrow(foo),
                            cases = round(foo$cases),
                            deaths = round(foo$deaths),
                            N = unique(foo$popData2018),
                            T_pred = 10,
                            cfr = 0.01),
                chains = 2)

df_pred <- gather_draws(fit,
                        cases_pred[t], deaths_pred[t], x_pred[t, d]) %>%
    group_by(.iteration, .chain, t, .variable) %>%
    arrange(d) %>%
    summarize(.value = ifelse(length(d) > 1,
                              .value[2] + .value[3],
                              .value))
df_pred %>%
    filter(.iteration < 50) %>%
    group_by(.iteration, .chain, .variable) %>%
    arrange(t) %>%
    mutate(.change = .value - lag(.value)) %>%
    ungroup() %>%
    mutate(.variable = map_chr(.variable,
                               ~ list(cases_pred = "Cases",
                                      deaths_pred = "Deaths",
                                      x_pred = "Model actual")[[.x]])) %>%
    ggplot(aes(t - 10, .value / unique(foo$popData2018),
               color = .variable,
               group = interaction(.variable, .iteration, .chain))) +
    geom_line(alpha = 0.4) +
    scale_y_log10() +
    scale_color_colorblind() +
    theme_tufte() +
    theme(text = element_text(size = 12)) +
    geom_line(aes(t - 10, value / unique(foo$popData2018),
                  group = series),
              data = foo %>%
                  mutate(t = 1:n()) %>%
                  gather(series, value,
                         cases, deaths) %>%
                  group_by(series) %>%
                  arrange(t) %>%
                  mutate(change = value - lag(value)),
              color = "red",
              size = 0.8) +
    labs(x = "Relative days",
         y = "Relative count",
         color = NULL)

ggsave("../../reports/figs/sir_pred_DEU.pdf")

spread_draws(fit, gamma, p_obs) %>%
    mutate(tau = rexp(n(), gamma)) %>%
    rename(`Reporting delay` = tau,
           `Observed fraction` = p_obs) %>%
    gather(param, est,
           `Reporting delay`, `Observed fraction`) %>%
    ggplot(aes(est)) +
    geom_density() +
    facet_wrap(~ param,
               scales = "free") +
    theme_tufte() +
    theme(text = element_text(size = 12)) +
    labs(x = NULL,
         y = "Posterior density")

ggsave("../../reports/figs/sir_est_DEU.pdf")
