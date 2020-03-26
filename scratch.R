## Run analysis on some data

library(tidyverse)
library(ggthemes)
library(bayesplot)
library(tidybayes)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Download data ... should be done each day to retrieve new cases
df_cases <- read_csv("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_confirmed_global.csv&filename=time_series_covid19_confirmed_global.csv")
df_deaths <- read_csv("https://data.humdata.org/hxlproxy/api/data-preview.csv?url=https%3A%2F%2Fraw.githubusercontent.com%2FCSSEGISandData%2FCOVID-19%2Fmaster%2Fcsse_covid_19_data%2Fcsse_covid_19_time_series%2Ftime_series_covid19_deaths_global.csv&filename=time_series_covid19_deaths_global.csv")
## df_recovered <- read_csv("")

df <- df_cases %>%
    gather(Date, Cases,
           ends_with("/20")) %>%
    full_join(df_deaths %>%
              gather(Date, Deaths,
                     ends_with("/20"))) %>%
    ## full_join(df_recovered %>%
    ##           gather(Date, Recovered,
    ##                  ends_with("/20"))) %>%
    mutate(Date = parse_date(Date, "%m/%d/%y"))

## Plot like in Wikipedia https://de.wikipedia.org/wiki/COVID-19-Pandemie_in_Deutschland#/media/Datei:Coronavirus_deaths.png

df_wiki <- df %>%
    filter(`Country/Region` %in% c("China", "Italy", "Germany", "US", "France", "Spain", "UK", "Iran", "Korea, South", "Austria", "Sweden", "Japan")) %>%
    gather(Series, Value,
           Cases, Deaths) %>%
    group_by(`Country/Region`, Date, Series) %>%
    summarize(Value = sum(Value)) %>%
    ungroup()

day_from_thresh <- function(date, val, thresh) {
    thresh_date <- date[which(val > thresh)[1]]
    date - thresh_date
}

df_wiki %>%
    filter(Series %in% c("Cases", "Deaths")) %>%
    ## Align days where count reaches 10 for the first time
    group_by(`Country/Region`, Series) %>%
    arrange(Date) %>%
    mutate(RelDay =
               as.numeric(day_from_thresh(Date, Value, ifelse(unique(Series) == "Cases", 1000, 10)))) %>%
    ggplot(aes(RelDay, Value,
               color = `Country/Region`)) +
    geom_line() +
    scale_color_colorblind() +
    scale_y_log10() +
    theme_tufte() +
    facet_wrap(~ Series)

## Fetch country population sizes
df_pop <- read_csv("population.csv")

coding <- tribble(~ `Country/Region`, ~ `Country Code`,
                  "China", "CHN",
                  "Italy", "ITA",
                  "Germany", "DEU",
                  "US", "USA",
                  "France", "FRA",
                  "Spain", "ESP",
                  "UK", "GBR",
                  "Iran", "IRN",
                  "Korea, South", "KOR",
                  "Austria", "AUT",
                  "Sweden", "SWE",
                  "Japan", "JPN")

df_wiki %>%
    left_join(coding) %>%
    left_join(df_pop %>%
              filter(Year == max(Year)) %>%
              rename(Population = Value)) %>%
    spread(Series, Value) %>%
    mutate(DorR = Deaths + Recovered) %>%
    gather(Series, Value,
           Cases, Deaths, Recovered, DorR) %>%
    mutate(RelValue = Value / Population) %>%
    ggplot(aes(Date, RelValue,
               color = Series)) +
    geom_line() +
    scale_color_colorblind() +
    scale_y_log10() +
    theme_tufte() +
    facet_wrap(~ `Country Code`)

df_wiki %>%
    ## Join with population size
    left_join(coding) %>%
    left_join(df_pop %>%
              filter(Year == max(Year)) %>%
              rename(Population = Value)) %>%
    ## Compute relative value (per 100 000)
    mutate(RelValue = Value / Population * 1e5) %>%
    gather(Type, Val,
           Value, RelValue) %>%
    ## Adjust relative date to deaths (absolute and relative)
    spread(Series, Val) %>%
    group_by(`Country/Region`, Type) %>%
    arrange(Date) %>%
    mutate(RelDay =
               as.numeric(day_from_thresh(Date, Deaths,
                                          ifelse(unique(Type) == "Value",
                                                 10,
                                                 1e-2)))) %>%
    ungroup() %>%
    gather(Series, Valu,
           Cases, Deaths) %>%
    ## Compute changes
    group_by(`Country/Region`, Series, Type) %>%
    arrange(Date) %>%
    mutate(Change = Valu - lag(Valu)) %>%
    ungroup() %>%
    gather(Tag, Val,
           Valu, Change) %>%
    ggplot(aes(RelDay, Val,
               color = `Country/Region`,
               linetype = Series)) +
    ## geom_line() +
    geom_smooth(se = FALSE) +
    scale_y_log10() +
    scale_color_viridis_d() +
    facet_grid(Type ~ Tag,
               scales = "free") +
    theme_tufte()

## Raw data instead of smoothed curve ... not changes are too noisy in this fashion
df_wiki %>%
    ## Join with population size
    left_join(coding) %>%
    left_join(df_pop %>%
              filter(Year == max(Year)) %>%
              rename(Population = Value)) %>%
    ## Compute relative value (per 100 000)
    mutate(RelValue = Value / Population * 1e5) %>%
    gather(Type, Val,
           Value, RelValue) %>%
    ## Adjust relative date to deaths (absolute and relative)
    spread(Series, Val) %>%
    group_by(`Country/Region`, Type) %>%
    arrange(Date) %>%
    mutate(RelDay =
               as.numeric(day_from_thresh(Date, Deaths,
                                          ifelse(unique(Type) == "Value",
                                                 10,
                                                 1e-2)))) %>%
    ungroup() %>%
    gather(Series, Val,
           Cases, Deaths) %>%
    ggplot(aes(RelDay, Val,
               color = `Country/Region`,
               linetype = Series)) +
    geom_line() +
    scale_y_log10() +
    scale_color_viridis_d() +
    facet_grid(Type ~ .,
               scales = "free") +
    theme_tufte()

## Try simple logistic growth model

model <- stan_model("growth.stan")

cases <- df_wiki %>%
    filter(Series == "Cases") %>%
    arrange(Date) %>%
    spread(`Country/Region`, Value) %>%
    select(- Date, -Series) %>%
    as.matrix()
deaths <- df_wiki %>%
    filter(Series == "Deaths") %>%
    arrange(Date) %>%
    spread(`Country/Region`, Value) %>%
    select(- Date, -Series) %>%
    as.matrix()

fit <- sampling(model,
                data = list(T = nrow(cases), C = ncol(cases), cases = t(cases), deaths = t(deaths)),
                chains = 2)
