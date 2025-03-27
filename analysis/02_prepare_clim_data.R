# This script prepares the climate data from Feni airport
# Based on Cliff's scripts

# Preamble ----------------------------------------------------------------

library(tidyverse) 


setNA <- function(x, na_val) {
  x[x == na_val] <- NA
  x
}

farToDeg <- function(x) {
  (x-32)/1.8
}

inTomm <- function(x) {
  x * 25.4
}

computeRH <- function(dewp, temp) {
  100 * (exp((17.625 * dewp)/(243.04 + dewp))/exp((17.625 * temp)/(243.04 + temp)))
}

computeSH <- function(relhum, stp, temp) {
  relhum / (0.263 * stp) * exp(17.67 * temp / (temp + 273.15 - 29.65))
}

# Load raw data -----------------------------------------------------------


chittaghon <- dir("data", 
                  pattern = "SHAH_AMANAT_INTL",
                  full.names = T) %>%
  map_df(function(x) {
    read_csv(x) %>% 
      janitor::clean_names() %>% 
      mutate(
        temp = map_dbl(temp, ~ as.numeric(.) %>% setNA(9999.9) %>% farToDeg()),
        dewp = map_dbl(dewp, ~ as.numeric(.) %>% setNA(9999.9) %>% farToDeg()),
        stp = map_dbl(stp, ~ as.numeric(.) %>% setNA(999.9)*100),
        min_temp = map_dbl(min, ~ as.numeric(.) %>% setNA(9999.9) %>% farToDeg()),
        max_temp = map_dbl(max, ~ as.numeric(.) %>% setNA(9999.9) %>% farToDeg()),
        prcp = map_dbl(prcp, ~ as.numeric(.) %>% setNA(99.99) %>% inTomm())
      ) %>% 
      rowwise() %>% 
      mutate(relhum = computeRH(dewp, temp),
             spechum = computeSH(relhum, dewp, temp)) %>% 
      ungroup() %>% 
      select(station, name, date, temp, min_temp, max_temp, prcp, relhum, spechum)
  }) %>% 
  mutate(date_num = as.numeric(date))

saveRDS(chittaghon, "generated_data/chittaghon_airport_weather_data.rds")

chittaghon %>% 
  pivot_longer(cols = c("temp", "min_temp", "max_temp", "prcp", "relhum", "spechum"),
               values_to = "value",
               names_to = "param") %>% 
  ggplot(aes(x = date, y = value)) +
  geom_line() +
  facet_wrap(~param, scales = "free") +
  theme_bw()

chittaghon_filled <- forecastML::fill_gaps(chittaghon, date_col = 3, 
                                           frequency = "1 days") %>% 
  as_tibble() %>% 
  mutate(date_num = as.numeric(date))


gam_fit <- mgcv::gam(temp ~ s(date_num, k = 20), data = chittaghon)
chittaghon_filled$temp_pred <- predict(gam_fit, chittaghon_filled)

chittaghon_filled <- chittaghon_filled %>% 
  mutate(temp = coalesce(temp, temp_pred))

saveRDS(chittaghon_filled, "generated_data/chittaghon_airport_weather_data_filled.rds")


# Precipitation data from IMERG -------------------------------------------

# Longitude and latitude of study area
lonlat <- data.frame(lon = 91.6809 * c(0.999, 1.001), lat = 22.6171 * c(0.999, 1.001))

dates <- c("2021-01-01", "2022-12-31")

imerg_precip <- chirps::get_imerg(lonlat, dates)

saveRDS(imerg_precip, "generated_data/imerg_precipitation_sitakunda.rds")

