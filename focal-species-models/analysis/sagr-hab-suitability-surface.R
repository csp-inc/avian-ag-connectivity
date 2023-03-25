#-----------------------------------------------------------
# PREP GREATER SAGE GROUSE HABITAT SUITABILITY SURFACE
# FROM RANDOM FOREST MODEL OUTPUT
#-----------------------------------------------------------

library(tidyverse)
library(sf)
library(spData)
library(gridExtra)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(raster)
library(rgdal)

# GET DATA / MODEL
# Load random forest model and associated data
load(file = 'output/sagr-random-forest-20220216.rda', verbose = TRUE)
# Read in prediction dataset
pred_covs <- read_csv(file = 'data/ebird/sagr-pred-covs.csv')

#-------------------------------
#-------------------------------
# DEAL WITH EFFORT/DETECTION COVS

# function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid
  rng <- range(data[[predictor]], na.rm = TRUE)
  x_grid <- seq(rng[1], rng[2], length.out = x_res)
  grid <- data.frame(covariate = predictor, x = x_grid, 
                     stringsAsFactors = FALSE)
  names(grid) <- c("covariate", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p <- predict(model, data = grid)
  
  # summarize
  pd <- grid[, c("covariate", predictor)]
  names(pd) <- c("covariate", "x")
  pd$pred <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, covariate, x) %>% 
    dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pd)
}

# Find peak time of day from partial dependence 
# (based on code in ebird best practices book)
pd_time <- calculate_pd("start_time_dec",
                        model = rf, 
                        data = brd_split$train,
                        # make estimates at 30 minute intervals
                        # using a subset of the training dataset
                        x_res = 2 * 24, n = 1000) %>% 
  transmute(start_time_dec = x, encounter_rate = pred)
# hours with at least 1% of checklists
search_hours <- brd_split$train %>% 
  mutate(hour = floor(start_time_dec)) %>%
  count(hour) %>% 
  mutate(pct = n / sum(n)) %>% 
  filter(pct >= 0.01)
# constrained peak time
t_peak <- pd_time %>% 
  filter(floor(start_time_dec) %in% search_hours$hour) %>% 
  top_n(1, wt = desc(start_time_dec)) %>% 
  pull(start_time_dec)
t_peak

# Find peak day of year from partial dependence 
# (based on code in ebird best practices book)
pd_day <- calculate_pd("day_of_year",
                        model = rf, 
                        data = brd_split$train,
                        # make estimates at 30 minute intervals
                        # using a subset of the training dataset
                        x_res = 365, n = 1000) %>% 
  transmute(day_of_year = x, encounter_rate = pred)
# hours with at least 1% of checklists
search_days <- brd_split$train %>% 
  mutate(day = floor(day_of_year)) %>%
  count(day) %>% 
  mutate(pct = n / sum(n)) %>% 
  filter(pct >= 0.001)
# constrained peak time
day_peak <- pd_day %>% 
  filter(floor(day_of_year) %in% search_days$day) %>% 
  top_n(1, wt = (encounter_rate)) %>% 
  pull(day_of_year) %>% floor()
day_peak

# Add detection/effort covs to pred_covs
pred_covs_eff <- pred_covs %>% 
  mutate(year = 2016,
         day_of_year = day_peak,
         start_time_dec = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1)

#-------------------------------
#-------------------------------
# PREDICT FROM RF MODEL
# predict
pred_rf <- predict(rf, data = pred_covs_eff, type = "response")
pred_rf <- pred_rf$predictions[, 2]
# apply calibration models
pred_rf_cal <- predict(calibration_model, 
                       data.frame(pred = pred_rf), 
                       type = "response")
# add to prediction surface
pred_er <- bind_cols(pred_covs_eff, encounter_rate = pred_rf_cal) %>% 
  select(latitude, longitude, encounter_rate) %>% 
  mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))

#-------------------------------
#-------------------------------
# CREATE HABITAT SUITABILITY SURFACE

# Convert prediction dataset to sf and then terra SpatVector
pred_sp <- pred_er %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(5070) 

# Create template raster
r_temp <- raster()
raster::extent(r_temp) <- raster::extent(pred_sp)
raster::res(r_temp) <- 250
raster::crs(r_temp) <- st_crs(pred_sp)$wkt

# free up some space
rm(list = c('pred_covs', 'pred_covs_eff', 'pred_rf', 'pred_rf_cal'))
gc()

# Rasterize
pred_rast <- raster::rasterize(pred_sp, r_temp, field = 'encounter_rate')

# save the raster
raster::writeRaster(pred_rast, 'output/sagr-hab-suitability.tif', overwrite = TRUE)
