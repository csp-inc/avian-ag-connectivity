#####################################################################
# FIT SPATIALLY SUB-SAMPLED RANDOM FOREST MODEL FOR SAGE GROUSE

# This script implements the spatial bias and class imbalance
# corrections described in Robinson et al. (2017, Div. Dist.)
# in which non-detections are spatially sub-sampled and all
# detections are retained (i.e., no sub-sampling of detections).
# Random forest modeling, accuracy testing, and variable importance
# visualization follow Strimas-Mackey et al. 2020 (eBird best 
# practices book).
# Key decision:
# - Sub-sampling grid size: 5km
# - Train/test set size: 80%/20%
# - RF parameters: 500 trees, default mtry, balanced model
#####################################################################

library(tidyverse)
library(sf)
library(spData)
library(gridExtra)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)

# GET DATA AND FUNCTIONS
# Get species range map for creating sampling grid
rmap <- st_read(dsn = 'data/focal-species-range-maps/sage-grouse-iucn-range',
                layer = 'sagr_union_range') %>% 
  st_transform(5070)

# Read in species data and transform to appropriate crs
bird = read.csv('data/ebird/ebd_sagr_2014_2018_wcovs.csv', header = T) %>% 
  st_as_sf(coords = c('longitude','latitude'), crs = st_crs(4326)) %>% 
  st_transform(5070)

# Some functions for preparing sampling grid and sampling within grid
makeGrid <- function(poly, spacing){
  gTemp <- st_make_grid(poly, square = T, cellsize = c(spacing, spacing)) %>% 
    st_sf() # make grid using range polygon as bounding box
  int <- st_intersects(gTemp, poly, sparse = TRUE) # clean up by only keeping grid cells that intersect range
  keep <- lengths(int) > 0
  grid <- gTemp[keep,] %>% 
    mutate(cellNum = 1:nrow(.))
  return(grid)
}

cellSample <- function(points, grid, n = 1){
  int <- st_intersects(points, grid, sparse = FALSE)
  cellNum <- unlist(apply(int, 1, which))
  points$cellNum <- cellNum
  pt_sample <- points %>% 
    group_by(cellNum) %>% 
    sample_n(size = n) %>% 
    ungroup() 
  return(pt_sample)
}

#-------------------------------
#-------------------------------
# SPATIAL SUB SAMPLING

# ** IMPORTANT ** Set spatial grid size for subsampling
gridSize = 10000

# Make sampling grid, based on 
sagr_grid <- makeGrid(rmap, gridSize)

# Divide data into detections and non-detections
# Keep all detections, Keep one detection per grid cell
brd_pres <- bird %>% filter(species_observed)
brd_abs_sample <- bird %>% filter(species_observed==F) %>% 
  cellSample(sagr_grid) %>% 
  select(-cellNum)
brd_ss <- rbind(brd_pres, brd_abs_sample)


#-------------------------------
#-------------------------------
# RUN RANDOM FOREST MODEL

# Select columns to use in analysis and ditch others
brd_split <- brd_ss %>% 
  select(species_observed,year,day_of_year,Ag_mean,Energy_mean,Transport_mean,Urban_mean,
         aspect_mean,elevation_mean,slope_mean,vrm_1km_mean,swe_mean, tmin_mean, tmax_mean, 
         prcp_mean, 
         # Ag_stdDev,Energy_stdDev,Transport_stdDev,Urban_stdDev,aspect_stdDev, 
         # elevation_stdDev,slope_stdDev,vrm_1km_stdDev,swe_stdDev, tmin_stdDev, tmax_stdDev, 
         # prcp_stdDev, 
         sage_pcov, start_time_dec, duration_minutes,
         effort_distance_km, number_observers) %>%
  st_drop_geometry() %>% 
  drop_na()

# Split into training and test set
brd_split <- brd_split %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(brd_split, nrow) # show breakdown
detection_freq <- mean(brd_split$train$species_observed) # calculate detection rate
print(paste('detection frequency:', detection_freq))
  
# fit RF
brd_split$train$species_observed <- factor(brd_split$train$species_observed)
rf <- ranger(formula =  species_observed ~ ., 
             data = brd_split$train,
             num.trees = 500,
             importance = "impurity",
             probability = TRUE,
             replace = TRUE, 
             sample.fraction = c(detection_freq, detection_freq),
             verbose = T
             )

#-------------------------------
#-------------------------------
# FIT CALIBRATION MODEL

# make predictions on training data
occ_pred <- rf$predictions[, 2]
# convert the observered response back to a numeric value from factor
occ_obs <- brd_split$train$species_observed %>% 
  as.logical() %>% 
  as.integer()
rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) %>% 
  drop_na()

# fit calibration model
calibration_model <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                          gamma = 1.4,
                          data = rf_pred_train)

# calculate the average observed encounter rates for different 
# categories of estimated encounter rates 
average_encounter <- rf_pred_train %>%
  mutate(pred_cat = cut(rf_pred_train$pred, breaks = seq(0, 1, by=0.02))) %>%
  group_by(pred_cat) %>%
  summarise(pred = mean(pred), obs = mean(obs), checklist_count = n()) %>%
  ungroup()

# plot
cal_pred <- tibble(pred = seq(0, 1, length.out = 100))
cal_pred <- predict(calibration_model, cal_pred, type = "response") %>% 
  bind_cols(cal_pred, calibrated = .)
ggplot(cal_pred) +
  aes(x = pred, y = calibrated) +
  geom_line() +
  geom_point(data = average_encounter, 
             aes(x = pred, y = obs, size = sqrt(checklist_count)),
             show.legend = FALSE, shape = 1) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model")

#-------------------------------
#-------------------------------
# ASSESS MODEL ACCURACY

# predict on test data using calibrated model
p_fitted <- predict(rf, data = brd_split$test, type = "response")
# extract probability of detection
p_fitted <- p_fitted$predictions[, 2]
# calibrate
p_calibrated <- predict(calibration_model, 
                        newdata = tibble(pred = p_fitted), 
                        type = "response")
rf_pred_test <- data.frame(id = seq_along(p_calibrated),
                           # actual detection/non-detection
                           obs = brd_split$test$species_observed,
                           # uncalibrated prediction
                           fit = p_fitted,
                           # calibrated prediction
                           cal = p_calibrated) %>%
  # constrain probabilities to 0-1
  mutate(cal = pmin(pmax(cal, 0), 1)) %>% 
  drop_na()

# mean squared error (mse)
mse_fit <- mean((rf_pred_test$obs - rf_pred_test$fit)^2, na.rm = TRUE)
mse_cal <- mean((rf_pred_test$obs - rf_pred_test$cal)^2, na.rm = TRUE)

# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")

# calculate accuracy metrics: auc, kappa, sensitivity, specificity,
metrics_fit <- rf_pred_test %>% 
  select(id, obs, fit) %>% 
  presence.absence.accuracy(threshold = opt_thresh$fit, 
                            na.rm = TRUE, 
                            st.dev = FALSE)
metrics_cal <- rf_pred_test %>% 
  select(id, obs, cal) %>% 
  presence.absence.accuracy(threshold = opt_thresh$cal, 
                            na.rm = TRUE, 
                            st.dev = FALSE)

# VIEW METRICS
rf_assessment <- tibble(
  model = c("RF", "Calibrated RF"),
  mse = c(mse_fit, mse_cal),
  sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
  specificity = c(metrics_fit$specificity, metrics_cal$specificity),
  auc = c(metrics_fit$AUC, metrics_cal$AUC),
  kappa = c(metrics_fit$Kappa, metrics_cal$Kappa)
)
validTable <- knitr::kable(rf_assessment, digits = 3)
validTable

#-------------------------------
#-------------------------------
# PREDICTOR IMPORTANCE AND PARTIAL DEPENDENCE PLOTS

# IMPORTANCE
pi <- enframe(rf$variable.importance, "predictor", "importance")
# plot
ggplot(pi) + 
  aes(x = fct_reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, size = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5))

# PARTIAL DEPENDENCE
top_pred <- pi %>% 
  filter(!predictor %in% c("year")) %>% 
  top_n(n = 16, wt = importance) %>% 
  arrange(desc(importance))
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

# calculate partial dependence for each predictor
# map is used to iteratively apply calculate_pd to each predictor
pd <- top_pred %>% 
  mutate(pd = map(predictor, calculate_pd, model = rf, 
                  data = brd_split$train),
         pd = map(pd, ~ .[, c(2, 3)]),
         pd = map(pd, set_names, nm = c("value",  "encounter_rate"))) %>% 
  unnest(cols = pd)

# calibrate predictions
pd$encounter_rate <- predict(calibration_model, 
                             newdata = tibble(pred = pd$encounter_rate), 
                             type = "response") %>% 
  as.numeric()

# plot
ggplot(pd) +
  aes(x = value, y = encounter_rate) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(predictor), nrow = 4, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

#-------------------------------
#-------------------------------
# SAVE MODEL
save(list = c('brd_split','rf','pi','validTable','calibration_model'), file = 'output/sagr-random-forest-20220216.rda')
