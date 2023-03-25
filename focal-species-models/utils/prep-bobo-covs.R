#-----------------------------------------------------------
# PREP BOBOLINK COVARIATE RASTERS FOR ANALYSIS
#-----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

# Get covariate raster file names
tifs <- list.files('data/ebird/cov-rasters') 
tifs <- tifs[str_detect(tifs, "bobo")]
# Create list of rasters (need to split up task to avoid errors)
m_names <- c("Ag_mean","Urban_mean","Transport_mean","Energy_mean","elevation_mean",
             "slope_mean","aspect_mean","vrm_1km_mean","tmin_mean","tmax_mean","prcp_mean")
sd_names <- c("Ag_stdDev","Urban_stdDev","Transport_stdDev","Energy_stdDev","elevation_stdDev",
             "slope_stdDev","aspect_stdDev","vrm_1km_stdDev","tmin_stdDev","tmax_stdDev","prcp_stdDev")
p_names <- c("water_pcov","forest_pcov","shrub_pcov","grass_pcov")

m_list <- list()
sd_list <- list()
p_list <- list()
for(i in 1:length(tifs)){
  r_path <- paste0('data/ebird/cov-rasters/', tifs[i])
  r_temp <- rast(r_path)
  m_list[[i]] <- subset(r_temp, m_names)
  sd_list[[i]] <- subset(r_temp, sd_names)
  p_list[[i]] <- subset(r_temp, p_names)
}
# Convert to SpatialRasterCollection 
m_src <- terra::src(m_list)
# sd_src <- terra::src(sd_list)
p_src <- terra::src(p_list)
# Merge collections
m_covs <- terra::merge(m_src)
# sd_covs <- terra::merge(sd_src)
p_covs <- terra::merge(p_src)
# Combine merged rasters
# covs <- c(m_covs, sd_covs, p_covs)
covs <- c(m_covs, p_covs)

# Get buffered range map
rbuff <- st_read(dsn = 'data/focal-species-range-maps/bobolink-iucn-range',
                 layer = 'bobo_union_range')
# Get cell numbers for raster cells overlapping range
pred_cell_nums <- terra::cells(covs, vect(rbuff))[,2]
# Extract covariate values from across the covs raster for use in pred surface
pred_cov_temp <- terra::extract(covs, pred_cell_nums) 
pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
  as.data.frame() %>% 
  rename(longitude = x, latitude = y) 
pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
  drop_na()
save(pred_covs, file = 'data/ebird/bobo-pred-covs.rda')

# Get checklist data
bird = read.csv('data/ebird/ebd_bobo_2014_2018.csv', header = T) %>% 
  arrange(checklist_id) 
# Extract covariate values at sagr checklist locations
brd_vec <- bird %>% # Convert to sf and then to terra's SpatVect format
  st_as_sf(coords = c('longitude','latitude'), crs = st_crs(4326)) %>% 
  terra::vect()
loc_covs <- terra::extract(covs, brd_vec) # Extract raster values
bird_out <- cbind(bird, loc_covs)
save(bird_out, file = "data/ebird/ebd_bobo_2014_2018_wcovs.rda")
# write_csv(bird_out, file = "data/ebird/ebd_bobo_2014_2018_wcovs.csv")
