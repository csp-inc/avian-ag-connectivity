#-----------------------------------------------------------
# PREP GREATER SAGE-GROUSE COVARIATE RASTERS FOR ANALYSIS
#-----------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

# Get checklist data
bird = read.csv('data/ebird/ebd_sagr_2014_2018.csv', header = T) %>% 
  arrange(checklist_id) 
# Get buffered range map
rbuff <- st_read(dsn = 'data/focal-species-range-maps/sage-grouse-iucn-range',
                 layer = 'sagr_buff_120km')

# Get covariate rasters
covs1 <- rast('data/ebird/cov-rasters/sagr-covs-250m-0000000000-0000000000.tiff')
covs2 <- rast('data/ebird/cov-rasters/sagr-covs-250m-0000000000-0000006656.tiff')
covs <- terra::merge(covs1, covs2) # Merge the two raster stacks
pred_cell_nums <- terra::cells(covs, vect(rbuff))[,2] # Get cell numbers for raster cells overlapping sagr range

# Extract covariate values at sagr checklist locations
brd_vec <- bird %>% # Convert to sf and then to terra's SpatVect format
  st_as_sf(coords = c('longitude','latitude'), crs = st_crs(4326)) %>% 
  terra::vect()
loc_covs <- terra::extract(covs, brd_vec) # Extract raster values
bird_out <- cbind(bird, loc_covs) %>% 
  mutate(sage_pcov = sage_sum/sage_count) %>% 
  write_csv(file = "data/ebird/ebd_sagr_2014_2018_wcovs.csv")

# Extract covariate values from across the covs raster for use in pred surface
pred_cov_temp <- terra::extract(covs, pred_cell_nums) %>% 
  mutate(sage_pcov = sage_sum/sage_count) %>% 
  dplyr::select(-sage_sum, -sage_count)
pred_xy <- terra::xyFromCell(covs, pred_cell_nums) %>% 
  as.data.frame() %>% 
  rename(longitude = x, latitude = y) 
pred_covs <- cbind(pred_cov_temp, pred_xy) %>% 
  drop_na()
write_csv(pred_covs, file = 'data/ebird/sagr-pred-covs.csv')
