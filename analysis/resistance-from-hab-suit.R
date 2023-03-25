#-----------------------------------------------------------
# TRANSFORM HABITAT SUITABILITY TO LANDSCAPE RESISTANCE
#-----------------------------------------------------------

library(tidyverse)
library(terra)

bird = "bobo"

# Read in habitat suitability raster
in1 = paste0("output/", bird, "-hab-suitability-full.tif")
hs <- rast(in1)

# Define resistance function (from Keeley et al. 2016 Landscape Ecol)
resist <- function(h, c){
  out <- 1000-999*((1-exp(-c*h))/(1-exp(-c)))
  return(out)
}

# Create three resistance surfaces based on scaling function above
out1 = paste0("output/", bird, "-resistance-linear-full.tif")
out2 = paste0("output/", bird, "-resistance-nonlin2-full.tif")
res1 <- resist(hs, 0.25) %>% 
  writeRaster(filename = out1, overwrite = TRUE) # roughly linear
res2 <- resist(hs, 8) %>% 
  writeRaster(filename = out2, overwrite = TRUE) # non-linear