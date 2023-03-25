#-----------------------------------------------------------
# DOWNLOAD HABITAT COVARIATE RASTERS FROM GOOGLE DRIVE TO LOCAL
#-----------------------------------------------------------

library(tidyverse)
library(googledrive)
library(lubridate)

# Set gdrive, download directory and save directory
gdrive_dir <- '~/ebd-covs' 
download_dir <- file.path("data/ebird/cov-rasters/")

# Authenticate in gdrive
options(httr_oob_default=TRUE)
drive_auth(email = 'YOUR-ADDRESS@gmail.com') # drive_find(n_max = 30)

# Get list of appropriate files
gdrive_files <- drive_ls(gdrive_dir)
gdrive_files_last <- gdrive_files %>% 
  hoist(drive_resource, "modifiedTime") %>% 
  group_by(name) %>% 
  arrange(desc(date(modifiedTime))) %>% 
  slice(1) %>%
  ungroup()

# Download 'em
mapply(function(id, name) {
  dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
  drive_download(
    as_id(id), 
    file.path(download_dir, name), 
    overwrite = TRUE
  )
}, id = gdrive_files_last$id, name = gdrive_files_last$name)


