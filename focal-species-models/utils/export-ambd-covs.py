#-----------------------------------------------------------
# PREP HABITAT COVARIATES FOR AMERICAN BLACK DUCK
# USING GOOGLE EARTH ENGINE
#-----------------------------------------------------------

# Import packages
import ee
import math

# Initialize ee
ee.Initialize()

# Decisions
kernel_rad = 2500 # Radius (in meters) of smoothing kernel around each eBird checklist
scale = 250 # Scale of reducer calcs

# Bring in buffered range map as 'region'
region = ee.FeatureCollection('projects/XXX/ambd_range_buff_250km').first()

#----------------------------------------
#----------------------------------------
# FUNCTIONS

# Focal mean
def focal_mean(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit),
                                    reducer = ee.Reducer.mean())

# Focal median
def focal_median(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit),
                                    reducer = ee.Reducer.median())
    
# Focal SD
def focal_sd(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit),
                                    reducer = ee.Reducer.stdDev())

# Focal sum
def focal_sum(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit, False),
                                    reducer = ee.Reducer.sum())

# Focal count
def focal_count(image, radius, unit):
    return image.reduceNeighborhood(kernel = ee.Kernel.circle(radius, unit, False),
                                    reducer = ee.Reducer.count())

# Vector ruggedness measure
def compute_vrm(slope_img, aspect_img, radius, units):
    slope_sine = slope_img.sin()
    x_sum_sq = focal_sum(slope_sine.multiply(aspect_img.sin()), radius, units).pow(2)
    y_sum_sq = focal_sum(slope_sine.multiply(aspect_img.cos()), radius, units).pow(2)
    z_sum_sq = focal_sum(slope_img.cos(), radius, units).pow(2)
    n = focal_sum(ee.Image(1), radius, units)
    r = x_sum_sq.add(y_sum_sq).add(z_sum_sq).sqrt()
    vrm_img = ee.Image(1).subtract(r.divide(n))
    return vrm_img

def toFloat(img):
    return img.float()
#----------------------------------------
#----------------------------------------
# COLLECT COVARIATES

# Get multiband land use intensity layer
Lmulti = ee.Image('projects/XXX/Land-use-intensity-multiband-focal-sp-250m-20220123').unmask()

# Create nlcd water and wetlands layers
nlcd = ee.ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD").filter(ee.Filter.eq('system:index', '2016')).first().select('landcover')
water = nlcd.eq(11)
wet_herb = nlcd.eq(95)
wet_all = nlcd.eq(95).Or(nlcd.eq(90))

# Get tidal flats layer (Murray et al. 2019 Nature)
intertidal = ee.ImageCollection('UQ/murray/Intertidal/v1_1/global_intertidal')
intertidal = intertidal.filter(ee.Filter.eq('system:index', '2014-2016')).first().select('classification')
tidal = intertidal.eq(1).unmask()

# create Percent Cover layers for water, wetlands, and tidal flats
water_fsum = focal_sum(water, kernel_rad, 'meters')
water_fcount = focal_count(water, kernel_rad, 'meters')
water_pcov = water_fsum.divide(water_fcount).rename('water_pcov')

wet_herb_fsum = focal_sum(wet_herb, kernel_rad, 'meters')
wet_herb_fcount = focal_count(wet_herb, kernel_rad, 'meters')
wet_herb_pcov = wet_herb_fsum.divide(wet_herb_fcount).rename('wet_herb_pcov')

wet_all_fsum = focal_sum(wet_all, kernel_rad, 'meters')
wet_all_fcount = focal_count(wet_all, kernel_rad, 'meters')
wet_all_pcov = wet_all_fsum.divide(wet_all_fcount).rename('wet_all_pcov')

tidal_fsum = focal_sum(tidal, kernel_rad, 'meters')
tidal_fcount = focal_count(tidal, kernel_rad, 'meters')
tidal_pcov = tidal_fsum.divide(tidal_fcount).rename('tidal_pcov')

# Create some topographic viables
dsm = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2") # Digital surface model. Native res 30m
projElev = dsm.first().select(0).projection()
elev = dsm.select("DSM").mosaic().setDefaultProjection(projElev).rename('elevation') # Elevation
slope = ee.Terrain.slope(elev).rename('slope') # Slope
aspect = ee.Terrain.aspect(elev).rename('aspect') # Aspect
# Get vector ruggedness metric
window_radius = 1000
slopeRad = slope.multiply(ee.Number(math.pi).divide(180))
aspectRad = slope.multiply(ee.Number(math.pi).divide(180))
vrm = compute_vrm(slopeRad, aspectRad, window_radius, "meters").rename('vrm_1km')

# Average annual values for climatic variables from DAYMET v4
# Choose date range to consider (winters for focal years)
dateFilter_winter = ee.Filter.Or(ee.Filter.date('2013-11-01', '2014-04-01'),
                                 ee.Filter.date('2014-11-01', '2015-04-01'),
                                 ee.Filter.date('2015-11-01', '2016-04-01'),
                                 ee.Filter.date('2016-11-01', '2017-04-01'),
                                 ee.Filter.date('2017-11-01', '2018-04-01'),
                                 ee.Filter.date('2018-11-01', '2019-04-01'))
dateFilter_migration = ee.Filter.Or(ee.Filter.date('2014-03-18', '2014-06-28'), ee.Filter.date('2014-10-05', '2014-12-18'),
                                    ee.Filter.date('2015-03-18', '2015-06-28'), ee.Filter.date('2015-10-05', '2015-12-18'),
                                    ee.Filter.date('2016-03-18', '2016-06-28'), ee.Filter.date('2016-10-05', '2016-12-18'),
                                    ee.Filter.date('2017-03-18', '2017-06-28'), ee.Filter.date('2017-10-05', '2017-12-18'),
                                    ee.Filter.date('2018-03-18', '2018-06-28'), ee.Filter.date('2018-10-05', '2018-12-18')
                                    )
dateFilter_all = ee.Filter.date('2014-01-01', '2018-12-31')
swe = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_winter).select('swe') # DAYMET weather data. Native res 1km
tmin = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_migration).select('tmin') # DAYMET weather data. Native res 1km
tmax = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_migration).select('tmax') # DAYMET weather data. Native res 1km
prcp = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_migration).select('prcp') # DAYMET weather data. Native res 1km
swe_mean = swe.reduce(ee.Reducer.mean()).rename('swe')
tmin_mean = tmin.reduce(ee.Reducer.mean()).rename('tmin')
tmax_mean = tmax.reduce(ee.Reducer.mean()).rename('tmax')
prcp_mean = prcp.reduce(ee.Reducer.mean()).rename('prcp')

#------------------------
# SMOOTH CONTINOUS COVS AND COLLECT ALL

# Collect all continuous covariates into a single image
cont_covs = Lmulti.addBands([elev,
                             slope,
                             aspect,
                             vrm,
                             swe_mean,
                             tmin_mean,
                             tmax_mean,
                             prcp_mean])
# Smooth 'em
covs_mean = focal_mean(cont_covs, kernel_rad, 'meters')
# covs_median = focal_median(cont_covs, kernel_rad, 'meters')
covs_sd = focal_sd(cont_covs, kernel_rad, 'meters')
# Pull them all together with percent cover covariates
covs_all = ee.Image([covs_mean, covs_sd, water_pcov, wet_all_pcov, wet_herb_pcov, tidal_pcov])
# Convert all to float-32
covs_all = covs_all.float().clip(region)

task = ee.batch.Export.image.toDrive(image = covs_all,
                                     folder = 'ebd-covs',
                                     description = 'ambd-covs-250m',
                                     scale = scale,
                                     region = region.geometry(),
                                     maxPixels = 1e13,
                                     crs = projElev)
task.start()