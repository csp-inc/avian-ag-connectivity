#-----------------------------------------------------------
# PREP AND EXPORT HABITAT COVARIATES FOR BOBOLINK
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
region = ee.FeatureCollection('projects/XXX/bobo-range').first()

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

# Create nlcd water
nlcd = ee.ImageCollection("USGS/NLCD_RELEASES/2019_REL/NLCD").filter(ee.Filter.eq('system:index', '2016')).first().select('landcover')
water = nlcd.eq(11)
# create Percent Cover layers for water
water_fsum = focal_sum(water, kernel_rad, 'meters')
water_fcount = focal_count(water, kernel_rad, 'meters')
water_pcov = water_fsum.divide(water_fcount).rename('water_pcov')

# Get NLCD forest and shrub
forest = nlcd.eq(41).Or(nlcd.eq(42)).Or(nlcd.eq(43))
shrub = nlcd.eq(52)
# create Percent Cover layers for forest and shrub
forest_fsum = focal_sum(forest, kernel_rad, 'meters')
forest_fcount = focal_count(forest, kernel_rad, 'meters')
forest_pcov = forest_fsum.divide(forest_fcount).rename('forest_pcov')

shrub_fsum = focal_sum(shrub, kernel_rad, 'meters')
shrub_fcount = focal_count(shrub, kernel_rad, 'meters')
shrub_pcov = shrub_fsum.divide(shrub_fcount).rename('shrub_pcov')

# Get comprehensive grasslands from NLCD and FUT pasture
fut = ee.Image('projects/GEE_CSP/AFT_FUT/LCU_2015_v6')
grass = ee.Image(0).where(nlcd.eq(71).Or(nlcd.eq(81)), 1).where(fut.eq(1102).Or(fut.eq(2)), 1)
# Percent cover layer
grass_fsum = focal_sum(grass, kernel_rad, 'meters')
grass_fcount = focal_count(grass, kernel_rad, 'meters')
grass_pcov = grass_fsum.divide(grass_fcount).rename('grass_pcov')

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
# Choose date range to consider (breeding season for focal years)
dateFilter_breeding = ee.Filter.Or(ee.Filter.date('2014-05-01', '2014-09-14'),
                                    ee.Filter.date('2015-05-01', '2015-09-14'),
                                    ee.Filter.date('2016-05-01', '2016-09-14'),
                                    ee.Filter.date('2017-05-01', '2017-09-14'),
                                    ee.Filter.date('2018-05-01', '2018-09-14'))
tmin = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_breeding).select('tmin') # DAYMET weather data. Native res 1km
tmax = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_breeding).select('tmax') # DAYMET weather data. Native res 1km
prcp = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_breeding).select('prcp') # DAYMET weather data. Native res 1km
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
                             tmin_mean,
                             tmax_mean,
                             prcp_mean])
# Smooth 'em
covs_mean = focal_mean(cont_covs, kernel_rad, 'meters')
covs_sd = focal_sd(cont_covs, kernel_rad, 'meters')
# Pull them all together with percent cover covariates
covs_all = ee.Image([covs_mean, covs_sd, water_pcov, forest_pcov, shrub_pcov, grass_pcov])
# Convert all to float-32
covs_all = covs_all.float().clip(region)

task = ee.batch.Export.image.toDrive(image = covs_all,
                                     folder = 'ebd-covs',
                                     description = 'bobo-covs-250m',
                                     scale = scale,
                                     region = region.geometry(),
                                     maxPixels = 1e13,
                                     crs = projElev)
task.start()