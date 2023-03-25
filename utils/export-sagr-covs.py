#-----------------------------------------------------------
# PREP AND EXPORT HABITAT COVARIATES FOR GREATER SAGE-GROUSE
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
Lmulti = ee.Image('projects/XXX/Land-use-intensity-multiband-focal-sp-250m-20220123')

# Create sagebrush presence/absence layer from landfire
lf = ee.Image('projects/XXX/landfire-v20-2016-existing-veg-type') # Landfire layer. Native res 30m
# UPDATE - REMOVE 7087 (BURSAGE)
sage_classes = ee.List([7064,7124,9030,7079,7080,7125,7126,7072]) # All classes tagged with sagebrush
sage = lf.remap(sage_classes, ee.List.repeat(1, sage_classes.size())).unmask().rename('sage')
# sage_area = sage.multiply(ee.Image.pixelArea()).rename('sage_area')
sage_fsum = focal_sum(sage, kernel_rad, 'meters')
sage_fcount = focal_count(sage, kernel_rad, 'meters')
sage_pcov = sage_fsum.divide(sage_fcount)

# Create some topographic iables
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
# UPDATE - INCLUDE "ALL" DATA FILTER FOR ALL COVS OTHER THAN SWE
dateFilter_winter = ee.Filter.Or(ee.Filter.date('2014-01-01', '2014-04-01'),
                          ee.Filter.date('2014-11-01', '2015-04-01'),
                          ee.Filter.date('2015-11-01', '2016-04-01'),
                          ee.Filter.date('2016-11-01', '2017-04-01'),
                          ee.Filter.date('2017-11-01', '2018-04-01'),
                          ee.Filter.date('2018-11-01', '2018-12-31'))
dateFilter_all = ee.Filter.date('2014-01-01', '2018-12-31')
swe = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_winter).select('swe') # DAYMET weather data. Native res 1km
tmin = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_all).select('tmin') # DAYMET weather data. Native res 1km
tmax = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_all).select('tmax') # DAYMET weather data. Native res 1km
prcp = ee.ImageCollection('NASA/ORNL/DAYMET_V4').filter(dateFilter_all).select('prcp') # DAYMET weather data. Native res 1km
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
# Pull them all together with sage covs
covs_all = ee.Image([covs_mean, covs_sd, sage_fsum, sage_fcount]) # EXCLUDES MEDIAN SMOOTH
# Convert all to float-32
covs_all = covs_all.float()

region = ee.FeatureCollection('projects/XXX/sagr-range-buff-120km').first()
task = ee.batch.Export.image.toDrive(image = covs_all,
                                     folder = 'ebd-covs',
                                     description = 'sagr-covs-250m',
                                     scale = scale,
                                     region = region.geometry(),
                                     maxPixels = 1e13,
                                     crs = projElev)
task.start()