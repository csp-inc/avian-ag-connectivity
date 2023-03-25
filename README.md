# avian-ag-connectivity
### Organization: Conservation Science Partners, Inc
### Contact: Dr. Justin Suraci (justin@csp-inc.org)
<br>
<br>

This repo contains the scripts necessary to reproduce the analyses described in:
<br>

### *Suraci, JP, TG Mozelewski, CE Littlefield, T Nogeire McRae, A Sorensen, BG Dickson. (2023) Management of U.S. Agricultural Lands Differentially Affects Avian Habitat Connectivity. Land*
<br>

This work utilized eBird data (https://ebird.org/home) for three focal species - American Black Duck (bird ID used in species-specific scripts = "ambd"), Bobolink ("bobo"), and Greater Sage-Grouse ("sagr") -  to model the effects of agricultural land management and other factors on bird species habitat suitability and connectivity. The scripts provided here for preparing eBird data and using it to model habitat suitability were heavily inspired by code developed by Strimas-Mackey et al. (2020) in their book [Best Practices for Using eBird Data](https://cornelllabofornithology.github.io/ebird-best-practices/). The analysis workflow is outlined below.  
<br>
**Abstract:** Despite frequently being implicated in species declines, agricultural lands may nonetheless play an important role in connecting wildlife populations by serving as movement corridors or stopover sites between areas of high-quality habitat. For many North American bird species, agricultural intensification over the past half century has substantially impacted populations, yet recent studies have noted the potential for supporting avian biodiversity on agricultural lands through the promotion of functional connectivity. To support avian conservation efforts on agricultural lands across the United States, we used publicly available data from eBird to quantify and map the effects of agriculture on habitat suitability (using random forest models) and functional connectivity (via circuit theory) for three focal species that have experienced agriculture-linked declines or range contractions in recent decades: Greater Sage-grouse (*Centrocercus urophasianus*), American Black Duck (*Anas rubripes*), and Bobolink (*Dolichonyx oryzivorus*). Our analysis drew on novel, remotely sensed estimates of agricultural management intensity to quantify the effects of management practices on avian habitat and movement, revealing complex, species-specific relationships between agriculture and habitat value for the three focal species. Rangelands and croplands exhibited relatively high connectivity values for Greater Sage-grouse and Bobolink, respectively, mirroring these species’ strong habitat preferences for open sagebrush and cultivated grasslands. By contrast, American Black Duck migratory connectivity was low on all agricultural cover types. Mapping our model results across each species’ geographic range in the U.S. revealed key areas for agricultural management action to preserve high-quality habitat and connectivity, and we link these spatial recommendations to government incentive programs that can be used to increase wildlife-friendly management on U.S. agricultural lands.
<br>
<br>
## *Analysis Workflow* <br>
<br>

## Step 1 - Prepare eBird data and environmental covariates
* After downloading the eBird Basic Dataset locally, data for each focal species are prepped for analysis using `utils/prep-ebird-data.r`
* Environmental covariate rasters are prepped in the cloud using Google Earth Engine and exported to google drive via `utils/export-{bird ID}-covs.py`. These scripts require that all environmental covariates are uploaded as Earth Engine assets. Once covariate rasters are exported, they can be downloaded from Google Drive to local using `utils/get-covs-gdrive.R`
* Further processings of covariate rasters (using locally stored files) is handled by `utils/prep-{bird ID}-covs.R`

## Step 2 - Random forest models of habitat suitability
* Random forest models are run and validation procedures performed using `analysis/{bird ID}-rand-forest-model.R`
* Model predictions are used to generate a range-wide habitat suitability surface via `analysis/{bird ID}-hab-suitability-surface.R`

## Step 3 - Connectivity models
* `resistance-from-hab-suit.R` is used to transform habitat suitability surfaces from Step 2 in to landscape resistance surfaces for use in connectivity modeling.
* The Omniscape.jl package in Julia (see [Landau et al. 2021](https://joss.theoj.org/papers/10.21105/joss.02829)) is used to run omnidirectional connectivity models via `omniscape/{bird ID}_omni.ini`



