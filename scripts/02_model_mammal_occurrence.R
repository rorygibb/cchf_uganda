



# ====================== MAMMAL HOST BRT DISTRIBUTION MODELS =============================


# fits, evaluates and predicts spatially from ensembles of boosted regression tree SDMs for mammal CCHF hosts in Africa
# currently Lepus spp.
# creates a folder in output/model_outputs for specific species/run


# ------------- dependencies and housekeeping --------------

# root dir
setwd("C:/Users/roryj/Documents/PhD/project_cchf_africa/analysis/")
#setwd("C:/Users/roryj/Dropbox/CrimeanCongo_Uganda/")

library(gbm)
library(rgbif)
library(countrycode)
library(data.table)
library(ggplot2)
library(maptools)
library(rgdal)
library(raster)
library(magrittr)
library(viridis)
library(rgeos)
#library(embarcadero)
library(sf)
library(ENMeval)
library(ggforce)
library(dplyr)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)



# ======================== MODELLING FUNCTIONS =============================



# ------------------------ data filtering and preparation -----------------------

### filter GBIF occurrences to remove erroneous/issues
#' @param occ dataframe of occurrences from GBIF

filterGBIFOccurrences = function(occ){
  
  # set id
  occ$localID = 1:nrow(occ)
  occ$countryName = countrycode(occ$countryCode, origin="iso2c", destination="country.name.en")
  
  ## 1: remove records with imperfect taxonomic match 
  occ = occ[ !grepl("TAXON", occ$issue), ]
  
  ## 2: only records with coordinate
  occ = occ[ occ$hasCoordinate == "true" & occ$hasGeospatialIssues == "false", ]
  
  # identify and remove records outside Africa using wrld_simpl
  data(wrld_simpl)
  shp = wrld_simpl[ wrld_simpl$REGION == 2, ]
  spdf = occ
  coordinates(spdf) = ~decimalLongitude+decimalLatitude
  proj4string(spdf) = proj4string(shp)
  yn = as.vector(over(spdf, as(shp, "SpatialPolygons")))
  occ = occ[ !is.na(yn), ]
  
  ## 3: remove records with no year and from before 1960
  occ = occ[ !is.na(occ$year), ]
  occ = occ[ occ$year >= 1930, ]
  
  ## 4: return coordinate precision (signif figures after decimal point) 
  precCheck = function(x){
    if(!is.numeric(x)) return(NA)
    if(!grepl("[.]", as.character(x))) return(0)
    y = strsplit(as.character(x), "[.]")
    return(nchar(y[[1]][2]))
  }
  
  ## 5. remove inaturalist observations: difficult with hare taxonomy
  occ = occ[ occ$institutionCode != "iNaturalist", ]
  
  # wrapper to exclude low precision coords
  precFilter = function(df, threshold=3){
    plat = sapply(df$decimalLatitude, precCheck)
    plon = sapply(df$decimalLongitude, precCheck)
    binary = ifelse(plat >= threshold | plon >= threshold, 1, 0)
    return( df[ binary > 0, ])
  }
  
  ## 5: exclude low precision coordinates (at least 2 decimal GPS points)
  occ1 = precFilter(occ, threshold=2)
  #plot(shp); points(occ1$decimalLongitude, occ1$decimalLatitude, cex=0.5, col="red", pch=16)
  
  ## 6. exclude duplicate lon-lat points
  occ1 = occ1[ !duplicated( occ1[ , c("decimalLongitude", "decimalLatitude")]), ]
  
  # 7. save
  gbif = occ1[ , c("class", "order", "family", "species", "countryCode", "countryName", "locality", "decimalLongitude", "decimalLatitude", "year", "basisOfRecord", "source", "datasetKey", "issue") ]
  names(gbif) = c("Class", "Order", "Family", "Species", "ISO2C", "Country", "Locality", "x", "y", "Year", "BasisofRecord", "RecordedBy", "datasetKey", "Issues")
  gbif$Source = "GBIF"
  gbif$Year = as.character(gbif$Year)
  return(gbif)
}


### thinOccurrences: thin occurrence points to reduce spatial autocorrelation
# remove all duplicate points falling within specified spatial resolution cells

#' @param occ spatial points dataframe of occurrences
#' @param res resolution of raster to use for thinning (i.e. how much to thin)
#' @param crs specify CRS of thinning raster (defaults to WGS84)
#' @param raster_ext specify extent of thinning raster (defaults to extent of SpatialPoints* object)

thinOccurrences = function(occ, res=0.0416, crs=NA, raster_ext=NA){
  
  # check class of occurrences
  if( ! grepl("SpatialPoints*", class(occ)) ) return("Error: occ must be a SpatialPoints* object")
  
  # set default CRS and raster ext
  if(is.na(crs)){ crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" }
  if(is.na(raster_ext)){ raster_ext = extent(occ) + 10  }
  
  # create raster
  r = raster(ext = raster_ext, crs = crs, resolution = res)
  #values(r) = sample(1:10, ncell(r), replace=T); plot(r); plot(occ, add=T)
  
  # find cell numbers in which each lat-lon point falls and remove duplicate points falling within same cell
  rcell = cellFromXY(r, occ)
  occ1 = occ[ which(!duplicated(rcell)), ]
  
  # return thinned occurrence points
  return(occ1)
}


### generatePseudos: generate pseudoabsences according to three possible weighting regimes:
# 1. no weighting: random across extent of study area
# 2. 'proximity': weighted towards areas with higher density of occurrences, with sigma controlling weighting bandwidth
# 3. 'distance': weighted away from areas with higher density of occurrences, with sigma controlling weighting bandwidth
# probably want to use larger sigmas for probability and smaller sigmas for distance (depending on choice)

#' @param occ spatial points dataframe of occurrences to use for weighting / extent generation
#' @param mask_raster raster to constrain points generation (i.e. land/ocean)
#' @param num_points number of pseudoabsences to generate (defaults to same number as occurrences)
#' @param weighting either "none" (randomly generated within bounds), "proximity" (weighted towards occurrences), "distance" (weighted away from occurrences)
#' @param ext extent object: to constrain points to specific region of mask_raster
#' @param sigma specifies bandwidth of kernel density: higher bandwidth reduces spatial clustering of pseudoabsences

generatePseudos = function(occ, mask_raster, num_points=NA, weighting="none", ext=NULL, sigma=NULL){
  
  # housekeeping
  if(!weighting %in% c("none", "proximity", "distance")) return("Error! Check weighting argument")
  if(is.na(num_points)) num_points = nrow(occ)
  if(!is.null(ext)) mask_raster = crop(mask_raster, ext) # if specified, crop raster to specified extent
  if(is.null(sigma)) sigma = 2.5 # default bandwidth: 2.5
  
  # generate points (unweighted)
  if(weighting == "none"){
    p = as.data.frame(dismo::randomPoints(mask_raster, n=num_points, p=occ))
    p$type = "pseudoabs_random"
    #plot(mask_raster); points(p$x, p$y, pch=16, cex=0.5, col="black")
  }
  
  # generate points (weighted towards sampled points)
  if(weighting == "proximity"){
    
    # estimate kernel density of points across raster
    pd = pointsDensityRaster(as.vector(occ@coords[, 1]), as.vector(occ@coords[, 2]), 
                             xrange = c(xmin(mask_raster), xmax(mask_raster)), 
                             yrange = c(ymin(mask_raster), ymax(mask_raster)), 
                             dimx = round(ncol(mask_raster)/2), round(dimy = nrow(mask_raster)/2), 
                             sigma=sigma, 
                             mask_raster=mask_raster)
    
    # scale to between zero and one (probability weighting) and set <0 values to
    pdx = pd / cellStats(pd, "max")
    values(pdx)[ values(pdx)<0 ] = 0
    p = as.data.frame(dismo::randomPoints(pdx, n=num_points, p=occ, prob=TRUE))
    p$type = "pseudoabs_proxweight"
    #plot(pdx); points(p$x, p$y, pch=16, cex=0.5, col="black")
  }
  
  # generate points (weighted away from sampled points)
  if(weighting == "distance"){
    
    # estimate kernel density of points across raster
    pd = pointsDensityRaster(as.vector(occ@coords[, 1]), as.vector(occ@coords[, 2]), 
                             xrange = c(xmin(mask_raster), xmax(mask_raster)), yrange = c(ymin(mask_raster), ymax(mask_raster)), 
                             dimx = ncol(mask_raster)/2, dimy = nrow(mask_raster)/2, 
                             sigma=sigma, 
                             mask_raster=mask_raster)
    
    # invert probability weighting (so points are spatially biased away from occurrences)
    pdx = 1- (pd / cellStats(pd, "max"))
    p = as.data.frame(dismo::randomPoints(pdx, n=num_points, p=occ, prob=TRUE))
    p$type = "pseudoabs_distweight"
    #plot(pdx); points(p$x, p$y, pch=16, cex=0.5, col="black")
  }
  
  return(p)
}


### extractEnviroData: function to extract environmental data for xy points from raster layer/stack/brick
# use to extract climatic data for fitting envelope

#' @param x vector of points on x dimension (longitude)
#' @param y vector of points on y dimension (latitude)
#' @param env_raster rasterlayer/stack/brick of environmental predictors to extract datafrom

extractEnviroData = function(x, y, env_raster){
  
  # check arguments and throw errors
  if(!is.null(env_raster) & !class(env_raster) %in% c("RasterLayer", "RasterStack", "RasterBrick")){ return("Error: env_raster must be of class raster") }
  if(length(x) != length(y)){ return("Error: x and y must be the same length") }
  if(length(is.na(x)) != length(is.na(y))){ return("Error: different number of NA vals in x and y") }
  
  # extract environmental data for all points
  coords = as.data.frame(cbind(x, y))
  coordinates(coords) = ~x+y
  projection(coords) = projection(env_raster)
  ed = as.data.frame(raster::extract(env_raster, coords))
  if(ncol(ed) == 1){ names(ed) = names(env_raster) }
  
  # return
  return(ed)
}


### processCHELSA: pre-processes CHELSA climatologies bioclimatic data 
# renames all layers for modelling

#' @param env_rasters RasterLayer/Stack/Brick from CHELSA climatologies data

processCHELSA = function(env_rasters){
  
  # data source and documentation: https://www.nature.com/articles/sdata2017122
  # http://chelsa-climate.org/
  
  # renames variables
  old_names = names(env_rasters)
  s = strsplit(old_names, "_")
  new_names = paste("bio", lapply(s, "[[", 3), "_", lapply(s, "[[", 5), "_", lapply(s, "[[", 6), sep="")
  names(env_rasters) = new_names
  
  # rescale bio03 and bio04 by another factor of 10 (*100 in original doc)
  env_rasters[[ grep("bio03", names(env_rasters)) ]] = env_rasters[[ grep("bio03", names(env_rasters)) ]]/10
  env_rasters[[ grep("bio04", names(env_rasters)) ]] = env_rasters[[ grep("bio04", names(env_rasters)) ]]/10
  
  # set precip seasonality values for areas with no precip (i.e. mean == 0) to max seasonality
  aa = env_rasters[[ grep("precip_annual", names(env_rasters)) ]]
  values(env_rasters[[ grep("precip_seas", names(env_rasters)) ]])[ which(values(aa) == 0) ] = cellStats(env_rasters[[ grep("precip_seas", names(env_rasters)) ]], max)
  env_rasters
}



# ---------------------- Functions to fit SDM ensembles ----------------------

### generateDataset: generates a list of datasets for fitting a bootstrap ensemble of SDMs 
# n pseudoabsence regimes (randomly selected) * k folds (either random or spatially structured)
# returns list of class 'nicheData' for fitting models

#' @param pres_data dataframe of occurrence points (presences)
#' @param num_pseudo specify number of pseudoabsence regimes
#' @param k specify k for kfold train-test if random
#' @param checkerboard specify if partitioning via checkerboard spatial approach
#' @param checker_type if checkerboard, eiher 1 (splits into 2 approx equal size checkerboard folds) or 2 (splits into 4 via nested checkerboards)
#' @param mask_raster raster to constrain pseudo points (e.g. to land)
#' @param num_pseudo_points specify number of pseudoabsence points per regime

generateDataset = function(pres_data, num_pseudo=5, k=NA, checkerboard=FALSE, checker_type=2, mask_raster, num_pseudo_points=200){
  
  # for storing results
  result = list()
  
  # create n sets of pseudoabsences, combine with occurrence data, and partition into 
  for(i in 1:num_pseudo){
    
    # return report
    print(paste("Generating pseudoabsences regime ", i, sep=""))
    
    # create random pseudoabsences and combine with presences
    pa = generatePseudos(pres_data, num_points=num_pseudo_points, mask_raster=mask_raster, weighting="none", ext=extent(mask_raster))
    occx = data.frame(x = pres_data@coords[, 1], y = pres_data@coords[, 2], type = "presence")
    df = rbind(occx, pa)
    df$occurrence = ifelse(df$type == "presence", 1, 0)
    df$pseud_regime = i
    
    # extract environmental variables
    ed = extractEnviroData(df$x, df$y, envx)
    df = cbind(df, ed)
    
    # if specified, random kfold partition
    # random kfold dataset into k training/test sets
    if(!is.na(k)){
      
      df$k = kfold(df, k=k)
      resx = vector("list", length=k)
      for(j in 1:k){
        dfx = df
        dfx$train = ifelse(dfx$k == j, 0, 1)
        resx[[j]] = dfx
      }
      # add to result
      result = append(result, resx)
    }
    
    # OR checkerboard kfold dataset
    # stricter rules of geographical partitioning to help account for spatial bias
    if(checkerboard==TRUE){
      
      # checkerboard partition into 4 datasets (aggregation factor of 2)
      # getcheckerboard2 function from ENMeval
      occ.pts = df[ df$occurrence == 1,]
      abs.pts = df[ df$occurrence == 0, ]
      if(checker_type == 1){ checker_points = ENMeval::get.checkerboard1(cbind(occ.pts$x, occ.pts$y), envx, cbind(abs.pts$x, abs.pts$y), aggregation.factor=2) }
      if(checker_type == 2){ checker_points = ENMeval::get.checkerboard2(cbind(occ.pts$x, occ.pts$y), envx, cbind(abs.pts$x, abs.pts$y), aggregation.factor=2) }
      occ.pts$k = checker_points$occ.grp
      abs.pts$k = checker_points$bg.grp
      df = rbind(occ.pts, abs.pts)
      
      # create train-test datasets
      resx = vector("list", length=max(df$k))
      
      for(j in 1:length(resx)){
        dfx = df
        dfx$train = ifelse(dfx$k == j, 0, 1)
        resx[[j]] = dfx
      }
      result = append(result, resx)
    }
    
  } # end of pseudos loop
  
  # dataset
  class(result) = "nicheData"
  return(result)
}

### fitModelEnsemble: combined function to fit full model ensemble of specified type

#' @param occ dataframe of occurrence points (presences)
#' @param species name of species being modelled
#' @param modeltype type of model (GLM, GAM, BRT or BioClim)
#' @param num_pseudo specify number of pseudoabsence regimes
#' @param mask_raster specify raster to mask pseudo points extent
#' @param k if random kfold: specify kfold for train-test
#' @param checkerboard if checkerboard cross-val: specify T/F
#' @param checker_type if checkerboard, specify type of checkerboard
#' @param learning_rate if BRT, specify learning rate for tree fitting (compelxity fixed at 5)
#' @param num_pseudo_points specify number of pseudoabsence points per regime
#' @param filepath directory to save ensemble model outputs

fitModelEnsemble = function(occ, species, modeltype, num_pseudo=5, mask_raster, k=NA, checkerboard=FALSE, checker_type=2, learning_rate=0.005, num_pseudo_points=200, filepath="./output/model_outputs/climate_ensemble/fitted_models/"){
  
  # generate datasets with specified regime
  print("Generating dataset...")
  sdm_data = generateDataset(occ, num_pseudo=num_pseudo, k=k, checkerboard=checkerboard, checker_type=checker_type, num_pseudo_points=num_pseudo_points, mask_raster=mask_raster)
  print(paste(class(sdm_data), "length", length(sdm_data)), sep=" ")
  
  # create directory for saving ensemble outputs
  identifier = sample(1:10^6, 1)
  dir_name = paste(species, modeltype, Sys.Date(), identifier, sep="_")
  model_path = paste(filepath, dir_name, "/", sep="")
  dir.create(model_path)
  
  # save occurrences
  dir.create(paste(model_path, "occurrence_points/", sep=""))
  write.csv(as.data.frame(occ), file=paste(model_path, "occurrence_points/", species, "_occurrences.csv", sep=""), row.names=FALSE)
  
  ### fitBRT: fits, evaluates and thresholds bernoulli boosted regression trees model
  fitBRT = function(x){
    
    # reporting
    cat(paste(x, "...", sep=""))
    
    # create unique identifier
    idx = sample(1:10^8, 1)
    
    # data (rearrange to name response as y)
    dfx = sdm_data[[x]]
    dfx2 = cbind(data.frame(y = dfx$occurrence), dfx[ , which(names(dfx) %in% predictors)])
    dfx_train = dfx2[ dfx$train == 1, ] 
    dfx_test = dfx2[ dfx$train == 0, ]
    
    # fit BRT model to training subset (suppress messages)
    mod = gbm.step(data = dfx_train, 
                   gbm.x = 2:ncol(dfx_train), 
                   gbm.y = 1, 
                   family = "bernoulli", 
                   tree.complexity=5, 
                   learning.rate=learning_rate, # have increased learning rate for ticks as so many occurrences
                   bag.fraction=0.5,
                   silent=TRUE)
    
    # evaluate and threshold based on test dataset
    e = evaluate(p = dfx_test[ dfx_test$y == 1,], a = dfx_test[ dfx_test$y == 0,], mod, n.trees=mod$gbm.call$best.trees)
    invlogit = function(x) exp(x)/(1+exp(x))
    t = invlogit(threshold(e, 'spec_sens'))
    
    # create and save model result list
    res = list(modeltype = "BRT", data = dfx, species=species, model = mod, auc = e@auc, threshold = t, unique_id = idx)
    save(res, file=paste(model_path, species, round(e@auc, 3), idx, ".R", sep="_"))
    rm(res)
    
    # create return df
    result_df = data.frame(date=Sys.Date(), ensemble_id = identifier, species=species, 
                           model_id = idx, modeltype = "BRT", 
                           kfold=k, checkerboard=checkerboard, checker_type=checker_type, 
                           mum_occurrences=nrow(occ), num_pseudo=num_pseudo_points,
                           num_train=nrow(dfx_train), num_test=nrow(dfx_test), 
                           auc=round(e@auc, 4), threshold=t, 
                           path=paste(model_path, species, round(e@auc, 3), idx, ".R", sep="_"))
    return(result_df)
  }
  
  # fit model ensemble
  print(paste("Fitting model ensemble of type", modeltype, sep=" "))
  if(modeltype == "GLM"){ ensemble = lapply(1:length(sdm_data), fitGLM) 
  } else if(modeltype == "GAM"){ ensemble = lapply(1:length(sdm_data), fitGAM)
  } else if(modeltype == "BRT"){ ensemble = lapply(1:length(sdm_data), fitBRT)
  } else if(modeltype == "Maxent"){ ensemble = lapply(1:length(sdm_data), fitMaxent)
  } else(return("Error: model must be either GLM, GAM, BRT, Maxent or BioClim"))
  
  # save and return result
  ens_metadata = do.call(rbind.data.frame, ensemble)
  write.csv(ens_metadata, paste(model_path, species, identifier, "ensemblemetadata.csv", sep="_"), row.names=FALSE)
  return(ens_metadata)
}




# ======================= LOAD DATA FOR MODELLING =========================


# --------- africa shp ----------

shp0 = st_read("./data/project_shapefiles/africa_admin0.shp")



# --------- environmental data -----------

# climatic predictors (already cropped to required extent)
env = raster::stack(list.files("./data/chelsa_climatologies/africa_rescaled_aggregated/", pattern=".tif", full.names=T))
env = processCHELSA(env)

# NDVI (Copernicus 1km original; derived from long-term statistics data 1992-2017)
nd = stack(list.files("./output/ndvi_copernicus/africa_summary/", full.names = TRUE, pattern=".tif"))
nd = resample(nd, env, "bilinear")
names(nd) = c("NDVI_AnnualMean", "NDVI_AnnualSD", "NDVI_AnnualMin")
env = stack(env, nd)

### create mask for pseudoabsences 
# geog_shp = rgdal::readOGR("./data/project_shapefiles/africa_shp.shp")
# 
# # combine geog shapefile into single polygon (ff needed)
# geog_shp2 = geog_shp; geog_shp2$xxx = 1
# geog_shp2 = maptools::unionSpatialPolygons(geog_shp2, ID=geog_shp2$xxx)

# create geographical raster for masking pseudoabsences 
# geo_mask = raster::raster(ext=extent(env), res=0.0416, crs = crs(env))
# geo_mask = raster::rasterize(geog_shp2, geo_mask)
# writeRaster(geo_mask, "./data/shapefiles/africa_maskraster/africa_mask.tif", format="GTiff", overwrite=T)

# reads in mask that was created using above code
geo_mask = raster::raster("./data/africa_maskraster/africa_mask.tif")

# 8 predictors for models
predictors = c("bio01_temp_annualmean",
               "bio04_temp_seasonality",
               "bio12_precip_annualtotal",
               "bio15_precip_seasonality",
               "bio05_temp_maxwarmestmonth",
               "bio06_temp_mincoldestmonth",
               "NDVI_AnnualMean", 
               "NDVI_AnnualSD",
               "NDVI_AnnualMin")

# subset rasters to only variables of interest
envx = env[[ which(names(env) %in% predictors)]]
rm(env)
envx = crop(envx, geo_mask)






# ----------------- load species occurrences --------------------

# read and filter for each species
h1 = filterGBIFOccurrences(read.delim("./data/host_occurrences/lepus_capensis/occurrence.txt"))
h2 = filterGBIFOccurrences(read.delim("./data/host_occurrences/lepus_microtis/occurrence.txt"))
h3 = filterGBIFOccurrences(read.delim("./data/host_occurrences/lepus_saxatilis/occurrence.txt"))

# combine
hosts = do.call(rbind.data.frame, list(h1, h2, h3))

# read extra data
ex = read.csv("./data/host_occurrences/rg_additional_literature.csv", stringsAsFactors = FALSE)
head(ex)
ex$Class = "Mammalia"
ex$Order = "Lagomorpha"
ex$Family = "Leporidae"
names(ex)[1] = "Species"
ex$ISO2C = "UGA"
ex$Country = "Uganda"
ex$Locality = NA
names(ex)[ names(ex) == "Longitude" ] = "x"
names(ex)[ names(ex) == "Latitude" ] = "y"
names(ex)[ names(ex) == "Year_collected" ] = "Year"
names(ex)[ names(ex) == "Reference" ] = "Source"
names(ex)[ names(ex) == "Notes" ] = "Issues"
ex = ex[ , -which(names(ex) == "URL")]
ex$BasisofRecord = "FieldObs"
ex$RecordedBy = "Thorne_2009"
ex$datasetKey = NA

# combine
hosts = rbind(hosts, ex)

# ggplot() +
#   geom_sf(data=shp0, colour="grey30", fill="grey95", size=0.25) +
#   theme_minimal() +
#   geom_point(data=hosts, aes(x, y, col=Species), pch=16, size=1) +
#   scale_color_viridis_d(begin=0, end=0.8, option="inferno") +
#   theme(axis.text = element_blank(),
#         axis.title=element_blank())

# species
species = c("Lepus capensis", "Lepus microtis")

# specify species to model (index)
nnn = 2
spp = species[nnn]

### make spdf 
occ = hosts[ hosts$Species == spp, ]
sp::coordinates(occ) = ~x+y

# if lsax, crop using iucn polygon (extremely range restricted)
# if(spp == "Lepus saxatilis"){
#   lsax = readOGR("./data/iucn_rangemaps/lepus/lsax.shp")
# }

# thin occurrence points to limit spatial autocorrelation 
occ = thinOccurrences(occ, res=0.05)






# ========================= FITTING MODELS ================================

# boosted regression trees refresher:
# https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1365-2656.2008.01390.x
# combines large numbers of simple tree models adaptively to optimise predictive performance
# origins in ML but can also be considered an advanced form of regression

# decision trees: partition predictor space into rectangles and identify regions w/ most homogeneous response to predictor variables
# fit a constant to each region (for regression trees, mean response for obs in that region, assuming gaussian errors)
# regions are terminal nodes (leaves), split points; chosen to minimise prediction errors via fitting algorithm
# via trees interactions are automatically modelled: response to an input variable depends on resposes to other input variables higher up in tree

# boosting: rather than bagging/stacking/model averaging, boosting is forward and stagewise, each successive tree depends on last stage
# essentially gradually increases emphasis on observations modelled poorly by existing trees
# boosting minimises the loss function by iteratively adding new terms (trees) to best reduce the loss function
# can see how this is sensitive to overfitting: hence holdout approaches
# loss function == binomial deviance (i.e. difference in log-lik between current and saturated model)
# first tree is fitted to maximally reduce the loss function (i.e. sum of residuals)
# subsequent trees are fitted to the residuals, and the fitted values are added to the current logit(p) pointwise, i.e. incrementally reducing deviance
# "shrinkwrapping"
# fitted value re-estimated at each iteration to reflect contribution of newly added tree
# final fitted values are computed by the sum of all trees, multiplied by the learning rate (i.e. sum of trees model)

# n.b. bagging: using a random subsample to fit each fitted tree (i.e. introducing stochasticity), bag fraction is proportion selected on each iteration
# in comparison to random forest, fits to random subsample of x% of all training data (whereas random forest == sampling w/ replacemnet)
# learning rate = contribution of each tree to final fitted model
# tree complexity = degree to which interactions are fitted, from tc1 (single decision stump, 2 terminal nodes), tc2 = 2 way interactions, etc
# these two parameters together determine the nt (number of trees) required for optimal prediction (i.e. prediction that minimises holdout deviance)



# ------------ fit models with species-specific calibration: ensemble of 200 per species --------------

# checkerboard = 2 (4-fold, spatially structured)
if(spp %in% c("Lepus capensis", "Lepus microtis")){
  s = Sys.time()
  brt_ens = fitModelEnsemble(occ, species=as.vector(spp), modeltype="BRT", num_pseudo=50, mask_raster=geo_mask, k=NA, checkerboard=TRUE, checker_type = 2, learning_rate=0.0005, num_pseudo_points=nrow(occ), filepath="./output/model_outputs/climate_ensemble/fitted_models/")
  hist(brt_ens$auc, main="AUC scores, BRT models")
  e = Sys.time()
  e-s
}




# =================== extract and save diagnostic info and partial dependencies for ensemble ======================

# getEnsembleParams: extracts for each model, variable importance and partial dependence plots, saves and summarises
#' @param filepath directory where model ensemble is stored
 
getEnsembleParams = function(filepath){
  
  # extract partial dependency plots and variable importance
  getGBMVariables = function(file){
    
    # loads file and creates data frame of variable importance
    load(file)
    #print(file)
    resx = res$model$contributions
    row.names(resx) = c()
    names(resx) = c("Variable", "RelativeInfluence")
    resx$Species = res$species
    resx$model_id = res$unique_id
    
    # get response curves
    partials = data.frame()
    for(i in 1:nrow(resx)){
      pdx = plot.gbm(res$model, i.var=i, return.grid=TRUE)
      pdx$var = names(pdx)[1]
      names(pdx)[1] = "var_value"
      partials = rbind(partials, pdx)
    }
    partials$Species = res$species
    partials$model_id = res$unique_id
    
    return(list(var_imp = resx,
                partials = partials))  
  }
  
  print("Loading models and extracting data...")
  files = list.files(filepath, pattern=".R", full.names = TRUE)
  vps = lapply(files, getGBMVariables)
  
  print("Summarising and saving...")
  var_imps = do.call(rbind.data.frame, lapply(vps, FUN=function(x) x$var_imp))
  partials = do.call(rbind.data.frame, lapply(vps, FUN=function(x) x$partials))
  
  # summarise variable importance to order variables in plot
  vi_summary = var_imps %>%
    group_by(Variable) %>%
    dplyr::summarise(var = unique(Variable),
                     MeanRelInfluence = mean(RelativeInfluence),
                     LowerRelInfluence = quantile(RelativeInfluence, 0.05),
                     UpperRelInfluence = quantile(RelativeInfluence, 0.95))
  vi_summary = vi_summary[ order(vi_summary$MeanRelInfluence, decreasing = TRUE), ]
  vi_summary$var2 = factor(vi_summary$var, levels=vi_summary$var, ordered=TRUE)
  
  # save variable importance
  write.csv(var_imps, paste(filepath, "variableimportance.csv", sep=""), row.names=FALSE)
  
  # combine with partials
  partials = left_join(partials, vi_summary[ , c("var", "var2", "MeanRelInfluence")])
  write.csv(partials, paste(filepath, "variable_partialdependencies.csv", sep=""), row.names=FALSE)
  
  # create plots
  var_imps$Variable = factor(var_imps$Variable, levels=as.vector(vi_summary$var), ordered=TRUE)
  vi_plot = ggplot(var_imps, aes(Variable, RelativeInfluence)) + 
    ggforce::geom_sina(width=1, pch=21, fill="grey70", col="grey20", alpha=0.15) +
    geom_boxplot(width=0.4, fill="skyblue4", alpha=0.5, outlier.shape = NA) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, size=10),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=13),
          plot.title = element_text(size=14, hjust=0.5)) +
    ylab("Relative variable influence (%)") +
    ggtitle(var_imps$Species[1])

  # plot
  pd_plot = ggplot(data=partials, aes(var_value, y, group=model_id)) + 
    geom_line(alpha=0.1, col="skyblue4") +
    #stat_summary(aes(group=model_id), fun=mean, colour="black", geom="line") +
    facet_wrap(~var2, scales="free_x") + 
    theme_classic() +
    theme(plot.title = element_text(size=14, hjust=0.5),
          strip.background = element_blank()) +
    geom_hline(yintercept=0, lty=2) +
    xlab("Covariate value") +
    ylab("Fitted function") +
    ggtitle(partials$Species[1])
  
  # return
  return(list(vi_data = vi_summary,
              vi_plot = vi_plot,
              pd_data = partials,
              pd_plot = pd_plot))
}

# run for just fitted model (search for unique identifier)
foo = list.files("./output/model_outputs/climate_ensemble/fitted_models/", full.names=TRUE)
filepathx = paste(foo[ grep(brt_ens$ensemble_id[1], foo) ], "/", sep="")
model_diag = getEnsembleParams(filepathx)
model_diag$vi_plot
model_diag$pd_plot






# ================== geographical predictions from model ensemble ===================

# predictSpDistribution: uses BRT to predict geographical species distributions for each ensemble model
#' @param filepath directory where ensemble is stored
#' @param env_rasters raster stack of new env data to predict on
#' @param aggregation_factor if !is.na, aggregates env rasters (mean) by a desired factor to reduce computing time

predictSpDistribution = function(filepath, env_rasters, aggregation_factor=NA){
  
  # aggregate if specified
  if(is.numeric(aggregation_factor) & !is.na(aggregation_factor)){
    env_rasters = aggregate(env_rasters, fact=aggregation_factor, fun=mean)
  }
  
  # list models
  files = list.files(filepath, pattern=".R", full.names = TRUE)
  
  # ensemble identifier
  ens_id = strsplit(filepath, "_")[[1]][length(strsplit(filepath, "_")[[1]])]
  ens_id = strsplit(ens_id, "/")[[1]][1]
  
  # create save dir
  save_dir = paste(filepath, "geo_preds/", sep="")
  dir.create(save_dir)
  
  # predict per model
  predModel = function(x){
    
    # loads file and creates data frame of variable importance
    cat(paste(x, "...", sep=""))
    load(files[x])
    
    # predict 
    px = predict(env_rasters, res$model, n.trees=res$model$gbm.call$best.trees, type="response")
    
    # create filename
    ras_name = paste(res$species, ens_id, res$unique_id, "AUC", round(res$auc, 3), "THRESH", round(res$threshold, 4), sep="_")
    names(px) = ras_name
    
    # save and remove raster
    writeRaster(px, filename=paste(save_dir, ras_name, ".tif", sep=""), format="GTiff", overwrite=TRUE)
    rm(px); rm(res)
  }
  
  # run predictions
  s = Sys.time()
  print("Predicting rasters:")
  lapply(1:length(files), predModel)
  e = Sys.time()
  print(e-s)
}

# run for models
predictSpDistribution(filepath = filepathx, envx, aggregation_factor = 2)





# ======================= save and plot mean prediction with points overlaid ============================

createMeanPrediction = function(filepath, speciesx){
  
  # get rasters
  rr = stack(list.files(paste(filepath, "geo_preds", sep=""), full.names=TRUE))
  
  # ensemble ID
  ens_id = strsplit(filepath, "_")[[1]][length(strsplit(filepath, "_")[[1]])]
  ens_id = strsplit(ens_id, "/")[[1]][1]
  
  # create dir
  dir.create(paste(filepath, "mean_geo/", sep=""))
  
  # calculate mean
  rr_mean = mean(rr)
  writeRaster(rr_mean, file=paste(filepath, "mean_geo/", speciesx, "_meanpred.tif", sep=""), format="GTiff", overwrite=TRUE)
  return(rr_mean)
}

rrm = createMeanPrediction(filepathx, speciesx=spp)
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlGnBu"))(200)
plot(rrm, col=colscalex)


# # plot mean prediction
# predx = mean(ss)
# shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
# ras = as.data.frame(predx, xy=TRUE)
# p1 = ggplot() + 
#   geom_raster(data=ras, aes(x, y, fill=layer)) +
#   theme_classic() +
#   theme(legend.title = element_blank(),
#         plot.title = element_text(hjust=0.5, vjust=-3),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line = element_blank()) + 
#   #scale_fill_viridis(option = "inferno", na.value = "white") +
#   scale_fill_viridis(option = "viridis", na.value = "white") +
#   geom_sf(data=shpx, colour="black", fill=NA, size=0.15) +
#   geom_point(data=as.data.frame(occ), aes(x, y), pch=21, col="white", cex=1)
# ggsave(p1, file=paste("./output/figures/draft/", spp, "_pred_viz.png", sep=""), device="png", dpi=600, width=10, height=12, units="in")
