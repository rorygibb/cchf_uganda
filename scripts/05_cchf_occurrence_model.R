

# ====================== CCHF BRT DISTRIBUTION MODELS =============================


# fits, evaluates and predicts spatially from ensembles of boosted regression tree SDM for the distribution of CCHF exposures in people



# ------------- dependencies and housekeeping --------------

# root dir
#setwd("D:/ResearchProjects/202002_ucl_congocrimean/analysis/")
#setwd("C:/Users/roryj/Dropbox/CrimeanCongo_Uganda/")
setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")

library(dplyr)
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
library(embarcadero)
library(sf)
library(ENMeval)
library(exactextractr)
library(dismo)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)

# extent
af_ext = extent(c(-20, 60, -40, 37))

# colscale
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)

# mapping theme
maptheme = theme_classic() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_text(hjust=0.5, size=14),
        legend.title = element_text(size=10), 
        strip.background = element_blank(),
        strip.text = element_text(size=16))

# shapefiles
shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
shpy = readOGR("./data/project_shapefiles/africa_admin0.shp")




# ======================== MODELLING FUNCTIONS =============================


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

#' @param occ dataframe of occurrence points (presences)
#' @param num_pseudo specify number of pseudoabsence regimes
#' @param k specify k for kfold train-test if random
#' @param checkerboard specify if partitioning via checkerboard spatial approach
#' @param checker_type if checkerboard, eiher 1 (splits into 2 approx equal size checkerboard folds) or 2 (splits into 4 via nested checkerboards)
#' @param mask_raster raster to constrain pseudo points (e.g. to land)
#' @param num_pseudo_points specify number of pseudoabsence points per regime

generateDataset = function(occ, env_layers, num_pseudo=5, k=NA, checkerboard=FALSE, checker_type=2, mask_raster, num_pseudo_points=200){
  
  # for storing results
  result = list()
  
  # create n sets of pseudoabsences, combine with occurrence data, and partition into 
  for(i in 1:num_pseudo){
    
    # return report
    print(paste("Generating pseudoabsences regime ", i, sep=""))
    
    # pres data
    pres_data = occ
    
    # create random pseudoabsences and combine with presences (num points = num occ)
    pa = generatePseudos(pres_data, num_points=200, mask_raster=mask_raster, weighting="none", ext=extent(mask_raster))
    occx = data.frame(x = pres_data@coords[, 1], y = pres_data@coords[, 2], type = "presence")
    df = rbind(occx, pa)
    df$occurrence = ifelse(df$type == "presence", 1, 0)
    df$pseud_regime = i
    
    # extract environmental variables
    # ed = extractEnviroData(df$x, df$y, env_layers)
    # df = cbind(df, ed)
    
    # extract within a 10km buffer
    occ_sf = st_as_sf(df, coords = c("x", "y") )
    st_crs(occ_sf) = crs(env_layers)
    buf = st_buffer(occ_sf, 0.1) # 111 km = 1 decimal degree/ 25 km = 25/111
    ed = exactextractr::exact_extract(env_layers, buf, fun='mean')
    names(ed) = substr(names(ed), 6, 200)
    df = cbind(df, ed)
    
    # remove any containing NaN or all NAs
    anynan = function(x){ any(is.nan(x)) }
    contains_nan = apply(ed, 1, anynan)
    allna = function(x){ all(is.na(x)) }
    all_nas = apply(ed, 1, allna)
    df = df[ !(contains_nan | all_nas), ]
    
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
      if(checker_type == 1){ checker_points = ENMeval::get.checkerboard1(cbind(occ.pts$x, occ.pts$y), env_layers, cbind(abs.pts$x, abs.pts$y), aggregation.factor=2) }
      if(checker_type == 2){ checker_points = ENMeval::get.checkerboard2(cbind(occ.pts$x, occ.pts$y), env_layers, cbind(abs.pts$x, abs.pts$y), aggregation.factor=2) }
      occ.pts$k = checker_points$occs.grp
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

fitModelEnsemble = function(occ, env_layers, species, modeltype, num_pseudo=5, mask_raster, k=NA, checkerboard=FALSE, checker_type=2, learning_rate=0.005, num_pseudo_points=200, filepath="./output/model_outputs/climate_ensemble/fitted_models/"){
  
  # generate datasets with specified regime
  print("Generating dataset...")
  sdm_data = generateDataset(occ, env_layers=env_layers, num_pseudo=num_pseudo, k=k, checkerboard=checkerboard, checker_type=checker_type, num_pseudo_points=num_pseudo_points, mask_raster=mask_raster)
  print(paste(class(sdm_data), "length", length(sdm_data)), sep=" ")
  
  # create directory for saving ensemble outputs
  identifier = sample(1:10^6, 1)
  dir_name = paste(species, modeltype, Sys.Date(), identifier, sep="_")
  model_path = paste(filepath, dir_name, "/", sep="")
  dir.create(model_path)
  
  # save occurrences
  dir.create(paste(model_path, "occurrence_points/", sep=""))
  write.csv(as.data.frame(occ), paste(model_path, "occurrence_points/", species, "_occurrences.csv", sep=""), row.names=FALSE)
  
  # predictors
  predictors = names(env_layers)
  
  # ### fitGLM: fits, evaluates and thresholds logistic GLM species distribution model
  # fitGLM = function(x){
  #   
  #   # create unique identifier
  #   identifier = sample(1:10^8, 1)
  #   
  #   # data
  #   dfx = sdm_data[[x]]
  #   dfx_train = dfx[ dfx$train == 1, ] 
  #   dfx_test = dfx[ dfx$train == 0, ]
  #   
  #   # fit logistic GLM model to training subset
  #   mod = glm(occurrence ~ bio01_temp_annualmean + bio04_temp_seasonality + 
  #               bio08_temp_meanwettestquarter + bio09_temp_meandriestquarter +
  #               bio12_precip_annualtotal + bio15_precip_seasonality +
  #               bio16_precip_wettestquarter + bio17_precip_driestquarter,
  #             family=binomial(link="logit"),
  #             data=dfx_train)
  #   
  #   # evaluate and threshold
  #   e = evaluate(p = dfx_test[ dfx_test$occurrence == 1,], a = dfx_test[ dfx_test$occurrence == 0,], mod)
  #   invlogit = function(x) exp(x)/(1+exp(x))
  #   t = invlogit(threshold(e, 'spec_sens'))
  #   
  #   # return result list
  #   res = list(modeltype = "GLM_logistic", data = dfx, model = mod, auc = e@auc, threshold = t, unique_id = identifier)
  #   return(res)
  # }
  # 
  # ### fitGAM: fits, evaluates and thresholds logistic GAM species distribution model
  # fitGAM = function(x){
  #   
  #   # create unique identifier
  #   identifier = sample(1:10^8, 1)
  #   
  #   ## http://environmentalcomputing.net/intro-to-gams/
  #   dfx = sdm_data[[x]]
  #   dfx_train = dfx[ dfx$train == 1, ] 
  #   dfx_test = dfx[ dfx$train == 0, ]
  #   
  #   # fit logistic GAM model to training subset
  #   mod = gam(occurrence ~ s(bio01_temp_annualmean) + s(bio04_temp_seasonality) + 
  #               s(bio08_temp_meanwettestquarter) + s(bio09_temp_meandriestquarter) +
  #               s(bio12_precip_annualtotal) + s(bio15_precip_seasonality) +
  #               s(bio16_precip_wettestquarter) + s(bio17_precip_driestquarter),
  #             family=binomial(link="logit"),
  #             data=dfx_train)
  #   #gam.check(mod)
  #   
  #   # evaluate and threshold based on test set
  #   e = evaluate(p = dfx_test[ dfx_test$occurrence == 1,], a = dfx_test[ dfx_test$occurrence == 0,], mod)
  #   invlogit = function(x) exp(x)/(1+exp(x))
  #   t = invlogit(threshold(e, 'spec_sens'))
  #   
  #   # return result list
  #   res = list(modeltype = "GAM_logistic", data = dfx, model = mod, auc = e@auc, threshold = t, unique_id = identifier)
  #   return(res)
  # }
  
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
    mod = dismo::gbm.step(data = dfx_train, 
                          gbm.x = 2:ncol(dfx_train), 
                          gbm.y = 1, 
                          family = "bernoulli", 
                          tree.complexity=4, 
                          learning.rate=learning_rate, 
                          bag.fraction=0.75,
                          silent=TRUE)
    
    # evaluate and threshold based on test dataset
    e = dismo::evaluate(p = dfx_test[ dfx_test$y == 1,], a = dfx_test[ dfx_test$y == 0,], mod, n.trees=mod$gbm.call$best.trees)
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
  
  ### fitMaxent: fits, evaluates and thresholds species occurrence model using Maxent
  # fitMaxent = function(x){
  #   
  #   # create unique identifier
  #   identifier = sample(1:10^8, 1)
  #   
  #   # data
  #   print(x)
  #   dfx = sdm_data[[x]]
  #   dfx_train = dfx[ dfx$train == 1, ] 
  #   dfx_test = dfx[ dfx$train == 0, ]
  #   
  #   # fit maxent model to training subset
  #   mx = maxent(x = dfx_train[ , 6:13], p = dfx_train$occurrence)
  #   #plot(maxent_sdm)
  #   #response(maxent_sdm)
  #   
  #   # evaluate and threshold based on test dataset
  #   e = evaluate(p = dfx_test[ dfx_test$occurrence == 1,], a = dfx_test[ dfx_test$occurrence == 0,], mx)
  #   t = threshold(e, 'spec_sens')
  #   # maxent_prediction = predict(env_layers, mx, ext=extent(env_layers), progress='')
  #   # plot(maxent_prediction > t)
  #   
  #   # return result list
  #   res = list(modeltype = "Maxent", data = dfx, model = mx, auc = e@auc, threshold = t, unique_id = identifier)
  #   return(res)
  # }
  
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



# # ================== prepare environmental layers ==================
# #
# #
# # --- 1a. Hyalomma tick distribution ---
# 
# # maximum occurrence for 4 species
# r1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/mean_geo/Hyalomma rufipes_meanpred.tif")
# r2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/mean_geo/Hyalomma truncatum_meanpred.tif")
# r3 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma dromedarii_BRT_2020-06-24_111728_wolepus/mean_geo/Hyalomma dromedarii_meanpred.tif")
# r4 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma impeltatum_BRT_2020-07-01_686854_wlepus/mean_geo/Hyalomma impeltatum_meanpred.tif")
# hya = stack(r1, r2, r3, r4)
# hyalomma = calc(hya, fun = max)
# hyalomma = crop(hyalomma, af_ext)
# names(hyalomma) = "Hyalomma_suitability"
# 
# 
 
# --- 1b. Amblyomma and Rhipicephalus distribution ----

# # 
# a1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Amblyomma variegatum_BRT_2022-05-17_856179/mean_geo/Amblyomma variegatum_meanpred.tif")
# 
# a2.1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Rhipicephalus appendiculatus_BRT_2022-05-18_96167/mean_geo/Rhipicephalus appendiculatus_meanpred.tif")
# a2.2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Rhipicephalus (Boophilus) decoloratus_BRT_2022-05-20_327891/mean_geo/Rhipicephalus (Boophilus) decoloratus_meanpred.tif")
# a2 = stack(a2.1, a2.2);
# a2 = calc(a2, fun = max)
# 
# ticks = stack(a1, a2)
# names(ticks) = c("Amblyomma_suitability", "Rhipicephalus_suitability")
# ticks = crop(ticks, af_ext)


# # --- 2. Lepus hare distribution ---
# 
# # maximum occurrence for 2 species
# l1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Lepus capensis_BRT_2020-06-29_236594/mean_geo/Lepus capensis_meanpred.tif")
# l2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Lepus microtis_BRT_2020-07-03_70644/mean_geo/Lepus microtis_meanpred.tif")
# lepus = stack(l1, l2)
# lepus = calc(lepus, fun=max)
# lepus = crop(lepus, af_ext)
# names(lepus) = "Lepus_suitability"
# 
# 
# # # --- 3. Cropland area (ESA-CCI) ---
# #
# # # --------------- prepare ESA land cover data -------------
# #
# # # lu1 = raster("./data/landcover_landuse/esa_cci/esa_africa_2015.tif")
# # # lu1 = crop(lu1, extent(shpy))
# # #
# # # aggregateESA = function(esa_ras, fact, classes){
# # #   propx = function(x, ...){ sum(x %in% classes) / sum(!is.na(x) & x != 210) }
# # #   rx = aggregate(esa_ras, fact=fact, propx)
# # #   rx
# # # }
# # #
# # # cropland = aggregateESA(lu1, 15, classes=10:20)
# # # cropmosaic = aggregateESA(lu1, 15, classes=10:40)
# # # #forest = aggregateESA(lu1, 15, classes=c(50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 160, 170))
# # # #shrub_grass = aggregateESA(lu1, 15, classes=c(130, 120:122))
# # #
# # # # save
# # # writeRaster(cropland, "./output/esa_cci_processed/Cropland_strict_Africa_2015.tif", format="GTiff")
# # # writeRaster(cropmosaic, "./output/esa_cci_processed/Cropland_Africa_2015.tif", format="GTiff")
# 
# # cropland
# cropland = raster("./output/esa_cci_processed/Cropland_Africa_2015.tif")
# cropland = crop(cropland, af_ext)
# cropland = aggregate(cropland, fact=2, fun=mean)
# cropland = resample(cropland, lepus, fun="ngb")
# names(cropland) = "Cropland_proportion"
# 
# 
# # --- 4. Livestock ---
# 
# # cattle (non downscaled)
# cattle = raster("./data/griddedlivestock_v3/cattle_2010/6_Ct_2010_Aw.tif")
# sheep = raster("./data/griddedlivestock_v3/sheep_2010/6_Sh_2010_Aw.tif")
# goats = raster("./data/griddedlivestock_v3/goats_2010/6_Gt_2010_Aw.tif")
# livestock = stack(cattle, sheep, goats)
# livestock = crop(livestock, af_ext)
# livestock = resample(livestock, lepus, fun="ngb")
# log_livestock_dens = log((livestock /  raster::area(livestock[[1]]))+1)
# names(log_livestock_dens) = c("Cattle_logDensity", "Sheep_logDensity", "Goats_logDensity")
# ruminant_dens = sum(livestock) / raster::area(livestock[[1]])
# log_ruminant_dens = log(ruminant_dens+1)
# names(log_ruminant_dens) = c("Ruminants_logDensity")
# 
# 
# # --- 5. Human density ---
# 
# human = raster("./data/ciesin_griddedpopworld_v4/gpw_v4_population_count_rev11_2015_2pt5_min.tif")
# human = crop(human, af_ext)
# human = aggregate(human, fact=2, fun=mean)
# human_dens = log((human / raster::area(human))+1)
# names(human_dens) = "HumanPop_logDensity"
# 
# 
# # # --- 6. Climate covariates --------
# 
# # climatic predictors (already cropped to required extent)
# env = raster::stack(list.files("./data/chelsa_climatologies/africa_rescaled_aggregated/", pattern=".tif", full.names=T))
# env = processCHELSA(env)
# 
# # NDVI (Copernicus 1km original; derived from long-term statistics data 1992-2017)
# nd = stack(list.files("./output/ndvi_copernicus/africa_summary/", full.names = TRUE, pattern=".tif"))
# nd = resample(nd, env, "bilinear")
# names(nd) = c("NDVI_AnnualMean", "NDVI_AnnualSD", "NDVI_AnnualMin")
# env = stack(env, nd)
# 
# # subset to covariates of interest
# predictors = c("bio01_temp_annualmean",
#                "bio05_temp_maxwarmestmonth",
#                "bio06_temp_mincoldestmonth",
#                "bio12_precip_annualtotal",
#                "bio13_precip_wettestmonth",
#                "bio14_precip_driestmonth",
#                "NDVI_AnnualMean",
#                "NDVI_AnnualSD",
#                "NDVI_AnnualMin")
# envx = env[[ which(names(env) %in% predictors)]]
# rm(env)
# 
# # crop to same extent and calculate mean
# envx = crop(envx, af_ext)
# envx = aggregate(envx, fact=2, fun="mean")
# 
# # combine all into a single stack
# cchf_predictors = stack(hyalomma, ticks, lepus, cropland, human_dens, log_ruminant_dens, log_livestock_dens, envx)
# plot(cchf_predictors[[16]], col=colscalex)
# 
# # save
# writeRaster(cchf_predictors, file=paste("./output/model_outputs/cchf_brt/predictors/", names(cchf_predictors), ".tif", sep=""), format="GTiff", bylayer=TRUE, overwrite=TRUE)



# ==================== cchf occurrence data =====================

# # messina 2015
# a1 = read.csv("./data/cchf_occurrences/messsina_db_2013/CCHF_1953_2012_Messina.csv", stringsAsFactors = FALSE) %>%
#   dplyr::rename("y"=LATITUDE, "x"=LONGITUDE)
# a1 = a1[ a1$REGION == "Africa" & !a1$COUNTRY %in% c("United Arab Emirates", "Oman", "Madagascar"), ]
# 
# a1$Decade = NA
# a1$Decade[ a1$YEAR > 1950 & a1$YEAR < 1980 ] = "1958-1979"
# a1$Decade[ a1$YEAR >= 1980 & a1$YEAR < 2000 ] = "1980-1999"
# a1$Decade[ a1$YEAR >= 2000 ] = "2000-2012"
# 
# a2 = a1[ , c("x", "y", "Decade")]
# a2$Source = "Messina 2015"
# a2$Type = "Human"
# 
# # RG geolocated from surveillance reports, pubmed, exclude Uganda becuase Balinandi et al 2022 provides specific locations
# a3 = read.csv("./data/cchf_occurrences/rg_db_/cchf_location_datasources_georef.csv", stringsAsFactors = FALSE) %>%
#   dplyr::filter(Country != "Uganda") %>%
#   dplyr::filter(Sero == FALSE | Sero_Type == "Diagnostic") %>%
#   dplyr::filter(Type == "Human") %>%
#   dplyr::filter(!is.na(Longitude))
# a4 = a3[ , c("Latitude", "Longitude", "Type")]
# names(a4) = c("y", "x", "Type")
# a4$Source = "RG 2020"
# a4$Decade = "2013-2019"
# 
# # Balinandi et al 2022
# a5 = read.csv("./data/cchf_occurrences/rg_db_/balinandi_et_al_2022/balinandietal_2022.csv") %>%
#   dplyr::rename("x"=Longitude, "y"=Latitude, "Decade"=Years) %>%
#   dplyr::select(-Country) %>%
#   dplyr::mutate(Type = "Human")
# 
# # combine
# cchf = rbind(a2, a4) %>% rbind(a5)
# cchf = cchf[ !is.na(cchf$y), ]
# write.csv(cchf, "./output/model_outputs/cchf_brt/occurrence_data/CCHF_Occurrences.csv", row.names=FALSE)
# 


# =================== read in data for models ===================

# occurrence data
occ = read.csv("./output/model_outputs/cchf_brt/occurrence_data/CCHF_Occurrences.csv", stringsAsFactors = FALSE)

# visualise
# px = ggplot() + 
#   geom_sf(data=shpx, fill=NA, col="grey30") + 
#   geom_point(data=occ %>% dplyr::rename("Epoch"=Decade), pch=21, col="black", aes(x, y, fill=Epoch), size=1.75, alpha=0.8) + 
#   maptheme +
#   theme(legend.text = element_text(size=12), legend.title = element_text(size=13))
# ggsave(px, file="./output/figures/paper/SuppFigure_CCHFdistribution.png", device="png", units="in", width=9, height=7, dpi=600)

# make spdf 
sp::coordinates(occ) = ~x+y

# thin occurrence points to limit spatial autocorrelation 
occ = thinOccurrences(occ)
occ$Species_current = "CCHF"
occ = occ[ occ$Type == "Human", ]

# read in environmental predictors
ff = list.files("./output/model_outputs/cchf_brt/predictors/", pattern=".tif", full.names=TRUE)
envx = stack(ff[ -grep("Amblyomma|Rhipi", ff) ])
envx = envx[[ -grep("bio01|bio12|bio13|bio14|Human|NDVI_AnnualMin|Ruminant|Sheep|Goat", names(envx)) ]] 
envy = stack(ff[ grep("Amblyomma|Rhipi", ff) ])
envy = crop(envy, envx)
envx = stack(envx, envy)
plot(envx, col=colscalex)

# mask for pseudos
geo_mask = raster::raster("./data/africa_maskraster/africa_mask.tif")

# pseudoabs to look at collinearity among vars: remove sheep and goat density
pp = generatePseudos(occ = occ, mask_raster = geo_mask, num_points = 1000)
coordinates(pp) = ~x+y
xx = raster::extract(envx, pp)
cc = cor(xx[ complete.cases(xx), ])
corrplot::corrplot(cc)



# =================== fit model ensemble ================

brt_ens = fitModelEnsemble(occ, env_layers=envx, species="cchf_10km_exp2", modeltype="BRT", num_pseudo=50, mask_raster=geo_mask, k=NA, checkerboard=TRUE, checker_type = 2, learning_rate=0.0007, num_pseudo_points=500, filepath="./output/model_outputs/cchf_brt/fitted_models/")
hist(brt_ens$auc, 20)




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
foo = list.files("./output/model_outputs/cchf_brt/fitted_models/", full.names=TRUE)
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
predictSpDistribution(filepath = filepathx, envx, aggregation_factor = NA)




# ======================= save and plot mean prediction with points overlaid ============================

createMeanPrediction = function(filepath, speciesx, stat="mean"){
  
  # get rasters
  ff = list.files(paste(filepath, "geo_preds", sep=""), full.names=TRUE)
  ff = ff[ -grep(".xml", ff) ]
  rr = stack(ff)
  
  # ensemble ID
  ens_id = strsplit(filepath, "_")[[1]][length(strsplit(filepath, "_")[[1]])]
  ens_id = strsplit(ens_id, "/")[[1]][1]
  
  # create dir
  dir.create(paste(filepath, "mean_geo/", sep=""))
  
  # calculate mean
  if(stat == "mean"){ rr_mean = mean(rr) }
  if(stat == "median"){ 
    rr_mean = raster::calc(rr, fun=median, na.rm=T)
    # dfx = as.data.frame(rr)
    # dfx$median = apply(dfx, 1, median, na.rm=T)
    # rr_mean = raster(rr)
    # values(rr_mean) = dfx$median
  }
  if(stat == "sd"){ 
    rr_mean = raster::calc(rr, fun=sd, na.rm=T)
    # dfx = as.data.frame(rr)
    # dfx$median = apply(dfx, 1, median, na.rm=T)
    # rr_mean = raster(rr)
    # values(rr_mean) = dfx$median
  }
  writeRaster(rr_mean, file=paste(filepath, "mean_geo/", speciesx, "_", stat, "pred.tif", sep=""), format="GTiff")
  return(rr_mean)
}

rrm = createMeanPrediction(filepathx, speciesx="CCHF", stat="median")
rrm = createMeanPrediction(filepathx, speciesx="CCHF", stat="mean")
rrs = createMeanPrediction(filepathx, speciesx="CCHF", stat="sd")
plot(rrm, col=viridisLite::mako(200))
#plot(rrs/sqrt(rrm), col=viridis::magma(200))





# ----------------------- create summary figures ----------------------------


# ------------------ 1. create composite plot of env predictors and response curves -----------------------

# ---------------- functions to summarise -----------------

getVariableImportance = function(filepath){
  read.csv(paste(filepath, "variableimportance.csv", sep=""), stringsAsFactors = FALSE)
}

getPartialDependencies = function(filepath){
  read.csv(paste(filepath, "variable_partialdependencies.csv", sep=""), stringsAsFactors = FALSE)
}

# function to create figures for each model
createModelFigure = function(filepath){
  
  # get data
  vi = getVariableImportance(filepath)
  pd = getPartialDependencies(filepath)
  
  # function to rename vars (vector)
  renameVars = function(x){
    x[ x == "bio12_precip_annualtotal" ] = "Precip annual mean (bio12)"
    x[ x == "bio01_temp_annualmean" ] = "Temp annual mean (bio01)"
    x[ x == "bio04_temp_seasonality" ] = "Temp seasonality (bio04)"
    x[ x == "bio15_precip_seasonality" ] = "Precip seasonality (bio04)"
    x[ x == "bio05_temp_maxwarmestmonth" ] = "Temp max warmest (bio05)"
    x[ x == "bio06_temp_mincoldestmonth" ] = "Temp min coldest (bio06)"
    x[ x == "NDVI_AnnualMean" ] = "NDVI annual mean"
    x[ x == "NDVI_AnnualSD" ] = "NDVI seasonality (sd)"
    x[ x == "Lepus_Suitability" ] = "Lepus spp. suitability"
    x[ x == "bio02_temp_meandiurnalrange" ] = "Temp diurnal range (bio02)"
    x[ x == "bio13_precip_wettestmonth" ] = "Precip mean wettest (bio13)"
    x[ x == "bio14_precip_driestmonth" ] = "Precip mean driest (bio14)"
    x[ x == "Sheep_logDensity" ] = "Sheep density (log)"
    x[ x == "Cattle_logDensity" ] = "Cattle density (log)"
    x[ x == "Cropland_proportion" ] = "Agricultural land proportion"
    x[ x == "Lepus_suitability" ] = "Lepus spp. suitability"
    x[ x == "Hyalomma_suitability" ] = "Hyalomma spp. suitability"
    x[ x == "Rhipicephalus_suitability" ] = "Rhipicephalus spp. suitability"
    x[ x == "Amblyomma_suitability" ] = "Amblyomma variegatum suitability"
    x[ x == "Goats_logDensity" ] = "Goats density (log)"
    
    x
  }
  vi$Variable = renameVars(vi$Variable)
  pd$var = renameVars(pd$var)
  
  # plot var importance
  vi_summary = vi %>%
    group_by(Variable) %>%
    dplyr::summarise(var = unique(Variable),
                     MedianRelInfluence = median(RelativeInfluence),
                     MeanRelInfluence = mean(RelativeInfluence),
                     LowerRelInfluence = quantile(RelativeInfluence, 0.05),
                     UpperRelInfluence = quantile(RelativeInfluence, 0.95))
  vi_summary = vi_summary[ order(vi_summary$MedianRelInfluence, decreasing = TRUE), ]
  vi$Variable = factor(vi$Variable, levels=as.vector(unique(vi_summary$var)), ordered=TRUE)
  vi_plot = ggplot(vi, aes(Variable, RelativeInfluence)) + 
    ggforce::geom_sina(width=1, pch=21, fill="grey70", col="grey20", alpha=0.1) +
    geom_boxplot(width=0.4, fill="skyblue4", alpha=0.5, outlier.shape = NA) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, size=11),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=11),
          axis.title.y = element_text(size=11),
          plot.title = element_text(size=15, hjust=0.5)) +
    ylab("Relative variable influence (%)") 
  
  # plot top 8 importance variables
  pd = pd[ pd$var %in% vi_summary$var[1:8], ]
  pd$var = factor(pd$var, levels=vi_summary$var[1:8], ordered=TRUE)
  pd_plot = ggplot(data=pd, aes(var_value, y, group=model_id)) + 
    geom_line(alpha=0.1, col="skyblue4") +
    #stat_summary(aes(group=model_id), fun=mean, colour="black", geom="line") +
    facet_wrap(~var, scales="free_x", nrow=2) + 
    theme_classic() +
    theme(plot.title = element_text(size=14, hjust=0.5),
          axis.text = element_text(size=11),
          axis.title = element_text(size=11),
          strip.text = element_text(size=11),
          strip.background = element_blank()) +
    geom_hline(yintercept=0, lty=2) +
    xlab("Covariate value") +
    ylab("Fitted function")
  
  #plot_composite = gridExtra::grid.arrange(pd_plot, vi_plot, nrow=2, heights=c(1, 1))
  plot_composite = gridExtra::grid.arrange(grobs=list(pd_plot, vi_plot), heights=c(1, 1), nrow=2)
  
  dir.create(paste(filepath, "summaryfigures", sep=""))
  ggsave(plot_composite, filename=paste(filepath, "summaryfigures/ensemble_model_summary.png", sep=""), device="png", units="in", dpi=300, width=10, height=9, scale=0.85)
  
  # combine
  return(plot_composite)
}

# save
createModelFigure(filepathx)



# -------------------- create -------------------

createModelMap = function(filepath){
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # read raster
  geo_pred = raster(paste(filepath, "mean_geo/CCHF_medianpred.tif", sep=""))
  
  # save
  ras = as.data.frame(geo_pred, xy=TRUE)
  names(ras)[3] = "val"
  colscalex = rev(viridisLite::mako(200))
  geo_plot = ggplot() +
    geom_raster(data=ras, aes(x, y, fill=val)) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c(0.25, 0.25),
          plot.title = element_text(hjust=0.5, size=18, vjust=-3),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(color="black", fill=NA)) +
    #scale_fill_viridis(option = "magma", na.value = "white") +
    scale_fill_gradientn(colors=colscalex, na.value = "white") +
    geom_sf(data=shpx, colour="grey25", fill=NA, size=0.25) +
    geom_point(data=occ, aes(x, y), pch=16, size=0.6, col="red", alpha=0.75)
  
  ggsave(geo_plot, filename=paste(filepath, "summaryfigures/riskmap.png", sep=""), device="png", units="in", dpi=600, width=9, height=9, scale=0.9)
}
createModelMap(filepathx)

# map of uncertainty
createUncertMap = function(filepath){
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # relative uncertainty: sd / sqrt(mean), because uncert expected to be higher when mean is higher
  geo_pred1 = raster(paste(filepath, "mean_geo/CCHF_sdpred.tif", sep=""))
  geo_pred2 = raster(paste(filepath, "mean_geo/CCHF_meanpred.tif", sep=""))
  ru_pred = geo_pred1 / sqrt(geo_pred2)
  
  # save
  ras = as.data.frame(ru_pred, xy=TRUE)
  names(ras)[3] = "val"
  colscalex = rev(viridisLite::mako(200))
  geo_plot = ggplot() +
    geom_raster(data=ras, aes(x, y, fill=val)) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = c(0.25, 0.25),
          plot.title = element_text(hjust=0.5, size=18, vjust=-3),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(color="black", fill=NA)) +
    #scale_fill_viridis(option = "magma", na.value = "white") +
    scale_fill_gradientn(colors=colscalex, na.value = "white") +
    geom_sf(data=shpx, colour="grey25", fill=NA, size=0.25) +
    geom_point(data=occ, aes(x, y), pch=16, size=0.6, col="red", alpha=0.75)
  
  ggsave(geo_plot, filename=paste(filepath, "summaryfigures/riskmap_sd.png", sep=""), device="png", units="in", dpi=600, width=9, height=9, scale=0.9)
}
createUncertMap(filepathx)



createUgandaMap = function(filepath){
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # read raster
  geo_pred = raster(paste(filepath, "mean_geo/CCHF_medianpred.tif", sep=""))
  
  colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)
  
  # plot subset for uganda
  shp_uga = st_read("./data/project_shapefiles/uga_admin1.shp")
  ux = extent(c(xmin=29, xmax=36, ymin=-2, ymax=4.5))
  cc2 = occ; coordinates(cc2) = ~x+y
  cc2 = as.data.frame(crop(cc2, ux))
  
  # crop hazard raster
  gp_uga = crop(geo_pred, ux)
  gp_uga = as.data.frame(gp_uga, xy=TRUE)
  
  uga_plot1 = ggplot() + 
    theme_classic() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust=0.5, vjust=-3),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + 
    geom_raster(data=gp_uga, aes(x, y, fill=CCHF_medianpred)) +
    geom_sf(data=shp_uga, colour="grey60", fill=NA, size=0.215) + 
    geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.5) +  
    geom_point(data=cc2, aes(x, y), col="red", fill="coral", size=2, pch=21, alpha=0.8) +
    scale_fill_gradientn(colors=colscalex, na.value = "white")
  
  ggsave(uga_plot1, filename=paste(filepath, "summaryfigures/uga_hyalommamap.png", sep=""), device="png", units="in", dpi=600, width=6, height=6, scale=0.9)
}

createUgandaMap(filepathx)
