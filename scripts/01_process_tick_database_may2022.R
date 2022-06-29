


# ------------- Reads and processes African tick occurrence dataset (Graeme Cumming) ----------------


setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")
#setwd("C:/Users/roryj/Dropbox/CrimeanCongo_Uganda/")

# ------------- dependencies --------------

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
library(sf)
library(ENMeval)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)


# ------------ Cumming tick database for Africa ---------------

# Reliability: index of the kind of data source and indication of how precise the coordinates may be
# A: coordinates given by author
# G: coordinates from a gazetteer
# M: coordinates from a map
# 1-3: 1 indicates high reliability, 3 indicates lower reliability
# exclude 3s; could also exlude 2s if not many,
tt = read.delim("./data/africanticks_cumming1998/Main_tick_database(.txt)/Main tick database (final).txt")

# get only records with 1/2 reliability
tt = tt[ grep("1|2", tt$Reliability), ]

# rename lat lon
names(tt)[5:6] = c("y", "x")




# --------- split apart year field into relevant dates -----------

# strsplit b comma and space
yyy = strsplit(as.vector(tt$year), " |,")

# write function to parse and extract dates
parse_date_field = function(x){
  
  # report
  if(x %in% c(1, seq(10, 10^6, by=10))){ print(x) }
  
  # item
  xx = yyy[[x]]
  
  # data frame for results
  result = data.frame(rec_id = x)
  
  # parse for before and remove "before" dates to avoid confusion with solo
  if("before" %in% xx){ 
    result$before_year = max(as.numeric(xx[ which(xx == "before")+1 ]))
    xx = xx[ -which(xx == "before")+1 ]
  } else{
    result$before_year = NA
  }
  
  # parse for year ranges
  if(any(grepl("-", xx))){
    year_ranges = xx[ grepl("-", xx) ]
    year_ranges = unlist(strsplit(year_ranges, "-"))
    result$yearrange_earliestyear = min(as.numeric(year_ranges))
    result$yearrange_latestyear = max(as.numeric(year_ranges))
  } else{
    result$yearrange_earliestyear = NA
    result$yearrange_latestyear = NA
  }
  
  # parse for numerics
  if(any(!is.na(as.numeric(xx)))){
    solo_years = as.numeric(xx[ !is.na(as.numeric(xx)) ])
    result$soloyear_earliest = min(solo_years)
    result$soloyear_latest = max(solo_years)
  } else{
    result$soloyear_earliest = NA
    result$soloyear_latest = NA
  }
  
  # get estimated decade of latest record
  years = as.numeric(as.vector(result[ 1, c(2, 3, 5, 6) ]))
  years = years[ !is.na(years) ]
  decade_est = max(years)
  decade_est = paste(substr(decade_est, 1, 3), "0", sep="")
  result$decade_est = as.numeric(decade_est)
  
  # save
  return(result)
}

# run example
date_parse = do.call(rbind.data.frame, lapply(1:length(yyy), parse_date_field))

# combine with tick data
tt = cbind(tt, date_parse[ , 2:ncol(date_parse)])

# subset to records from 1930 onwards
tt = tt[ tt$decade_est >= 1930, ]

# save
#write.csv(tt, "./data/africanticks_cumming1998/Main_tick_database(.txt)/RG_YearProcessed_23062020.csv", row.names = FALSE)
tt = read.csv("./data/africanticks_cumming1998/Main_tick_database(.txt)/RG_YearProcessed_23062020.csv")



# ------------------ estimate density of tick sampling across all of Africa ------------------------

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


### pointsDensityRaster: create raster of kernel density of sampled points across a grid
# using kernel density estimation from spatstat assuming gaussian kernel with specified bandwidth (sigma)
# additional arguments allow for specifying geographical extent of density estimation (e.g. across larger area than points)

#' @param x vector of points on x dimension (longitude)
#' @param y vector of points on y dimension (latitude)
#' @param xrange vector specifying x extent of kernel density raster (xmin, xmax; defaults to extent of points
#' @param yrange vector specifying y extent of kernel density raster (ymin, ymax; defaults to extent of points)
#' @param dimx number of grid cells along x dimension
#' @param dimy number of grid cells along y dimension
#' @param sigma bandwidth for kernel density estimation (see spatstat::density.ppp, default 2.5)
#' @param mask_raster raster of study area used to mask density raster

pointsDensityRaster = function(x, y, xrange, yrange, dimx, dimy, sigma=2.5, mask_raster=NULL){
  
  # check arguments and set defaults
  if(!is.null(mask_raster) & class(mask_raster) != "RasterLayer"){ return("Error: mask_raster must be of class RasterLayer") }
  if(length(x) != length(y)){ return("Error: x and y must be the same length") }
  if(length(is.na(x)) != length(is.na(y))){ return("Error: different number of NA vals in x and y") }
  #if(is.null(dimx) | is.null(dimy)){ }
  if(is.null(xrange) | is.null(yrange)){ xrange = c(min(x), max(x)); yrange = c(min(y), max(y))}
  
  # specify geographical extent (window) for density estimation 
  window = spatstat::owin(xrange=xrange, yrange=yrange)
  
  # create spatial point pattern, and calculate density, specifying ncells along x and y axes, and specifying bandwidth
  print("Creating point pattern...")
  pts = spatstat::ppp(x, y, window=window)
  dens = spatstat::density.ppp(pts, dimyx=c(dimy, dimx), sigma=sigma)
  densr = raster::raster(dens)
  
  # if mask raster is provided, resample to same resolution, and mask raster
  # currently this part is more computationally intensive: resampling across large areas is time-consuming
  # there may be a better way to do this
  print("Masking raster...")
  if(!is.null(mask_raster)){
    projection(densr) = projection(mask_raster)
    densr = resample(densr, mask_raster, 'ngb')
    densr = mask(densr, mask_raster)
  }
  return(densr)
}


# raster for masking
geo_mask = raster::raster("./data/africa_maskraster/africa_mask.tif")

# thin occurrence points so remove duplicates falling within same cell (so looking at sampling over space)
tt2 = tt[ !is.na(tt$x), ]
coordinates(tt2) = ~ x+y
tt2 = as.data.frame(thinOccurrences(tt2, res=0.05))

# calculate raster
dens_ras = pointsDensityRaster(x = tt2$x,
                               y = tt2$y,
                               xrange = c(xmin(geo_mask), xmax(geo_mask)),
                               yrange = c(ymin(geo_mask), ymax(geo_mask)),
                               dimx = dim(geo_mask)[1],
                               dimy = dim(geo_mask)[2],
                               sigma = 1,
                               mask_raster = geo_mask)
plot(dens_ras, col=viridis::viridis(200))
#writeRaster(dens_ras, file="./data/africanticks_cumming1998/AllTicks_SamplingDensity_sigma1.tif", format="GTiff")
dens_ras = raster("./data/africanticks_cumming1998/AllTicks_SamplingDensity_sigma1.tif")



# --------------------- subset to Hyalomma and Rhipicephalus spp only ---------------------

# tick taxonomy as used here

# Hyalomma aegyptium - Hyalomma aegyptium aegyptium; 
# Hyalomma marginatum - Hyalomma marginatum marginatum
# Hyalomma marginatum rufipes - Hyalomma aegyptium var. impressum; 
# Hyalomma marginatum turanicum - Hyalomma glabrum; Hyalomma turanicum
# Hyalomma rufipes - Hyalomma impressum rufipes; Hyalomma rufipes rufipes
# Hyalomma truncatum - Hyalomma impressum transiens; Hyalomma nitidum; Hyalomma transiens

# New taxonomy:

# Hyalomma marginatum - Hyalomma marginatum
# Hyalomma rufipes - Hyalomma rufipes; Hyalomma marginatum rufipes, Hyalomma impressum
# Hyalomma turanicum - Hyalomma marginatum turanicum
# Hyalomma truncatum - Hyalomma truncatum

# also Boophilus (Rhipicephalus) decoloratus


# --------- subset -----------

# updated taxonomic designations for each relevant species
new_taxon = data.frame(
  Species_current = c("Hyalomma marginatum", 
                      "Hyalomma rufipes",
                      "Hyalomma rufipes", 
                      "Hyalomma rufipes", 
                      "Hyalomma turanicum", 
                      "Hyalomma truncatum",
                      "Hyalomma dromedarii",
                      "Hyalomma aegyptium",
                      "Hyalomma impeltatum",
                      "Boophilus decoloratus",
                      "Rhipicephalus appendiculatus",
                      "Amblyomma variegatum"),
  Species = c("Hyalomma marginatum", 
              "Hyalomma rufipes", 
              "Hyalomma marginatum rufipes",
              "Hyalomma impressum", 
              "Hyalomma marginatum turanicum",
              "Hyalomma truncatum",
              "Hyalomma dromedarii",
              "Hyalomma aegyptium",
              "Hyalomma impeltatum",
              "Boophilus decoloratus", 
              "Rhipicephalus appendiculatus",
              "Amblyomma variegatum")
)


# subset data to hyalomma
tth = tt[ grep("Hyalomma|Boophilus decoloratus|Rhipicephalus appendiculatus|Amblyomma variegatum", tt$Species), ]
tth = dplyr::left_join(tth, new_taxon)
tth = tth[ !is.na(tth$Species_current), ]

# extract points density metric: probability of densest points being excluded is either 0.75 or 0.9
coordinates(tth) = ~x+y
tth$PointsDensity_90 = raster::extract(dens_ras / (cellStats(dens_ras, max)*(1/0.9)), tth, method="simple")
tth$PointsDensity_75 = raster::extract(dens_ras / (cellStats(dens_ras, max)*(1/0.75)), tth, method="simple")

# save
write.csv(as.data.frame(tth), "./data/africanticks_cumming1998/Main_tick_database(.txt)/RG_Hyalommasubset_23062020.csv", row.names = FALSE)



# # --------- view -----------

# shp0 = st_read("./data/project_shapefiles/africa_admin0.shp")
# px = ggplot() +
#   geom_sf(data=shp0, colour="grey30", fill="grey95", size=0.25) +
#   theme_minimal() +
#   geom_point(data=as.data.frame(tt[tt$Species == "Boophilus decoloratus", ]), aes(x, y, col=Species), pch=16, size=0.8, alpha=0.6, col="red") +
#   #scale_color_viridis_d(begin=0, end=0.8, option="inferno") +
#   theme(axis.text = element_blank(),
#         axis.title=element_blank())
# px = ggplot() +
#   geom_sf(data=shp0, colour="grey30", fill="grey95", size=0.25) +
#   theme_minimal() +
#   geom_point(data=as.data.frame(tt[tt$Species == "Amblyomma variegatum", ]), aes(x, y, col=Species), pch=16, size=0.8, alpha=0.6, col="red") +
#   #scale_color_viridis_d(begin=0, end=0.8, option="inferno") +
#   theme(axis.text = element_blank(),
#         axis.title=element_blank())



# -------------- data cleaning to exclude records that are likely to be erroneous or very imprecise coordinates --------------------

# data frame
ticks = as.data.frame(tth)


# 1. H rufipes: remove northbound outliers (likely to be h marginatum rufipes which is a subsp. of h marginatum)
ticks = ticks[ -which(ticks$Species_current == "Hyalomma rufipes" & ticks$y > 31.5), ]

# 2. H dromedarii: remove southbound outlier (may be an introduction)
ticks = ticks[ -which(ticks$Species_current == "Hyalomma dromedarii" & ticks$y < -16), ]

# 3. Rhipicephalus appendiculatus

#As per Walker et al Western most distribution of Rhi appen controversial but reliably identified from Western DRC and Western Zambia, 
#more Westernly may be misindentification e.g. nitens, punctus, zambesensis, dutonni
# Most Northern South Sudan 
# Walker Ticks of Domestic Animals in Africa: a Guide to Identification of Species 2014 http://www.alanrwalker.com/assets/PDF/tickguide-africa.pdf
# Walker et al The Genus Rhipicephalus (Acari, Ixodidae) A Guide to the Brown Ticks of the World 2011) 
# https://www-cambridge-org.libproxy.ucl.ac.uk/core/books/genus-rhipicephalus-acari-ixodidae/accounts-of-individual-species-occurring-in-the-afrotropical-region-pages-278-to-490/611578D1DFD0C426C94FFBC02664AD4F
# Therefore remove records to West of 20 East and North of 10N as unclear if definitely Rhi Appen and may be misindetification. 
ticks = ticks[ -which(ticks$Species_current == "Rhipicephalus appendiculatus" & ticks$x < 20), ]
ticks = ticks[ -which(ticks$Species_current == "Rhipicephalus appendiculatus" & ticks$y > 10), ]

# 4. Amblyomma variegatum

# Limit of distrubition unclear and overlap with other Amb species, likely limited in S distribution in S Africa by hebraem and overlap with pompsum in Angola area. 
# Overlap with lepidum in NE Africa
# Walker et al do not include occurrences in S Africa and W Angola, others have.

# https://www.researchgate.net/publication/43346733_Use_of_geographic_information_systems_to_identify_areas_at_risk_of_introducing_Amblyomma_variegatum_and_A_hebraeum_to_Italy
# Walker 1987 https://pubmed.ncbi.nlm.nih.gov/3329325/
# Walker Ticks of Domestic Animals in Africa: a Guide to Identification of Species 2014 http://www.alanrwalker.com/assets/PDF/tickguide-africa.pdf
# https://www.researchgate.net/publication/6458868_Using_invaded_range_data_to_model_the_climate_suitability_for_Amblyomma_variegatum_Acari_Ixodidae_in_the_New_World
# In GC dataset only a handful of occurrences in S Africa/W.Angola therefore remove as possible misidentification and unlikely to have large effect on SDM. 
# Remove occurrences in S Africa as potentially misidentified and actually Amb hebraeum (only 2 records)
# Remove occurences from W. Angola as potentially misidenitifed and actually pompsum (only 2 records)
ticks = ticks[ -which(ticks$Species_current == "Amblyomma variegatum" & ticks$Record %in% c(3349, 3345, 4392, 4393))]

# viz: all looks fine
#ggplot() +
    # geom_sf(data=shp0, colour="grey30", fill="grey95", size=0.25) +
    # theme_minimal() +
    # geom_point(data=as.data.frame(ticks[ticks$Species == "Amblyomma variegatum", ]), aes(x, y, col=Species), pch=16, size=0.8, alpha=0.6, col="red") +
    # #scale_color_viridis_d(begin=0, end=0.8, option="inferno") +
    # theme(axis.text = element_blank(),
    #       axis.title=element_blank())

write.csv(as.data.frame(ticks), "./data/africanticks_cumming1998/Main_tick_database(.txt)/RGLO_AllTicks_May2022.csv", row.names = FALSE)






# =============== Cleaning GBIF ===============


  