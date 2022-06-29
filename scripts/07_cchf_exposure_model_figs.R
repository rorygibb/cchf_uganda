

# ================== Composite plot of CCHF exposure models for manuscript ====================

setwd("D:/ResearchProjects/202002_ucl_crimeancongo/analysis/")
#setwd("C:/Users/roryj/Dropbox/CrimeanCongo_Uganda/")
setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")

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
library(exactextractr)
library(dplyr)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)

# extent
af_ext = extent(c(-20, 60, -40, 37))

# colscale
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)

# shapefiles
shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
shpy = readOGR("./data/project_shapefiles/africa_admin0.shp")

# africa overall border
shpz = sf::st_union(shpx)



# ----------------------- create figures ----------------------------

filepathx = "./output/model_outputs/cchf_brt/fitted_models/cchf_10km_exp_BRT_2022-05-20_926065/"



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
    x[ x == "bio05_temp_maxwarmestmonth" ] = "Temp max warmest"
    x[ x == "bio06_temp_mincoldestmonth" ] = "Temp min coldest"
    x[ x == "NDVI_AnnualMean" ] = "NDVI annual mean"
    x[ x == "NDVI_AnnualSD" ] = "NDVI seasonality"
    x[ x == "Lepus_Suitability" ] = "Lepus suitability"
    x[ x == "bio02_temp_meandiurnalrange" ] = "Temp diurnal range (bio02)"
    x[ x == "bio13_precip_wettestmonth" ] = "Precip mean wettest (bio13)"
    x[ x == "bio14_precip_driestmonth" ] = "Precip mean driest (bio14)"
    x[ x == "Sheep_logDensity" ] = "Sheep density (log)"
    x[ x == "Cattle_logDensity" ] = "Cattle density (log)"
    x[ x == "Cropland_proportion" ] = "Agricultural land prop."
    x[ x == "Lepus_suitability" ] = "Lepus suitability"
    x[ x == "Hyalomma_suitability" ] = "Hyalomma suitability"
    x[ x == "Rhipicephalus_suitability" ] = "Rhipicephalus suitability"
    x[ x == "Amblyomma_suitability" ] = "Amblyomma suitability"
    x[ x == "Goats_logDensity" ] = "Goats density (log)"
    
    x
  }
  
  renameVars2 = function(x){
    x[ x == "bio12_precip_annualtotal" ] = "Precip mean"
    x[ x == "bio01_temp_annualmean" ] = "Temp mean"
    x[ x == "bio04_temp_seasonality" ] = "Temp seas"
    x[ x == "bio15_precip_seasonality" ] = "Precip seas"
    x[ x == "bio05_temp_maxwarmestmonth" ] = "Tmax warmest "
    x[ x == "bio06_temp_mincoldestmonth" ] = "Tmin coldest"
    x[ x == "NDVI_AnnualMean" ] = "NDVI mean"
    x[ x == "NDVI_AnnualSD" ] = "NDVI seas"
    x[ x == "Lepus_Suitability" ] = "Lepus"
    x[ x == "bio02_temp_meandiurnalrange" ] = "Tdrange"
    x[ x == "bio13_precip_wettestmonth" ] = "Precip wettest"
    x[ x == "bio14_precip_driestmonth" ] = "Precip driest"
    x[ x == "Sheep_logDensity" ] = "Sheep"
    x[ x == "Cattle_logDensity" ] = "Cattle"
    x[ x == "Cropland_proportion" ] = "Agri prop."
    x[ x == "Lepus_suitability" ] = "Lepus"
    x[ x == "Hyalomma_suitability" ] = "Hyalomma"
    x[ x == "Rhipicephalus_suitability" ] = "Rhipicephalus"
    x[ x == "Amblyomma_suitability" ] = "Amblyomma"
    x[ x == "Goats_logDensity" ] = "Goats"
    
    x
  }
  
  vi$Variable1 = renameVars(vi$Variable)
  vi$Variable2 = renameVars2(vi$Variable)
  pd$var = renameVars(pd$var)
  
  # plot var importance
  vi_summary = vi %>%
    dplyr::group_by(Variable1) %>%
    dplyr::summarise(var = unique(Variable1),
                     Variable2 = unique(Variable2),
                     MedRelInfluence = median(RelativeInfluence),
                     MeanRelInfluence = mean(RelativeInfluence),
                     LowerRelInfluence = quantile(RelativeInfluence, 0.05),
                     UpperRelInfluence = quantile(RelativeInfluence, 0.95))
  vi_summary = vi_summary[ order(vi_summary$MedRelInfluence, decreasing = TRUE), ]
  vi$Variable1 = factor(vi$Variable1, levels=as.vector(unique(vi_summary$var)), ordered=TRUE)
  vi$Variable2 = factor(vi$Variable2, levels=as.vector(unique(vi_summary$Variable2)), ordered=TRUE)
  vi_plot = ggplot(vi, aes(Variable2, RelativeInfluence)) + 
    ggforce::geom_sina(width=1, pch=21, fill="grey70", col="grey20", alpha=0.1) +
    geom_boxplot(width=0.4, fill="skyblue4", alpha=0.5, outlier.shape = NA) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, size=12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=11),
          axis.title.y = element_text(size=12),
          plot.title = element_text(size=15, hjust=0.5)) +
    ylab("Relative variable importance (%)") 
  
  # plot top 8 importance variables
  pd = pd[ pd$var %in% vi_summary$var[1:8], ]
  pd$var = factor(pd$var, levels=vi_summary$var[1:8], ordered=TRUE)
  
  #
  mean_func = pd %>% 
    dplyr::group_by(var) %>%
    dplyr::mutate(var_bin = INLA::inla.group(var_value, n = 50)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(var, var_bin) %>%
    dplyr::summarise(mean_func = mean(y))
    
  pd_plot = ggplot() + 
    geom_line(data=pd, aes(var_value, y, group=model_id), alpha=0.05, col="skyblue4") +
    geom_line(data=mean_func, aes(var_bin, mean_func), col="black", size=1.65) + 
    geom_line(data=mean_func, aes(var_bin, mean_func), col="yellow", size=0.5) + 
    #stat_summary(aes(group=model_id), fun=mean, colour="black", geom="line") +
    facet_wrap(~var, scales="free_x", nrow=2) + 
    theme_classic() +
    theme(plot.title = element_text(size=14, hjust=0.5),
          axis.text.y = element_text(size=11),
          axis.text.x = element_text(size=10),
          axis.title = element_text(size=12),
          strip.text = element_text(size=11),
          strip.background = element_blank()) +
    #geom_hline(yintercept=0, lty=2) +
    xlab("Covariate value") +
    ylab("Fitted function")
  
  #plot_composite = gridExtra::grid.arrange(pd_plot, vi_plot, nrow=2, heights=c(1, 1))
  plot_composite = gridExtra::grid.arrange(grobs=list(pd_plot, vi_plot), heights=c(1.5, 1), nrow=2)
  
  dir.create(paste(filepath, "summaryfigures", sep=""))
  ggsave(plot_composite, filename=paste(filepath, "summaryfigures/ensemble_model_summary.png", sep=""), device="png", units="in", dpi=300, width=10, height=9, scale=0.85)
  
  # combine
  return(plot_composite)
}

# save
f1 = createModelFigure(filepathx)



# -------------------- create -------------------

createModelMap = function(filepath){
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # read raster
  geo_pred = raster(paste(filepath, "mean_geo/CCHF_medianpred.tif", sep=""))
  
  #colscalex = rev(viridisLite::mako(200))
  colscalex = rev(viridisLite::magma(200))
  
  # save
  ras = as.data.frame(geo_pred, xy=TRUE)
  names(ras)[3] = "val"
  geo_plot = ggplot() +
    geom_raster(data=ras, aes(x, y, fill=val)) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust=0.5, size=18, vjust=-3),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          #panel.border = element_rect(color="grey70", fill=NA),
          panel.border = element_blank(),
          axis.line = element_blank()) +
    #scale_fill_viridis(option = "magma", na.value = "white") +
    scale_fill_gradientn(colors=colscalex, na.value = "white", limits=c(0, 1)) +
    #geom_point(data=occ, aes(x, y), pch=16, size=0.75, col="red") +
    geom_sf(data=shpz, colour="grey25", fill=NA, size=0.8) +
    geom_sf(data=shpx[ shpx$NAME_0 == "Uganda", ], colour="blue", fill=NA, size=0.8)
  
  
  #ggsave(geo_plot, filename=paste(filepath, "summaryfigures/riskmap.png", sep=""), device="png", units="in", dpi=600, width=9, height=9, scale=0.9)
  return(geo_plot)
}

# createUncertMap = function(filepath){
#   
#   # read in species occurrences
#   occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
#   
#   # relative uncertainty: sd / sqrt(mean), because uncert expected to be higher when mean is higher
#   geo_pred1 = raster(paste(filepath, "mean_geo/CCHF_sdpred.tif", sep=""))
#   geo_pred2 = raster(paste(filepath, "mean_geo/CCHF_meanpred.tif", sep=""))
#   ru_pred = geo_pred1 / sqrt(geo_pred2)
#   
#   # save
#   ras = as.data.frame(ru_pred, xy=TRUE)
#   names(ras)[3] = "val"
#   geo_plot = ggplot() +
#     geom_raster(data=ras, aes(x, y, fill=val)) +
#     theme_classic() +
#     theme(legend.title = element_blank(),
#           legend.position = c(0.25, 0.25),
#           plot.title = element_text(hjust=0.5, size=18, vjust=-3),
#           axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.line = element_blank(),
#           panel.border = element_rect(color="black", fill=NA)) +
#     #scale_fill_viridis(option = "magma", na.value = "white") +
#     scale_fill_gradientn(colors=colscalex, na.value = "white") +
#     geom_sf(data=shpx, colour="grey25", fill=NA, size=0.25) +
#     geom_point(data=occ, aes(x, y), pch=16, size=0.75, col="red")
#   
#   ggsave(geo_plot, filename=paste(filepath, "summaryfigures/riskmap_sd.png", sep=""), device="png", units="in", dpi=600, width=9, height=9, scale=0.9)
# }
# 
# createUncertMap(filepathx)

createUgandaMap = function(filepath){
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # read raster
  geo_pred = raster(paste(filepath, "mean_geo/CCHF_medianpred.tif", sep=""))
  
  # plot subset for uganda
  shp_uga = st_read("./data/project_shapefiles/uga_admin1.shp")
  ux = extent(c(xmin=29, xmax=36, ymin=-2, ymax=4.5))
  cc2 = occ; coordinates(cc2) = ~x+y
  cc2 = as.data.frame(crop(cc2, ux))
  
  # crop hazard raster
  gp_uga = crop(geo_pred, ux)
  gp_uga = as.data.frame(gp_uga, xy=TRUE)
  
  #colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="GnBu"))(200)
  colscalex = rev(viridisLite::magma(200))
  #colscalex = viridisLite::mako(200)
  
  uga_plot1 = ggplot() + 
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=11),
          plot.title = element_text(hjust=0.5, vjust=-3),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + 
    geom_raster(data=gp_uga, aes(x, y, fill=CCHF_medianpred)) +
    #geom_sf(data=shp_uga, colour="grey60", fill=NA, size=0.215) + 
    geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.8) +  
    geom_point(data=cc2, aes(x, y), col="lightblue", fill="cyan2", size=1.8, pch=21, alpha=0.7) +
    scale_fill_gradientn(colors=colscalex, na.value = "white", limits=c(0, 1))
  
  #ggsave(uga_plot1, filename=paste(filepath, "summaryfigures/uga_hyalommamap.png", sep=""), device="png", units="in", dpi=600, width=6, height=6, scale=0.9)
  return(uga_plot1)
}

p1 = createModelMap(filepathx)
p2 = createUgandaMap(filepathx)

maps = gridExtra::grid.arrange(p2, p1, nrow=2, heights=c(0.9, 1))


# combine
# p_comb = gridExtra::grid.arrange(f1, maps, nrow=1, widths=c(1, 0.6))
# p_comb2 = ggpubr::as_ggplot(p_comb)  +
#   cowplot::draw_plot_label(label = c("a", "b", "c", "d"),
#                            fontface = "bold", size = 24,
#                            x = c(0.005, 0.04, 0.66, 0.68), y = c(0.99, 0.42, 0.99, 0.48))
# ggsave(p_comb2, file="./output/figures/paper/Figure4_exposuremodels.jpeg", units="in", scale=0.9, width=14, height=8.5, dpi=600)


p_comb = gridExtra::grid.arrange(maps, f1, nrow=1, widths=c(0.6, 1))
p_comb2 = ggpubr::as_ggplot(p_comb)  +
  cowplot::draw_plot_label(label = c("a", "b", "c", "d"),
                           fontface = "bold", size = 24,
                           x = c(-0.005, 0.01, 0.36, 0.35), y = c(0.99, 0.47, 0.99, 0.43))
ggsave(p_comb2, file="./output/figures/paper/Figure4_exposuremodels.jpeg", units="in", scale=0.9, width=14, height=8.5, dpi=600)

