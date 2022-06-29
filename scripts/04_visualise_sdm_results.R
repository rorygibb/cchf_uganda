

# ====================== VISUALISE SDMS =============================



# ------------- dependencies and housekeeping --------------

# root dir
setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")
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

# africa shapefile
shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
shpy = readOGR("./data/project_shapefiles/africa_admin0.shp")



# --------------------- paths for models ----------------------------

# path where models are stored
model_path = "./output/model_outputs/climate_ensemble/fitted_models/"
model_path = "D:/ResearchProjects/202002_ucl_congocrimean/analysis/output/model_outputs/climate_ensemble/fitted_models/"

# paths where models are
# filepaths = paste(model_path, c("Hyalomma rufipes_BRT_2020-06-24_743236",
#               "Hyalomma truncatum_BRT_2020-06-24_531097",
#               "Hyalomma dromedarii_BRT_2020-06-24_111728",
#               "Lepus capensis_BRT_2020-06-29_236594" ,
#               "Lepus microtis_BRT_2020-06-30_359946"), "/", sep="")

# filepaths = paste(model_path, c("Hyalomma truncatum_BRT_2020-06-30_666918",
#                                 "Hyalomma truncatum_BRT_2020-07-01_873699",
#                                 "Hyalomma rufipes_BRT_2020-06-30_523116",
#                                 "Hyalomma rufipes_BRT_2020-07-01_479942",
#                                 "Hyalomma impeltatum_BRT_2020-07-01_686854"), "/", sep="")

filepaths = paste(list.files(model_path, full.names=TRUE), "/", sep="")
#filepath="./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-02_841864/"


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
    x[ x == "bio14_precip_driestmonth" ] = "Precip mean driest (bio14"
    x
  }
  vi$Variable = renameVars(vi$Variable)
  pd$var = renameVars(pd$var)
  
  # plot var importance
  vi_summary = vi %>%
    group_by(Variable) %>%
    dplyr::summarise(var = unique(Variable),
                     MeanRelInfluence = mean(RelativeInfluence),
                     LowerRelInfluence = quantile(RelativeInfluence, 0.05),
                     UpperRelInfluence = quantile(RelativeInfluence, 0.95))
  vi_summary = vi_summary[ order(vi_summary$MeanRelInfluence, decreasing = TRUE), ]
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
  
  # plot top 4 importance variables
  pd = pd [ pd$var %in% vi_summary$var[1:4], ]
  pd$var = factor(pd$var, levels=vi_summary$var[1:4], ordered=TRUE)
  pd_plot = ggplot(data=pd, aes(var_value, y, group=model_id)) + 
    geom_line(alpha=0.1, col="skyblue4") +
    #stat_summary(aes(group=model_id), fun=mean, colour="black", geom="line") +
    facet_wrap(~var, scales="free_x") + 
    theme_classic() +
    theme(plot.title = element_text(size=14, hjust=0.5),
          axis.text = element_text(size=11),
          strip.text = element_text(size=11),
          strip.background = element_blank()) +
    geom_hline(yintercept=0, lty=2) +
    xlab("Covariate value") +
    ylab("Fitted function") +
    ggtitle(pd$Species[1])
  
  # read in species occurrences
  occ = read.csv(list.files(paste(filepath, "occurrence_points", sep=""), pattern=".csv", full.names=TRUE), stringsAsFactors = FALSE)
  
  # read raster
  ff = list.files(paste(filepath, "mean_geo", sep=""), pattern="meanpred", full.names=TRUE)
  ff = ff[ -grep(".aux.xml", ff)]
  geo_pred = raster(ff)
  
  # save
  ras = as.data.frame(geo_pred, xy=TRUE)
  names(ras)[3] = "val"
  #colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlGnBu"))(200)
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
          axis.line = element_blank()) +
    #scale_fill_viridis(option = "magma", na.value = "white") +
    scale_fill_gradientn(colors=colscalex, na.value = "white") +
    geom_sf(data=shpx, colour="grey25", fill=NA, size=0.25) +
    geom_point(data=occ, aes(x, y), pch=16, size=0.25, col="red", alpha=0.6)
  
  #
  #plot_composite = gridExtra::grid.arrange(pd_plot, vi_plot, nrow=2, heights=c(1, 1))
  plot_composite = gridExtra::grid.arrange(grobs=list(pd_plot, vi_plot, geo_plot), heights=c(1, 1), widths=c(1,1.5), layout_matrix=rbind(c(1, 3), c(2, 3)))
  dir.create(paste(filepath, "summaryfigures", sep=""))
  ggsave(plot_composite, filename=paste(filepath, "summaryfigures/ensemble_summary.png", sep=""), device="png", units="in", dpi=300, width=18, height=9, scale=0.9)
  
  # combine
  return(plot_composite)
}


# -------------- create figure for each model ------------

# create figures
# createModelFigure("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_658520_lepus/")
# createModelFigure("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/")
# createModelFigure("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-04_233612_lepus/")
# createModelFigure("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/")

# f2 = createModelFigure(filepaths[2])
f3 = createModelFigure(filepaths[3])
# f4 = createModelFigure(filepaths[5])

# run for all files
lapply(filepaths, createModelFigure)

