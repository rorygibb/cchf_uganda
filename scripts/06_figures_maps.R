
# ====================== MAPPING OF FIGURES/RESULTS =============================



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



# ----------------- plot updated human CCHF occurrences -------------------

a1 = read.csv("./data/cchf_occurrences/messsina_db_2013/CCHF_1953_2012_Messina.csv", stringsAsFactors = FALSE)
a1 = a1[ a1$REGION == "Africa" & !a1$COUNTRY %in% c("United Arab Emirates", "Oman", "Madagascar"), ]

a1$Decade = NA
a1$Decade[ a1$YEAR > 1950 & a1$YEAR < 1980 ] = "1958-1979"
a1$Decade[ a1$YEAR >= 1980 & a1$YEAR < 2000 ] = "1980-1999"
a1$Decade[ a1$YEAR >= 2000 ] = "2000-2012"

a2 = a1[ , c("LATITUDE", "LONGITUDE", "Decade")]
names(a2) = c("y", "x", "Decade")
a2$Source = "Messina"
a2$Type = "Human"

a3 = read.csv("./data/cchf_occurrences/rg_db_/cchf_location_datasources_georef.csv", stringsAsFactors = FALSE)
a4 = a3[ , c("Latitude", "Longitude", "Type")]
names(a4) = c("y", "x", "Type")
a4$Source = "RG"
a4$Decade = "2013-2018"

# combine
cchf = rbind(a2, a4)
cchf = cchf[ !is.na(cchf$y), ]
write.csv(cchf, "./output/")

#
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlGnBu"))(200)
cchf_plot = ggplot() + 
  geom_sf(data=shpx, colour="grey30", fill="grey98", size=0.2) +  
  geom_point(data=cchf[cchf$Type == "Human",], aes(x, y, fill=Decade), col="black", size=1.2, pch=21, alpha=0.8) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position=c(0.9, 0.1), 
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_brewer(palette="YlOrRd")

ggsave(cchf_plot, file="./output/figures/report/FigX_CCHF_occurrences_updated_human.png",
       device="png", dpi=600, height=8, width=8, units="in", scale=0.9)


cchf_plot = ggplot() + 
  geom_sf(data=shpx, colour="grey70", fill="grey90", size=0.2) +  
  geom_point(data=cchf[cchf$Type == "Human",], aes(x, y), fill="red", col="black", size=1.5, pch=21, alpha=0.8) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position=c(0.9, 0.1), 
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) 

ggsave(cchf_plot, file="./output/figures/occurrences.png",
       device="png", dpi=600, height=8, width=8, units="in", scale=0.9)


# plot subset for uganda
shp_uga = st_read("./data/project_shapefiles/uga_admin1.shp")
ux = extent(c(xmin=29.5, xmax=36, ymin=-2, ymax=4))
cc2 = cchf; coordinates(cc2) = ~x+y
cc2 = as.data.frame(crop(cc2, ux))

uga_plot = ggplot() + 
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  geom_sf(data=shp_uga, colour="grey60", fill="grey98", size=0.215) + 
  geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.5) +  
  geom_point(data=cc2, aes(x, y, fill=Decade), col="black", size=2, pch=21, alpha=0.8) +
  scale_fill_brewer(palette="YlOrRd")

ggsave(uga_plot, file="./output/figures/report/FigX_CCHF_occurrences_updated_Uganda.png",
       device="png", dpi=600, height=6, width=6, units="in", scale=0.9)





# ----------------------- mapping of ecological hazard for CCHF (hyalomma) -------------------------

r1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/mean_geo/Hyalomma rufipes_meanpred.tif")
r2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/mean_geo/Hyalomma truncatum_meanpred.tif")
# r1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_658520_lepus/mean_geo/Hyalomma rufipes_meanpred.tif")
# r2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-04_233612_lepus/mean_geo/Hyalomma truncatum_meanpred.tif")
r3 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma dromedarii_BRT_2020-06-24_111728_wolepus/mean_geo/Hyalomma dromedarii_meanpred.tif")
r4 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma impeltatum_BRT_2020-07-01_686854_wlepus/mean_geo/Hyalomma impeltatum_meanpred.tif")
haz = stack(r1, r2, r3, r4)
haz = calc(haz, fun = max)

rx = as.data.frame(haz, xy=TRUE)
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlGnBu"))(200)
# colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlOrRd"))(200)
# colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)

haz_plot = ggplot() + 
  geom_raster(data=rx, aes(x, y, fill=layer)) +
  geom_sf(data=shpx, colour="grey30", fill=NA, size=0.15) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position=c(0.9, 0.1), 
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscalex, na.value = "white")

ggsave(haz_plot, file="./output/figures/report/FigX_CCHF_hazard_model.png",
       device="png", dpi=600, height=8, width=8, units="in", scale=0.9)

hp2 = haz_plot + 
  geom_point(data=cchf, aes(x, y), col="red", size=1.2, pch=21, fill=NA, alpha=0.8)

ggsave(hp2, file="./output/figures/report/FigX_CCHF_hazard_model_points_env.png",
       device="png", dpi=600, height=8, width=8, units="in", scale=0.9)


# ------------------------ Uganda specific --------------------------

# read raster
geo_pred1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/mean_geo/Hyalomma rufipes_meanpred.tif")
geo_pred2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/mean_geo/Hyalomma truncatum_meanpred.tif")

# occ points
occ1 = read.csv("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/occurrence_points/Hyalomma rufipes_occurrences.csv", stringsAsFactors = FALSE)
occ2 = read.csv("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/occurrence_points/Hyalomma truncatum_occurrences.csv", stringsAsFactors = FALSE)
occ = rbind(occ1, occ2)

colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="YlGnBu"))(200)

# plot subset for uganda
shp_uga = st_read("./data/project_shapefiles/uga_admin1.shp")
ux = extent(c(xmin=29, xmax=36, ymin=-2, ymax=4.5))
cc2 = occ; coordinates(cc2) = ~x+y
cc2 = as.data.frame(crop(cc2, ux))

# 1. for rufipes
gp_uga = crop(geo_pred1, ux)
gp_uga = as.data.frame(gp_uga, xy=TRUE)
names(gp_uga)[3] = "var"
uga_plot1 = ggplot() + 
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  ggtitle("Hyalomma rufipes") +
  geom_raster(data=gp_uga, aes(x, y, fill=var)) +
  #geom_sf(data=shp_uga, colour="grey60", fill=NA, size=0.215) + 
  geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.5) +  
  geom_point(data=cc2[ cc2$Species_current == "Hyalomma rufipes", ], aes(x, y), col="red", fill="coral", size=2, pch=21, alpha=0.8) +
  scale_fill_gradientn(colors=colscalex, na.value = "white")

# 2. for truncatum
gp_uga = crop(geo_pred2, ux)
gp_uga = as.data.frame(gp_uga, xy=TRUE)
names(gp_uga)[3] = "var"
uga_plot2 = ggplot() + 
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  geom_raster(data=gp_uga, aes(x, y, fill=var)) +
  ggtitle("Hyalomma truncatum") +
  #geom_sf(data=shp_uga, colour="grey60", fill=NA, size=0.215) + 
  geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.5) +  
  geom_point(data=cc2[ cc2$Species_current == "Hyalomma truncatum", ], aes(x, y), col="red", fill="coral", size=2, pch=21, alpha=0.8) +
  scale_fill_gradientn(colors=colscalex, na.value = "white")

ppp = grid.arrange(uga_plot1, uga_plot2, ncol=2)

ggsave(ppp, filename="./output/figures/report/Hyalomma_ugandamap.png", device="png", units="in", dpi=600, width=9, height=4.5, scale=0.9)



colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)

# plot subset for uganda
shp_uga = st_read("./data/project_shapefiles/uga_admin1.shp")
ux = extent(c(xmin=29.5, xmax=36, ymin=-2, ymax=4))
cc2 = cchf; coordinates(cc2) = ~x+y
cc2 = as.data.frame(crop(cc2, ux))

# crop hazard raster
haz_uga = crop(haz, ux)
haz_uga = as.data.frame(haz_uga, xy=TRUE)

uga_plot1 = ggplot() + 
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  geom_raster(data=haz_uga, aes(x, y, fill=layer)) +
  #geom_sf(data=shp_uga, colour="grey60", fill=NA, size=0.215) + 
  geom_sf(data=st_crop(shpx, ux), colour="grey10", fill=NA, size=0.5) +  
  geom_point(data=cc2, aes(x, y), col="red", size=2, pch=21, alpha=0.8) +
  scale_fill_gradientn(colors=colscalex, na.value = "white")


#
r1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/mean_geo/Hyalomma rufipes_meanpred.tif")
r2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/mean_geo/Hyalomma rufipes_sdpred.tif")
plot(r2, col=colscalex)
plot(crop(r2 / sqrt(r1), ux) )


r1 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/mean_geo/Hyalomma truncatum_meanpred.tif")
r2 = raster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/mean_geo/Hyalomma truncatum_sdpred.tif")
plot(crop(r2 / sqrt(r1), ux) )

r1 = raster("./output/model_outputs/cchf_brt/fitted_models/cchf_10km_BRT_2020-07-07_563371/mean_geo/CCHF_meanpred.tif")
r2 = raster("./output/model_outputs/cchf_brt/fitted_models/cchf_10km_BRT_2020-07-07_563371/mean_geo/CCHF_sdpred.tif")
plot(r2 / sqrt(r1), col=colscalex )
