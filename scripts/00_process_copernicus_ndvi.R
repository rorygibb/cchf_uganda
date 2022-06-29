



setwd("C:/Users/roryj/Dropbox/CrimeanCongo_Uganda/")

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
library(embarcadero)
library(sf)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)

# africa_ext
af_ext = extent(c(-20, 60, -40, 40))

# shapefile and mask for raster
ss= readOGR("./data/project_shapefiles/africa_admin0.shp")
gm = raster("./data/africa_maskraster/africa_mask.tif")


# -------------

output_loc = "./output/ndvi_copernicus/africa_lts/mean/"
file_locs = list.files("./data/ndvi_copernicus/longtermstatistics_ndvi/", pattern=".nc", full.names = TRUE)
file_names = list.files("./data/ndvi_copernicus/longtermstatistics_ndvi/", pattern=".nc", full.names = FALSE)

# process files and save
for(i in 1:length(file_locs)){
  
  #reporting
  cat(paste(i, "..."))
  
  # mean NDVI for time period
  loc = file_locs[i]
  rasx = stack(loc, varname="min")
  #print(rasx)
  
  # crop and reclassify ocean cells
  rasx = crop(rasx, af_ext)
  rasx = reclassify(rasx, cbind(0.936, 1, NA), right=FALSE)
  
  #aggregate by factor of 4 // resample to africa mask 
  #rasx = aggregate(rasx, fact=4, fun="mean")
  rasx = resample(rasx, gm, method="bilinear")
  
  # updated name
  namex = strsplit(file_names[i], "_")[[1]][4]
  namex = strsplit(namex, "-")[[1]][3]
  names(rasx) = paste("Copernicus_NDVImin-LTS_1992-2017_AfricaResamp_", namex, sep="")
  
  # save raster
  output_name = paste(output_loc, names(rasx), ".tif", sep="")
  writeRaster(rasx, file=output_name, format="GTiff", overwrite=TRUE)
}


# -------------- create summary metrics: mean and SD --------------

# read in ndvi stack
rr = stack(list.files(output_loc, full.names = TRUE, pattern=".tif"))

# pixelwise mean and sd
ndvi_mean = mean(rr)
ndvi_sd = calc(rr, fun=sd)

# save
writeRaster(ndvi_mean, "./output/ndvi_copernicus/africa_summary/NDVI_Copernicus_LTS19922017_DekadMean.tif", format="GTiff")
writeRaster(ndvi_sd, "./output/ndvi_copernicus/africa_summary/NDVI_Copernicus_LTS19922017_DekadSD.tif", format="GTiff")


# to calculate pixelwise minimum
ndvi_min = calc(rr, fun=min)
writeRaster(ndvi_min, "./output/ndvi_copernicus/africa_summary/NDVI_Copernicus_LTS19922017_MinimumDekad_longterm.tif", format="GTiff", overwrite=TRUE)


# plot and save
ss = st_read("./data/project_shapefiles/africa_admin0.shp")
nm = as.data.frame(ndvi_min, xy=TRUE)
p1 = ggplot() + 
  geom_raster(data=nm, aes(x, y, fill=layer)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis(option = "inferno", na.value = "white") +
  geom_sf(data=ss, colour="black", fill=NA, size=0.25) + 
  ggtitle("NDVI annual minimum (long term mean)") 
ggsave(p1, file="./output/figures/draft/Africa_longtermNDVI_minimum_copernicus2.tif", device="tiff", dpi=600, width=10, height=12, units="in")

nsd = as.data.frame(ndvi_sd, xy=TRUE)
p2 = ggplot() + 
  geom_raster(data=nsd, aes(x, y, fill=layer)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis(option = "inferno", na.value = "white") +
  geom_sf(data=ss, colour="black", fill=NA, size=0.25) + 
  ggtitle("NDVI long term SD") 
ggsave(p2, file="./output/figures/draft/Africa_longtermNDVI_copernicus_sd.tif", device="tiff", dpi=600, width=10, height=12, units="in")

