

# ============ Preliminary processing of spatial/shapefile data ============

# wd
setwd("C:/Users/roryj/Documents/PhD/project_cchf_africa/analysis/")

# dependencies
library(raster)
library(rgdal)
library(countrycode)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)



# ============ prepare and save shapefiles for Africa, Uganda at country, admin1 and admin2levels ==============

# ---------- Africa ------------

# GADM level 0: countries
shp0 = st_read("./data/gadm_36/gadm36_levels_shp/gadm36_0.shp")
shp0$Continent = countrycode(as.vector(shp0$GID_0), origin="iso3c", destination="continent")
shp0 = shp0[ shp0$Continent == "Africa" & !is.na(shp0$Continent), ]

# GADM admin level 1 (states/provinces)
shp1 = st_read("./data/gadm_36/gadm36_levels_shp/gadm36_1.shp")
shp1$Continent = countrycode(as.vector(shp1$GID_0), origin="iso3c", destination="continent")
shp1 = shp1[ shp1$Continent == "Africa" & !is.na(shp1$Continent), ]

# GADM admin level 2: local government/districts
shp2 = st_read("./data/gadm_36/gadm36_levels_shp/gadm36_2.shp")
shp2$Continent = countrycode(as.vector(shp2$GID_0), origin="iso3c", destination="continent")
shp2 = shp2[ shp2$Continent == "Africa" & !is.na(shp2$Continent), ]

# save all
st_write(shp0, dsn="./data/project_shapefiles/africa_admin0.shp", layer="shp0", driver="ESRI Shapefile")
st_write(shp1, dsn="./data/project_shapefiles/africa_admin1.shp", layer="shp1", driver="ESRI Shapefile")
st_write(shp2, dsn="./data/project_shapefiles/africa_admin2.shp", layer="shp2", driver="ESRI Shapefile")



# ----------- Uganda -------------

# spatial extent:
ux = extent(c(xmin=25, xmax=39, ymin=-5, ymax=9))
ux = extent(c(xmin=27, xmax=37, ymin=-3, ymax=6))

# subset 
uga_admin0 = shp0[ shp0$NAME_0 == "Uganda", ]
uga_admin1 = shp1[ shp1$GID_0 == "UGA", ]
uga_admin2 = shp2[ shp2$GID_0 == "UGA", ]

# plot uganda
ggplot() +
  geom_sf(data=uga_admin2, colour="coral", fill="grey96") +
  geom_sf(data=uga_admin1, colour="darkred", fill=NA) +
  geom_sf(data=uga_admin0, colour="black", fill=NA) +
  theme_classic()
  
# plot uganda and surrounds
ggplot() +
  geom_sf(data=uga_admin1, colour="coral2", fill=NA) +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA) +
  theme_classic()

# save
st_write(uga_admin0, dsn="./data/project_shapefiles/uga_admin0.shp", layer="uga_admin0", driver="ESRI Shapefile")
st_write(uga_admin1, dsn="./data/project_shapefiles/uga_admin1.shp", layer="uga_admin1", driver="ESRI Shapefile")
st_write(uga_admin2, dsn="./data/project_shapefiles/uga_admin2.shp", layer="uga_admin2", driver="ESRI Shapefile")




# ================= preliminary mapping and visualisation of key env/socioec variables ==============

# spatial extent for uganda
ux = extent(c(xmin=27, xmax=37, ymin=-3, ymax=5.5))
#ux = extent(c(xmin=20, xmax=45, ymin=-10, ymax=14))
colscalex = colorRampPalette(brewer.pal(9, name="YlGnBu"))(200)

# specify shapefile border size
shp_size = 0.7


# ------------ africa CCHF occurrences ----------------

cchf = read.csv("./data/cchf_occurrences/messsina_db_2013/CCHF_1953_2012_Messina.csv", stringsAsFactors = FALSE)
cchf= cchf[ cchf$REGION == "Africa" & !cchf$COUNTRY %in% c("United Arab Emirates", "Oman"), ]
cchf$Decade = NA
cchf$Decade[ cchf$YEAR > 1950 & cchf$YEAR < 1980 ] = "1958-1979"
cchf$Decade[ cchf$YEAR >= 1980 & cchf$YEAR < 2000 ] = "1980-1999"
cchf$Decade[ cchf$YEAR >= 2000 ] = "2000-2013"

cchf_plot = ggplot() + 
  geom_sf(data=shp0, colour="grey30", fill=NA, size=0.25) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.position=c(0.9, 0.1), 
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  geom_point(data=cchf, aes(LONGITUDE, LATITUDE, fill=Decade), pch=21)
  
ggsave(cchf_plot, file="./output/figures/draft/CCHF_distribution.png",
       device="png", dpi=600, height=8, width=8, units="in", scale=0.9)



# ------------ climate variables (CHELSA) ----------------

# create dataframe to spatially plot raster
cc = stack(list.files("./data/chelsa_climatologies/africa_rescaled_aggregated/", pattern=".tif", full.names=TRUE))
cc = crop(cc, ux)
ccx = as.data.frame(cc, xy=TRUE)

p1 = ggplot() + 
  geom_raster(data=ccx, aes(x, y, fill=CHELSA_bio10_01_land_temp_annualmean_rescaled)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscalex, na.value = "white") +
  #geom_sf(data=uga_admin1, colour="black", fill=NA) +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Annual mean temperature (C)") 

p2 = ggplot() + 
  geom_raster(data=ccx, aes(x, y, fill=CHELSA_bio10_12_land_precip_annualtotal)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscalex, na.value = "white") +
  #geom_sf(data=uga_admin1, colour="black", fill=NA) +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Annual mean precipitation") 

p3 = ggplot() + 
  geom_raster(data=ccx, aes(x, y, fill=CHELSA_bio10_15_land_precip_seasonality)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscalex, na.value = "white") +
  #geom_sf(data=uga_admin1, colour="black", fill=NA) +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Precipitation seasonality") 



# ================ livestock ==============

l1 = raster("./data/griddedlivestock_v3/cattle_2010/5_Ct_2010_Da.tif")
l2 = raster("./data/griddedlivestock_v3/goats_2010/5_Gt_2010_Da.tif")
l3 = raster("./data/griddedlivestock_v3/pigs_2010/5_Pg_2010_Da.tif")
ll = stack(l1, l2, l3)
ll = stack(ll, area(ll[[1]]))
names(ll) = c("Cattle", "Goats", "Pigs", "Area")
ll = crop(ll, ux)
llx = as.data.frame(ll, xy=TRUE)

p4 = ggplot() + 
  geom_raster(data=llx, aes(x, y, fill=Cattle/Area)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option="magma", direction = -1, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Cattle density") 

p5 = ggplot() + 
  geom_raster(data=llx, aes(x, y, fill=Goats/Area)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option="magma", direction = -1, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Goats density") 



# =============== Demographic/socioeconomic indicators ===============

# human density
hh = raster("./data/ciesin_griddedpopworld_v4/gpw_v4_population_count_rev11_2015_2pt5_min.tif")
hh = crop(hh, ux)
hhx = as.data.frame(hh, xy=TRUE)
hhx$gpw_v4_population_count_rev11_2015_2pt5_min[ hhx$gpw_v4_population_count_rev11_2015_2pt5_min == 0 ] = NA

p6 = ggplot() + 
  geom_raster(data=hhx, aes(x, y, fill=log(gpw_v4_population_count_rev11_2015_2pt5_min))) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option="magma", direction = -1, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Human population 2015 (log)") 


# poverty
pp = raster("./data/poverty_inequality/worldpop_poverty2011_uganda/uga11povmpi.tif")
pp = crop(pp, ux)
ppx = as.data.frame(pp, xy=TRUE)

p7 = ggplot() + 
  geom_raster(data=ppx, aes(x, y, fill=uga11povmpi)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_viridis_c(option="magma", direction = -1, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Proportion in poverty (Uganda, MPI)") 



# ==================== land cover / land use =========================

colscaley = colorRampPalette(brewer.pal(9, name="BuGn"))(200)


# ------------- ESA-CCI -----------------

lu1 = raster("./data/landcover_landuse/esa_cci/esa_africa_2015.tif")
lu1 = crop(lu1, extent(shp0))

aggregateESA = function(esa_ras, fact, classes){
  propx = function(x, ...){ sum(x %in% classes) / sum(!is.na(x)) }
  rx = aggregate(esa_ras, fact=fact, propx)
  rx
}

# crop = aggregateESA(lu1, 6, classes=10:20)
# mosaic = aggregateESA(lu1, 6, classes=30:40)
cropmosaic = aggregateESA(lu1, 15, classes=10:40)
forest = aggregateESA(lu1, 15, classes=c(50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100, 160, 170))
shrub_grass = aggregateESA(lu1, 15, classes=c(130, 120:122))



lcras = stack(crop, mosaic, forest, shrub_grass)
names(lcras) = c("Cropland", "Mosaic", "Forest", "Shrub/grassland")
lcx = as.data.frame(lcras, xy=TRUE)

p8 = ggplot() + 
  geom_raster(data=lcx, aes(x, y, fill=Cropland)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscaley, na.value = "white") +
  #geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Cropland (proportion grid cell)") 

p9 = ggplot() + 
  geom_raster(data=lcx, aes(x, y, fill=Mosaic)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscaley, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Crop-natural mosaic (proportion grid cell)") 

p10 = ggplot() + 
  geom_raster(data=lcx, aes(x, y, fill=Forest)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscaley, na.value = "white") +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Forest (proportion grid cell)") 



# ------------- Hoskins land use ----------------

# pasture
p1 = raster("./data/landcover_landuse/hoskins_landuse_downscaled/PAS_1km_2005_0ice.bil")
crs(p1) = crs(pp)
p1 = crop(p1, ux)
px = as.data.frame(p1, xy=TRUE)

p11 = ggplot() + 
  geom_raster(data=px, aes(x, y, fill=PAS_1km_2005_0ice)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5, vjust=-3),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_fill_gradientn(colors=colscaley, na.value = "white", limits=c(0, 1)) +
  geom_sf(data=st_crop(shp0, ux), colour="black", fill=NA, size=shp_size) + 
  ggtitle("Rangeland (proportion grid cell)") 



# combine
ras_plot = grid.arrange(p1, p2, p3, p4, p6, p7, p8, p9, p10, ncol=3)
ggsave(ras_plot, file="./output/figures/draft/Environmental_spatial_variables.png",
       device="png", dpi=600, height=13, width=16, units="in", scale=0.9)

