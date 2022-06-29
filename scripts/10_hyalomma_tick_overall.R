
setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")
library(raster)
library(rgdal)
library(sf)
library(ggplot2)

# shapefiles
shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
shpy = sf::st_union(shpx)

# uganda
shp_uga = shpx[ shpx$NAME_0 == "Uganda", ]

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


# ================== Observed and predicted Hyalomma distribution ==================

# read tick records
# read and filter for each species
ticks = read.csv("./data/africanticks_cumming1998/Main_tick_database(.txt)/RG_Hyalommasubset_23062020.csv", stringsAsFactors=FALSE)
ticks = ticks[ ticks$decade_est >= 1930, ]
ticks$PointsDensity = ticks$PointsDensity_90

# combine with additional occurrences from literature
ticks_rg = read.csv("./data/africanticks_rg/rg_additional_literature.csv", stringsAsFactors = FALSE)
ticks_rg$PointsDensity = 0.001
names(ticks_rg)[1:3] = c("Species_current", "x", "y")
names(ticks_rg)[7] = "Year1"

# filter down
ticks = rbind(ticks[ , which(names(ticks) %in% c("Species_current", "x", "y", "Country", "PointsDensity")) ], 
              ticks_rg[ , which(names(ticks_rg) %in% c("Species_current", "x", "y", "Country", "PointsDensity")) ])
ticks = ticks[ !is.na(ticks$y), ]

# species
ticks = ticks[ ticks$Species_current %in% c("Hyalomma impeltatum", "Hyalomma truncatum", "Hyalomma rufipes", "Hyalomma dromedarii"), ]

# all ticks
all_ticks = read.csv("./data/africanticks_cumming1998/Main_tick_database(.txt)/RG_YearProcessed_23062020.csv")



# ================= predictions ==================

# predicted hyalomma suitability
rr = stack(list.files("output/figures/manuscript/tick_suitability_layers/", pattern="meanpred.tif", full.names=TRUE))[[2:5]]
rr = calc(rr, fun = max)





# ============= plots ==============

sampling_plot = ggplot() + 
  geom_point(data = all_ticks, aes(x, y), col="grey70", size=0.03) +
  geom_point(data=ticks, aes(x, y), col="red", size=0.04) + 
  geom_sf(data=shpy, fill=NA, col="black", size=0.75) + 
  geom_sf(data=shp_uga, fill=NA, col="black", size=0.75) + 
  maptheme 

ras_dat = as.data.frame(rr, xy=TRUE)
suitability_plot = ggplot() + 
  geom_raster(data = ras_dat, aes(x, y, fill=layer)) +
  geom_sf(data=shpy, fill=NA, col="black", size=0.75) + 
  geom_sf(data=shp_uga, fill=NA, col="black", size=0.75) + 
  #scale_fill_viridis_c(option="magma", direction=-1, na.value = "white", limits=c(0, 1)) +
  scale_fill_gradientn(colors=rev(viridisLite::mako(200)), na.value="white", limits=c(0,1)) +
  maptheme +
  theme(legend.title=element_blank())

pc = gridExtra::grid.arrange(sampling_plot, suitability_plot, nrow=1, widths=c(1, 1.2))
pc2 = ggpubr::as_ggplot(pc)  +
  cowplot::draw_plot_label(label = c("a", "b"),
                           fontface = "bold", size = 24,
                           x = c(0.03, 0.5), y = c(0.97, 0.97))
ggsave(pc2, file="output/figures/manuscript/Hyalomma_SI_comb.jpeg", units="in", scale=0.9, width=11, height=5.5, dpi=900)
