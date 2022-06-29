

# ============================ Produces spatial distribution figures for Hyalomma truncatum and rufipes for Africa and Uganda ==============================

setwd("D:/ResearchProjects/202002_ucl_crimeancongo/analysis/")
#setwd("C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/")
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
library(embarcadero)
library(sf)
library(ENMeval)
library(exactextractr)

# raster temp files - remove after 0.5 hours
removeTmpFiles(0.5)

# extent
af_ext = extent(c(-20, 60, -40, 37))

# colscale
colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, name="PuBuGn"))(200)

# shapefiles
shpx = st_read("./data/project_shapefiles/africa_admin0.shp")
shpy = readOGR("./data/project_shapefiles/africa_admin0.shp")

# shpy
shpz = sf::st_union(shpx)

# read in predicted Hyalomma maps
hruf1 = raster("./output/figures/manuscript/tick_suitability_layers/Hyalomma rufipes_meanpred.tif")
hruf2 = raster("./output/figures/manuscript/tick_suitability_layers/Hrufipes_mean_uga_pred.tif")
htru1 = raster("./output/figures/manuscript/tick_suitability_layers/Hyalomma truncatum_meanpred.tif")
htru2 = raster("./output/figures/manuscript/tick_suitability_layers/Htruncatum_mean_uga_pred.tif")

# read in predicted tick maps (africawide)
ticks_af = stack(list.files("./output/figures/manuscript/tick_suitability_layers/", pattern="meanpred.tif", full.names=TRUE))
ticks_af = ticks_af[[ -c(2:3) ]]
names(ticks_af) = c("Amblyomma variegatum", "Hyalomma rufipes", "Hyalomma truncatum", "Rhipicephalus decoloratus", 
                    "Rhipicephalus appendiculatus")

# ticks for Uganda
ticks_uga = stack(list.files("./output/figures/manuscript/tick_suitability_layers/", pattern="uga_pred.tif", full.names=TRUE))
ticks_uga = ticks_uga[[ -c(2) ]]
names(ticks_uga) = c("Amblyomma variegatum", "Hyalomma rufipes", "Hyalomma truncatum", "Rhipicephalus decoloratus", 
                    "Rhipicephalus appendiculatus")

# shapefile for uganda
shp_uga = st_crop(shpx, ticks_uga)

# # occurrence points
# hruf_oc = read.csv("./output/figures/manuscript/hyalomma_suitability_layers/Hyalomma rufipes_occurrences.csv")
# coordinates(hruf_oc) = ~x+y
# hruf_oc = crop(hruf_oc, shp_uga)
# htru_oc = read.csv("./output/figures/manuscript/hyalomma_suitability_layers/Hyalomma truncatum_occurrences.csv")
# coordinates(htru_oc) = ~x+y
# htru_oc = crop(htru_oc, shp_uga)

# decide on a colour scheme
#colscalex = colorRampPalette(RColorBrewer::brewer.pal(9, "PuBu"))(200)
colscalex = rev(viridisLite::mako(200))



# ============= Uganda maps ================

ras2 = as.data.frame(ticks_uga, xy=TRUE) %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(
    variable = gsub("[.]", "\n", variable),
    variable = replace(variable, variable == "Rhipicephalus\ndecoloratus", "Rhipicephalus (Boophilus)\ndecoloratus"),
    variable = factor(variable, levels=unique(variable)[c(3, 2, 5, 4, 1)], ordered=TRUE)
  )

p1 = ggplot() +
  geom_raster(data=ras2, aes(x, y, fill=value)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust=0.5, size=13, face="italic"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_rect(color="grey70", fill=NA),
        panel.border = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12)) +
  #scale_fill_viridis(option = "magma", na.value = "white") +
  scale_fill_gradientn(colors=colscalex, na.value = "white", limits=c(0, 1)) +
  #scale_fill_viridis_c(option="magma", direction=-1, na.value = "white", limits=c(0, 1)) +
  #geom_point(data=occ, aes(x, y), pch=16, size=0.75, col="red") +
  geom_sf(data=shp_uga, colour="grey10", fill=NA, size=0.8) +  
  facet_wrap(~variable, nrow=1)


# ================ Africa wide maps =================

ras1 = as.data.frame(ticks_af, xy=TRUE) %>%
  reshape2::melt(id.vars = 1:2) %>%
  dplyr::mutate(
    variable = gsub("[.]", "\n", variable),
    variable = replace(variable, variable == "Rhipicephalus\ndecoloratus", "Rhipicephalus (Boophilus)\ndecoloratus"),
    variable = factor(variable, levels=unique(variable)[c(3, 2, 5, 4, 1)], ordered=TRUE)
  )

p2 = ggplot() +
  geom_raster(data=ras1, aes(x, y, fill=value)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust=0.5, size=13, face="italic"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_rect(color="grey70", fill=NA),
        panel.border = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  #scale_fill_viridis(option = "magma", na.value = "white") +
  scale_fill_gradientn(colors=colscalex, na.value = "white", limits=c(0, 1)) +
  #scale_fill_viridis_c(option="magma", direction=-1, na.value = "white", limits=c(0, 1)) +
  #geom_point(data=occ, aes(x, y), pch=16, size=0.75, col="red") +
  geom_sf(data=shpz, colour="grey25", fill=NA, size=0.8) +
  geom_sf(data=shpx[ shpx$NAME_0 == "Uganda", ], colour="red", fill=NA, size=0.8) + 
  facet_wrap(~variable, nrow=1)




# combined plot

pc1 = gridExtra::grid.arrange(p1, p2, nrow=2, heights=c(1.2, 0.9))
pc2 = ggpubr::as_ggplot(pc1)  +
  cowplot::draw_plot_label(label = c("a", "b"),
                           fontface = "bold", size = 25,
                           x = c(0.01, 0.01), y = c(0.99, 0.47))
ggsave(pc2, file="./Figure3_tickmodels_v2.jpeg", units="in", scale=0.9, width=15, height=6, dpi=600)




# ================ composite plot ==================

# pc = gridExtra::grid.arrange(uga_plot1, uga_plot2, hruf_af, htru_af, ncol=2, heights=c(1, 0.9))
# pc2 = ggpubr::as_ggplot(pc)  +
#   cowplot::draw_plot_label(label = c("a", "b", "c", "d"),
#                            fontface = "bold", size = 24,
#                            x = c(0.02, 0.02, 0.52, 0.52), y = c(0.99, 0.48, 0.99, 0.48))
# ggsave(pc2, file="./Figure3_tickmodels.jpeg", units="in", scale=0.9, width=8, height=8, dpi=600)


# plot v2
pc1 = gridExtra::grid.arrange(hruf_af, uga_plot1, htru_af, uga_plot2, ncol=2, widths=c(0.9, 1))
pc2 = ggpubr::as_ggplot(pc1)  +
  cowplot::draw_plot_label(label = c("a", "c", "b", "d"),
                           fontface = "bold", size = 24,
                           x = c(0.03, 0.03, 0.54, 0.54), y = c(0.97, 0.47, 0.97, 0.47))
ggsave(pc2, file="./Figure3_tickmodels_v2.jpeg", units="in", scale=0.9, width=8, height=8, dpi=600)



