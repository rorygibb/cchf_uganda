
# ================== Panel plot with respective areas of high and low suitability for Hyalomma and CCHF ======================

# for MS SI

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

shpz_simpl = shpz %>% 
  rmapshaper::ms_simplify(keep = 0.01) %>%
  smoothr::smooth(method = "chaikin")

# read in predicted Hyalomma maps
hya = stack(list.files("./output/figures/manuscript/hyalomma_suitability_layers/", full.names=TRUE, pattern="meanpred.tif")) 
hya = raster::aggregate(hya, fact=3, fun="mean")

# read in and calculate thresholds
t1 = read.csv("./output/figures/manuscript/ensemble_metadata/_Hyalomma truncatum_563139_ensemblemetadata.csv") %>% dplyr::summarise(mean(threshold))
t2 = read.csv("./output/figures/manuscript/ensemble_metadata/_Hyalomma rufipes_824893_ensemblemetadata.csv") %>% dplyr::summarise(mean(threshold))
t3 = read.csv("./output/figures/manuscript/ensemble_metadata/_Hyalomma impeltatum_686854_ensemblemetadata.csv") %>% dplyr::summarise(mean(threshold))
t4 = read.csv("./output/figures/manuscript/ensemble_metadata/_Hyalomma dromedarii_111728_ensemblemetadata.csv") %>% dplyr::summarise(mean(threshold))

# threshold maps
h1 = hya[[ 1 ]] >= as.numeric(t1)
h2 = hya[[ 2 ]] >= as.numeric(t2)
h3 = hya[[ 3 ]] >= as.numeric(t3)
h4 = hya[[ 4 ]] >= as.numeric(t4)

# combine into hyalomma layer
hya_comb = h1 + h2 + h3 + h4
hya_comb = hya_comb > 0 
plot(hya_comb)

# CCHF presence/absence
cchf = raster("./output/figures/manuscript/cchf_suitability_layers/CCHF_meanpred.tif")
cchf = raster::aggregate(cchf, fact=3, fun="mean")
c1 = read.csv("./output/figures/manuscript/ensemble_metadata/_cchf_10km_563371_ensemblemetadata.csv") %>% dplyr::summarise(mean(threshold))
cchf = cchf >= as.numeric(c1)

# create rasters of low low etc
h2 = crop(hya_comb, cchf)

low_low = (h2 + cchf) == 0
high_high = (h2 + cchf) == 2

v1 = values(h2); v2 = values(cchf)
hh_lc = raster(cchf)
values(hh_lc) = v1 == 1 & v2 == 0

v1 = values(h2); v2 = values(cchf)
lh_hc = raster(cchf)
values(lh_hc) = v1 == 0 & v2 == 1


# stack and save
ss = stack(low_low, high_high, hh_lc, lh_hc)
names(ss) = c("LowHyalomma_LowCCHF", "HighHyalomma_HighCCHF", "HighHyalomma_LowCCHF", "LowHyalomma_HighCCHF")
dd = as.data.frame(ss, xy=TRUE) %>%
  reshape2::melt(id.vars=1:2) %>%
  dplyr::mutate(HyalommaSuitability = ifelse(grepl("LowHyalomma", variable), "Low Hyalomma", "High Hyalomma"),
                CCHFSuitability = ifelse(grepl("LowCCHF", variable), "Low CCHF", "High CCHF"),
                HyalommaSuitability = factor(HyalommaSuitability, levels=c("Low Hyalomma", "High Hyalomma"), ordered=TRUE),
                CCHFSuitability = factor(CCHFSuitability, levels=c("Low CCHF", "High CCHF"), ordered=TRUE))

# plot
bivar_facet = dd %>%
  ggplot() +
  geom_raster(aes(x=x, y=y, fill=value)) + 
  facet_grid(HyalommaSuitability ~ CCHFSuitability, switch = "y") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14, angle=90)) +
  scale_fill_viridis_d(begin=0, end=0.85, na.value="white") +
  geom_sf(data=shpz_simpl, fill=NA, color="black", size=0.3)
ggsave(bivar_facet, file="./SIFIgure_BivarPlot.jpeg", units="in", height=8, width=8, dpi=600)


# =================== alternative approach: do it properly ===================

setwd("D:/ResearchProjects/202002_ucl_congocrimean/analysis/")

calculateConsensusRaster = function(filepath){
  
  ss = stack(list.files(filepath, pattern=".tif", full.names=TRUE))
  
  # get thresholds
  threshs = unlist(lapply(strsplit(names(ss), "_"), "[", 8))
  
  # threshold them
  res_stack = stack()
  for(i in 1:nlayers(ss)){
    rr = raster::aggregate(ss[[i]], fact=3, fun="mean")
    rr = rr >= as.numeric(threshs[i])
    res_stack = stack(res_stack, rr)
  }
  
  cr = modal(res_stack)
  return(cr)
}

# save consensus rasters
cchf_cons = calculateConsensusRaster("./output/model_outputs/cchf_brt/fitted_models/cchf_10km_BRT_2020-07-07_563371/geo_preds/")
h1_cons = calculateConsensusRaster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma truncatum_BRT_2020-07-05_563139_envonly/geo_preds/")
h2_cons = calculateConsensusRaster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma impeltatum_BRT_2020-07-01_686854_wlepus/geo_preds/")
h3_cons = calculateConsensusRaster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma rufipes_BRT_2020-07-05_824893_envonly/geo_preds/")
h4_cons = calculateConsensusRaster("./output/model_outputs/climate_ensemble/fitted_models/Hyalomma dromedarii_BRT_2020-06-24_111728_wolepus/geo_preds/")

# save
writeRaster(cchf_cons, "C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/output/figures/manuscript/consensus_suitable/cchf_consensus.tif", format="GTiff")
writeRaster(h1_cons, "C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/output/figures/manuscript/consensus_suitable/htruncatum_consensus.tif", format="GTiff")
writeRaster(h2_cons, "C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/output/figures/manuscript/consensus_suitable/himpeltatum_consensus.tif", format="GTiff")
writeRaster(h3_cons, "C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/output/figures/manuscript/consensus_suitable/hrufipes_consensus.tif", format="GTiff")
writeRaster(h4_cons, "C:/Users/roryj/Documents/PhD/202002_ucl_crimeancongo/analysis/output/figures/manuscript/consensus_suitable/hdromedarii_consensus.tif", format="GTiff")


hya_comb = h1_cons + h2_cons + h3_cons + h4_cons
hya_comb = hya_comb > 0 
plot(hya_comb)

cchf = cchf_cons

h2 = crop(hya_comb, cchf)

low_low = (h2 + cchf) == 0
high_high = (h2 + cchf) == 2

v1 = values(h2); v2 = values(cchf)
hh_lc = raster(cchf)
values(hh_lc) = v1 == 1 & v2 == 0

v1 = values(h2); v2 = values(cchf)
lh_hc = raster(cchf)
values(lh_hc) = v1 == 0 & v2 == 1


# stack and save
ss = stack(low_low, high_high, hh_lc, lh_hc)
names(ss) = c("LowHyalomma_LowCCHF", "HighHyalomma_HighCCHF", "HighHyalomma_LowCCHF", "LowHyalomma_HighCCHF")
dd = as.data.frame(ss, xy=TRUE) %>%
  reshape2::melt(id.vars=1:2) %>%
  dplyr::mutate(HyalommaSuitability = ifelse(grepl("LowHyalomma", variable), "Low Hyalomma", "High Hyalomma"),
                CCHFSuitability = ifelse(grepl("LowCCHF", variable), "Low CCHF", "High CCHF"),
                HyalommaSuitability = factor(HyalommaSuitability, levels=c("Low Hyalomma", "High Hyalomma"), ordered=TRUE),
                CCHFSuitability = factor(CCHFSuitability, levels=c("Low CCHF", "High CCHF"), ordered=TRUE))

# plot
bivar_facet = dd %>%
  ggplot() +
  geom_raster(aes(x=x, y=y, fill=value)) + 
  facet_grid(HyalommaSuitability ~ CCHFSuitability, switch = "y") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14, angle=90)) +
  scale_fill_viridis_d(begin=0, end=0.85, na.value="white") +
  geom_sf(data=shpz_simpl, fill=NA, color="black", size=0.3)
ggsave(bivar_facet, file="./SIFIgure_BivarPlot.jpeg", units="in", height=8, width=8, dpi=600)
