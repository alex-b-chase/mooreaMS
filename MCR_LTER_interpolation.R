rm(list=ls())
require(sf)
require(tidyverse)
require(mapview)
library(kriging)
setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/publicdata/")

nutb = read_csv("MCR_LTER_Adam_EcolAp_Nut_boundary_20200331.csv")
nuts = read_csv("burkepile_dataPNAS/MCR_LTER_Burkepile_Adams_2016_Coral_Bleaching_Survey_20191023.csv")
extrapoints = read_csv("extra_boundaryGPS.csv")

samples = read.table("../moorea2020/16S_analysis/MO18-metadata.txt", header = T, sep = '\t', comment.char = '')
samples2 <- samples[samples$BarcodeSequence != "categorical",]
nuts_samp = samples2 %>% st_as_sf(coords = c("longitude", "latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
nuts_samp_fil = nuts_samp %>% select(siteID, geometry) %>% distinct() %>% na.omit()
nuts_samp_fil_rpj = nuts_samp_fil %>% st_transform(crs = 2976)


nuts_ext = extrapoints %>% st_as_sf(coords = c("long", "lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")
nuts_ext_fil_rpj = nuts_ext %>% st_transform(crs = 2976)

nuts_sf = nuts %>% st_as_sf(coords = c("Longitude", "Latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84")

nuts_sf_fil = nuts_sf %>% select(avgTotN, geometry) %>% distinct() %>% na.omit()
nuts_sf_fil_rpj = nuts_sf_fil %>% st_transform(crs = 2976)

nutb = nutb %>% st_as_sf(coords = c("X", "Y"), crs = 2976)

write_csv(nuts_ext_fil_rpj, file = "temp.csv")
mapview(nutb) + mapview(nuts_sf_fil_rpj)


x = (nuts_sf_fil_rpj %>% st_coordinates())[, 1]
y = (nuts_sf_fil_rpj %>% st_coordinates())[, 2]
z = nuts_sf_fil_rpj$avgTotN

bound = list(nutb %>% st_coordinates() %>% as.data.frame())

krig1 <- kriging(x=x, y=y, response=z, pixels=1000, polygons=bound) 

##image(krig1, xlim = extendrange(x), ylim = extendrange(y), col = rev(rainbow(12)))

mat = krig1$map
krigLLdf <- as.data.frame(mat)

p1 <- ggplot() +
  theme_bw() +
  geom_raster(data=krigLLdf, aes(x=x, y=y, fill=pred)) +
  geom_point(aes(x=202359.6, y=8065270), size = 3) +
  geom_point(aes(x=203557.2, y=8065515), size = 3) +
  geom_point(aes(x=204394.7, y=8065669), size = 3) +
  geom_point(aes(x=200898, y=8064911), size = 3) +
  geom_point(aes(x=204915.4, y=8065676), size = 3) +
  geom_point(aes(x=202492.5, y=8064906), size = 3) +
  geom_point(aes(x=203435.8, y=8065090), size = 3) +
  geom_point(aes(x=204492.8, y=8065501), size = 3) +
  scale_fill_gradientn(colours = rev(rainbow(10))) +
  coord_sf()
p1
mat2 = mat  %>% st_as_sf(coords = c("x", "y"), crs = 2976)
totalpred <- st_join(nuts_samp_fil_rpj, mat2, join=st_nearest_feature)
totalpred

### now do heat stress
heat_sf_fil = nuts_sf %>% select(cumstress, geometry) %>% distinct() %>% na.omit()
heat_sf_fil_rpj = heat_sf_fil %>% st_transform(crs = 2976)

x = (heat_sf_fil_rpj %>% st_coordinates())[, 1]
y = (heat_sf_fil_rpj %>% st_coordinates())[, 2]
z = heat_sf_fil_rpj$cumstress

krig2 <- kriging(x=x, y=y, response=z, pixels=1000, polygons=bound) 

mat = krig2$map
krigLLdf <- as.data.frame(mat)

p2 <- ggplot() +
  theme_bw() +
  geom_raster(data=krigLLdf, aes(x=x, y=y, fill=pred)) +
  geom_point(aes(x=202359.6, y=8065270), size = 3) +
  geom_point(aes(x=203557.2, y=8065515), size = 3) +
  geom_point(aes(x=204394.7, y=8065669), size = 3) +
  geom_point(aes(x=200898, y=8064911), size = 3) +
  geom_point(aes(x=204915.4, y=8065676), size = 3) +
  geom_point(aes(x=202492.5, y=8064906), size = 3) +
  geom_point(aes(x=203435.8, y=8065090), size = 3) +
  geom_point(aes(x=204492.8, y=8065501), size = 3) +
  scale_fill_gradientn(colours = rev(rainbow(10))) +
  coord_sf()
p2

library(ggmap)
library(ggrepel)
library(gridExtra)
register_google(key = "AIzaSyDT8zC9_PlR5ROJnv8vGtDp7aZQe5QnJuI")
map <- get_googlemap(c(-137.01, 20), zoom = 4,
                     color = "bw")
ggmap(map)
bw_map3 <- get_googlemap(c(-149.84, -17.53), zoom = 12,
                         color = "bw", maptype = "satellite",
                         style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature:administrative|visibility:off")

p3 <- ggmap(bw_map3)

pdf("moorea_map.pdf", height = 4, width =10)
grid.arrange(p1,p2,p3, nrow=1)
dev.off()

mat2 = mat  %>% st_as_sf(coords = c("x", "y"), crs = 2976)
totalpred <- st_join(nuts_samp_fil_rpj, mat2, join=st_nearest_feature)
totalpred

### now do bleaching - need to average bleaching by location (geometry) 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=T,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

nuts2 <- na.omit(nuts)
aggregate(nuts2$avgTotN, by = list(nuts2$LTER_site), mean)




### this is for C:N ratios - don't know what to make of it...
#### get abiotic factors from Moorea LTER site
library(reshape2)
abiotics = read_csv("MCR_LTER_Macroalgal_CHN_2005_to_2020_20220322.csv")
abiotics2 <- abiotics[, c("Site", "Year", "CN_ratio")]
abiotics2$Year <- as.character(abiotics2$Year)
abioticsm <- melt(abiotics2, value.name = "CN_ratio")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summtot <- summarySE(abioticsm, measurevar = "CN_ratio", groupvars = c("Site", "Year"))

ggplot(summtot, aes(x = Year, y = CN_ratio, group = Site, color = Site)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = CN_ratio - se, ymax = CN_ratio + se),
                width = 0,                    # Width of the error bars
                position = position_dodge(.1)) +
  geom_smooth(se = F) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) 

summtot2 <- summarySE(abioticsm, measurevar = "CN_ratio", groupvars = c("Site"))

cn2018 <- abioticsm[abioticsm$Year == "2018",]
summtot2018 <- summarySE(cn2018, measurevar = "CN_ratio", groupvars = c("Site"))

ggplot(summtot2, aes(x = Site, y = CN_ratio, fill = "")) +
  stat_summary(fun = "mean", geom = "bar", aes(colour="black")) +
  geom_errorbar(aes(ymin = CN_ratio - se, ymax = CN_ratio + se),
                width = 0,                    # Width of the error bars
                position = position_dodge(1)) +
  geom_point(data = summtot2018, aes(x = Site, y = CN_ratio), fill = "red", size = 4, shape = 23) +
  geom_errorbar(data = summtot2018, aes(ymin = CN_ratio - se, ymax = CN_ratio + se), color = "red",
                width = 0,                    # Width of the error bars
                position = position_dodge(.1)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "C:N Ratio", x = "") +
  scale_color_manual(values = c("black"), guide = F) +
  scale_fill_manual(values = c("lightgray"), guide = F)

