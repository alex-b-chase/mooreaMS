# network analysis for population structure
setwd('/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/bigscape/BGCnetwork/') 

library(tidyverse)

###############################################################
###############################################################
# BigScape calls a GCF at c=0.6 
###############################################################
###############################################################

gephiGCF <- read.table('../bigscape_out/network_files/2022-09-21_14-59-54_hybrids_glocal/mix/mix_clustering_c0.80.tsv', sep = '\t', header = F)
colnames(gephiGCF) <- c("genome", "gcfnum")
gephiGCF$GCFFAM <- paste("GCF", gephiGCF$gcfnum, sep = "")
gephiGCF$gcfnum <- NULL

library(splitstackshape)
gephiGCFtemp <- cSplit(gephiGCF, "genome", ".")

library(reshape2)

GCFABD <- dcast(gephiGCFtemp, genome_1 ~ GCFFAM)
rownames(GCFABD) <- GCFABD$genome_1
GCFABD$genome_1 <- NULL

# these should match the total BGC counts minus singletons
rowSums(GCFABD)
colSums(GCFABD)

write.table(GCFABD, "bigscape_c0.30.tsv", sep = '\t', row.names = T, quote = F)

library(gplots)

pdf('GCF-abundance.pdf', height = 15, width = 20)

my_palette <- colorRampPalette(c("white", "darkgray", "black"))(n = 200)
heatmap.2(as.matrix(GCFABD),
          main = "PA of BGCs", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace = "row",
          tracecol = "lightgray",
          margins =c(2,1),     # widens margins around plot
          col=my_palette       # use on color palette defined earlier
)  

dev.off()

site6 <- c("MO_167", "MO_170", "MO_180", "MO_184", "MO_188", "MO_192")
GCFABD2 <- GCFABD[ ! rownames(GCFABD) %in% site6, ]
GCFABD2 <- GCFABD2[, which(colSums(GCFABD2) != 0)]

metadata <- read.table("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/mgenomeID.txt", header = T, sep = '\t', comment.char = "", row.names = 1)
tokeep <- rownames(GCFABD2)
metadataGCF <- metadata[which(rownames(metadata) %in% tokeep), ]

total.info <- merge(metadataGCF, GCFABD2, by = 0, all = F)
row.names(total.info) <- total.info$Row.names
total.info$site <- factor(total.info$site)
total.info$reef <- factor(total.info$reef)
total.info$sitecolor <- factor(total.info$sitecolor)

library(vegan)
library(ggplot2)

sol <- metaMDS(GCFABD2, distance = "bray", k = 2, trymax = 5000)

#Make a new data frame, and put country, latrine, and depth information there, to be useful for coloring, and shape of points
NMDS = data.frame(x = sol$point[,1], y = sol$point[,2], 
                  site = as.factor(total.info$site), 
                  reef = as.factor(total.info$reef),
                  sitecolor = as.factor(total.info$sitecolor))

plot.new()
ord <- ordiellipse(sol, total.info$site, display = "sites", kind ="sd", label = T)
dev.off()


veganCovEllipse <-
  function(cov, center = c(0,0), scale = 1, npoints = 100){
    ## Basically taken from the 'car' package: The Cirlce
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    ## scale, center and cov must be calculated separately
    Q <- chol(cov, pivot = TRUE)
    ## pivot takes care of cases when points are on a line
    o <- attr(Q, "pivot")
    t(center + scale * t(Circle %*% Q[,o]))
  }

df_ell <- data.frame()
for(g in levels(NMDS$site)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$site==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center)))
                                ,site=g))
}

NMDS.mean=aggregate(NMDS[,1:2], list(site = NMDS$site), mean)
staxcolors <- setNames(as.character(NMDS$sitecolor), NMDS$site)

pdf('mdsbray-GCFbigscape.pdf', height = 15, width = 20)

ggplot(data = NMDS, aes(x, y, color = site)) + 
  geom_point(size = 10, stroke = 3, alpha = 0.8) +
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, color = site), size = 3) +
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y, label = NMDS.mean$site, size = 8) + 
  theme_bw() +
  scale_colour_manual(values = staxcolors) +
  scale_fill_manual(values = staxcolors)

dev.off()

# distance matrix for BGC across
m <- as.matrix(dist(GCFABD2))

library(ecodist)
library(vegan)

dist.mat <- vegdist(GCFABD2, method = "bray", diag = T)

clust.res <- hclust(dist.mat)
plot(clust.res)

#### change from Mantel test to PERMANOVA
perm.log <- adonis2(dist.mat ~ reef * site, data=metadataGCF, permutations = 999, method = "bray")
perm.log

##            Df SumOfSqs      R2      F Pr(>F)    
##  reef      1   0.6983 0.06093 2.0350  0.002 ** 
##  site      5   2.8701 0.25043 1.6728  0.001 ***
##  Residual 23   7.8923 0.68864                  
##  Total    29  11.4607 1.00000    

write.csv(as.matrix(dist.mat),'braycurt_BGCdist_nosite6.csv',quote=F)

