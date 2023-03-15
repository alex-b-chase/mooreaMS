rm(list=ls())
setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/metabolomics")

library(vegan)
library(ggplot2)
library(ecodist)
library(RFmarkerDetector) # https://rdrr.io/cran/RFmarkerDetector/
# RFmarkerDetector: Multivariate Analysis of Metabolomics Data using Random Forests

totalMS <- read.table(file = "allMS1data_sedimentV2.txt", header = T, sep = '\t', row.names = 1)
totalMS1 <- totalMS[,!grepl("duration", colnames(totalMS))]
totalMS2 <- totalMS1[,!grepl("RT", colnames(totalMS1))]
totalMS3 <- totalMS2[,!grepl("m.z", colnames(totalMS2))]
totalMS4 <- totalMS3[,!grepl("Solvent", colnames(totalMS3))]
totalMS2 <- totalMS4[,!grepl("area", colnames(totalMS4))]
colnames(totalMS2) = gsub(".Peak.height", "", colnames(totalMS))
baseMSdf <- totalMS2

# read in metadata
metadata <- read.table("../Metabolomics_Metadata_sediment.txt", header = T, sep = '\t', row.names = 1, comment.char = "")
staxcolors <- setNames(as.character(metadata$sitecolor), metadata$ATTRIBUTE_Site)

### remove singleton/doubleton masses from the feature table
#totalMS <- totalMS2[!rowSums(totalMS2 == 0) >= 124, , drop = FALSE]

totalMS <- totalMS2
totalMS[totalMS <= 1E3] <- NA
# remove all rows with NA - meaning only found in media
totalMS2 <- totalMS[rowSums(is.na(totalMS)) != ncol(totalMS), ]

my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
totalMS2[is.na(totalMS2)] <- my.min(totalMS2) / 10

tflex <- t(totalMS2)
logtflex <- log(tflex, 2)
paretotflex <- paretoscale(logtflex)
scaletflex <- scale(logtflex)

dist.mat <- vegdist(paretotflex, method = "euclidean", diag = T)
clust.res <- hclust(dist.mat)
plot(clust.res)

hist(t(logtflex))
hist(t(scaletflex))
hist(t(paretotflex))

################################################
################################################
#### Principal Component Analysis (PCoA)
################################################
################################################

# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
pcatflex2 <- merge(metadata, paretotflex, by = 0)
rownames(pcatflex2) <- pcatflex2$Row.names

pcatflex <- pcatflex2[pcatflex2$ATTRIBUTE_Coordinates != "seawater",]

#pcatflex$ATTRIBUTE_Site <- droplevels(pcatflex$ATTRIBUTE_Site)

# PERMANOVA stats analysis
perm.pareto <- adonis2(dist.mat ~ ATTRIBUTE_Site, data=pcatflex, permutations = 999, method = "euclidean")
perm.pareto


# PCoA analysis for MS data
pca <- prcomp(t(paretotflex), center = F, scale = F)
pcaresults <- summary(pca)
pcaresults$importance[3,1:3] 

pcadata <- as.data.frame(pcaresults$rotation)
pcadataPC12 <- pcadata[, c(1:3)]
 
pcameta <- merge(metadata, pcadataPC12, by = 0, all = F)
staxcolors <- setNames(as.character(pcameta$sitecolor), pcameta$ATTRIBUTE_Site)

pdf("allMS1data.pdf", height = 10, width = 12)
ggplot(pcameta, aes(PC1, PC2)) +
  geom_point(aes(color = ATTRIBUTE_Site), size = 4) +
  geom_text(aes(label = ATTRIBUTE_SampleID)) +
  theme_bw() +
  scale_colour_manual(values = staxcolors)

dev.off()

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




summPC1 <- summarySE(pcameta, measurevar = "PC1", groupvars = c("ATTRIBUTE_Site", "sitecolor"))
summPC1$PC1sd <- summPC1$sd
summPC2 <- summarySE(pcameta, measurevar = "PC2", groupvars = c("ATTRIBUTE_Site", "sitecolor"))
summPC2$PC2sd <- summPC2$sd

totalPC12temp <- merge(summPC1, summPC2, by = c("ATTRIBUTE_Site", "sitecolor"))
totalPC12 <- totalPC12temp[, c("ATTRIBUTE_Site", "PC1", "PC2", "PC1sd", "PC2sd", "sitecolor")]
staxcolors <- setNames(as.character(totalPC12$sitecolor), totalPC12$ATTRIBUTE_Site)

pdf("MS1data-avgsite.pdf", height = 10, width = 12)
ggplot(totalPC12, aes(PC1, PC2)) +
  geom_errorbarh(aes(xmin = PC1 - PC1sd, xmax = PC1 + PC1sd), height = 0) +
  geom_errorbar(aes(ymin = PC2 - PC2sd, ymax = PC2 + PC2sd), width = 0) +
  geom_point(aes(color = ATTRIBUTE_Site), size = 6, shape = 15) +
  theme_bw() +
  scale_colour_manual(values = staxcolors) +
  scale_fill_manual(values = staxcolors)

dev.off()

################################################
################################################
#### remove site 6 from analysis and re-run
################################################
################################################

site6 <- c("MO_167.mzXML", "MO_170.mzXML", "MO_180.mzXML", "MO_184.mzXML",
           "MO_188.mzXML", "MO_192.mzXML")

totalMS <- baseMSdf[ , !(names(baseMSdf) %in% site6)]

totalMS[totalMS <= 1E3] <- NA
## need to remove features only found in site 6
totalMS2 <- totalMS[rowSums(is.na(totalMS)) != ncol(totalMS), ]

my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
totalMS2[is.na(totalMS2)] <- my.min(totalMS2) / 10

tflex <- t(totalMS2)
logtflex <- log(tflex, 2)
paretotflex <- paretoscale(logtflex)
scaletflex <- scale(logtflex)

hist(logtflex)
hist(paretotflex)
hist(scaletflex)


# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
pcatflex2 <- merge(metadata, paretotflex, by = 0)
rownames(pcatflex2) <- pcatflex2$Row.names

# PERMANOVA stats analysis
dist.mat <- vegdist(paretotflex, method = "euclidean", diag = T)
perm.pareto <- adonis2(dist.mat ~ ATTRIBUTE_Reef*ATTRIBUTE_Site, data=pcatflex2, permutations = 999, method = "euclidean")
perm.pareto
write.csv(as.matrix(dist.mat), file = "MS1distmat_sedimentV2.csv", quote = F)
write.csv(as.matrix(paretotflex), file = "MS1paretoscale_sedimentV2.csv", quote = F)
write.csv(as.matrix(totalMS2), file = "MS1nosite6_sedimentV2.csv", quote = F)



# PCoA analysis for MS data
pca <- prcomp(t(paretotflex), center = F, scale = F)
pcaresults <- summary(pca)
pcaresults$importance[3,1:3] 

pcadata <- as.data.frame(pcaresults$rotation)
pcadataPC12 <- pcadata[, c(1:3)]

pcameta <- merge(metadata, pcadataPC12, by = 0, all = F)
staxcolors <- setNames(as.character(pcameta$sitecolor), pcameta$ATTRIBUTE_Site)

pdf("allMS1data_nosite6.pdf", height = 10, width = 12)
ggplot(pcameta, aes(PC1, PC2)) +
  geom_point(aes(color = ATTRIBUTE_Site), size = 4) +
  geom_text(aes(label = ATTRIBUTE_SampleID)) +
  theme_bw() +
  scale_colour_manual(values = staxcolors)

dev.off()

################################################
### plot PCoA distances
################################################

library("FactoMineR")
library("factoextra")


site.pca <- PCA(pcatflex2[,c(9:length(colnames(pcatflex2)))], graph = FALSE)
fviz_contrib(site.pca, choice = "var", axes = 1:2, top = 250)

## alexB curated the features and we agreed to remove some that were "not real"
## remove features 392, 625, 652 (PEG contaminants)
## remove features 10 and 1988 ( same molecule with RT slightly off )

pdf("PCA-biplot_nofeatures.pdf", height = 10, width = 8 )
fviz_pca_biplot(site.pca, geom.ind = "point",
                fill.ind = pcatflex2$ATTRIBUTE_Site, col.ind = "black",
                pointshape = 24, pointsize = 6,
                palette = c("#597343", "#8435BB", "#36CDFF", "#DC61E5", 
                            "#DAC75D", "#E15314", "#FF7983", "#d1db4b"), 
                addEllipses = TRUE, ellipse.level = 0.75, 
                select.var = list(name = c("feature772")),
                label = "var", 
                repel = TRUE, mean.point = FALSE,
                #alpha.var = "contrib", 
                col.var = "contrib",
                legend.title = "MS Features") 
dev.off()

pdf("PCA-biplot_wfeatures.pdf", height = 14, width = 10.88)
fviz_pca_biplot(site.pca, geom.ind = "point",
                fill.ind = pcatflex2$ATTRIBUTE_Site, col.ind = "black",
                pointshape = 24, pointsize = 10,
                palette = c("#597343", "#8435BB", "#36CDFF", "#DC61E5", 
                            "#DAC75D", "#E15314", "#FF7983", "#d1db4b"), 
                addEllipses = TRUE, ellipse.level = 0.75, 
                select.var = list(name = c("feature672", "feature709", "feature248", "feature788", 
                                           "feature615", "feature1984", "feature801", "feature489", 
                                           "feature772", "feature815", "feature460", "feature1913", 
                                           "feature1620", "feature183", "feature702", "feature370")),
                label = "var", 
                repel = TRUE, mean.point = FALSE,
                #alpha.var = "contrib", 
                col.var = "contrib",
                legend.title = "MS Features") 
dev.off()


featurelist <- read.table("featurelist_sedV2.txt", header = T, sep = '\t', row.names = 1)
## possibly brominated compound == mz=645.207 rt=21.4 and 21.65 (isomers)
featurelist[grep("^645.", featurelist$mz),]
## feature902 and feature1446
## cyanobacterial peptide == mz=759.4311 rt=16.57
featurelist[grep("^759.", featurelist$mz),]
featurelist$featureID <- rownames(featurelist)
featurelist[featurelist$featureID == "feature616",]
## feature1382

# subset the top 50 features along axes 1 and 2
contribFEAT <- as.data.frame(site.pca$var$contrib)
contribFEAT$PCA12 <- rowMeans(contribFEAT[1:2], na.rm=TRUE)
contribFEAT$featureID <- rownames(contribFEAT)
contribFEAT[contribFEAT$featureID == "feature616",]
## definitely a dropoff at 12% - subset by that instead
topFEAT <- contribFEAT[contribFEAT$PCA12 >= 0.12,]

barplot(sort(contribFEAT$PCA12), ylim = c(0, max(contribFEAT$PCA12)),
       xlim = c(0, NROW(contribFEAT$PCA12)), col = "Blue", ylab = "Contribution", xlab = "Feature")

topfeatures <- merge(topFEAT, featurelist, by = 0)
topfeatures2 <- topfeatures[,!grepl("Dim.", colnames(topfeatures))]
write.csv(topfeatures2, file = "topfeaturesPCOA.csv", quote = F)

################################################
################################################
### convert MS data to P/A
################################################
################################################

site6 <- c("MO_167.mzXML", "MO_170.mzXML", "MO_180.mzXML", "MO_184.mzXML",
           "MO_188.mzXML", "MO_192.mzXML")

totalMS <- baseMSdf[ , !(names(baseMSdf) %in% site6)]

totalMS[totalMS <= 1E3] <- NA
## need to remove features only found in site 6
totalMS2 <- totalMS[rowSums(is.na(totalMS)) != ncol(totalMS), ]

totalMS2[is.na(totalMS2)] <- 0

tflexJ <- t(totalMS2)

### nMDS
Jdist.mat <- vegdist(tflexJ, method = "jaccard", diag = T, binary = T)
sol <- metaMDS(Jdist.mat, k = 3, trymax = 500, autotransform = F)
plot(sol$points, col = metadata2$sitecolor)

tokeep <- rownames(sol$points)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]

NMDS = data.frame(MDS1 = sol$points[,1], MDS2=sol$points[,2], group = metadata2$ATTRIBUTE_Site)
NMDS.mean = aggregate(NMDS[,1:2], list(group = NMDS$group), mean)
NMDS$group <- as.factor(NMDS$group)


#To make a more sophisticated plot, we will merge the stress scores with the metadata.
tokeep <- rownames(sol$points)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]
nmds_plus_metadata <- merge(NMDS, metadata, by = 0)

samplecolors <- setNames(as.character(nmds_plus_metadata$sitecolor), nmds_plus_metadata$ATTRIBUTE_Site)

en = envfit(sol, tflexJ, permutations = 999, na.rm = TRUE)
en2 = envfit(sol, metadata2, permutations = 999, na.rm = TRUE)

## extract pvalues and R2 values from feature table
pval <- as.data.frame(en$vectors$pvals)
pval$featureID <- rownames(pval)
colnames(pval) <- c("pvalue", "featureID")
sigpval <- pval[pval$pvalue < 0.05, ]

r2 <- as.data.frame(en$vectors$r)
r2$featureID <- rownames(r2)
colnames(r2) <- c("r2", "featureID")
tokeeppvals <- rownames(sigpval)
r2_pval <- r2[which(rownames(r2) %in% tokeeppvals), ]
totalvector <- merge(sigpval, r2_pval, by = "featureID")

### get the coordinates for the top ones
en_coord_cont <- as.data.frame(scores(en, display = "vectors"))
en_coord_cont <- cbind(en_coord_cont, Species = rownames(en_coord_cont))

### now, only subset the top 20 features
topfeatures <- totalvector[with(totalvector,order(-r2)),]
topfeatures <- topfeatures[1:20,]
top20features <- topfeatures$featureID

en_coord_cont <- en_coord_cont[which(rownames(en_coord_cont) %in% top20features), ]

barplot(rev(sort(totalvector$r2)), ylim = c(0, max(totalvector$r2)),
        xlim = c(0, NROW(totalvector$r2)), col = "Blue", ylab = "Contribution", xlab = "Feature")

## extract centroids for back and fringe reef
en_coord_cat = as.data.frame(scores(en2, "factors"))
fringestatus <- c("ATTRIBUTE_ReefBack_Reef", "ATTRIBUTE_ReefFringing_Reef")
en_coord_cat <- en_coord_cat[which(rownames(en_coord_cat) %in% fringestatus), ]

pdf("nMDS_jaccardMS.pdf", height = 12, width = 14)

ggplot(nmds_plus_metadata, aes(MDS1, MDS2)) +
  stat_ellipse(geom = "polygon",
               aes(fill = ATTRIBUTE_Site), 
               alpha = 0.25, level = 0.85) +
  coord_fixed() +
  geom_point(aes(color = ATTRIBUTE_Site), size = 6) +
  geom_segment(data = en_coord_cont,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  geom_text(aes(label = ATTRIBUTE_SampleID)) +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  scale_fill_manual(values = samplecolors)

dev.off()

# Use scale = TRUE if your variables are on different scales (e.g. for abiotic variables).
pcatflex2 <- merge(metadata, tflexJ, by = 0)
rownames(pcatflex2) <- pcatflex2$Row.names

# PERMANOVA stats analysis
Jdist.mat <- vegdist(tflexJ, method = "jaccard", diag = T, binary = T)
perm.pareto <- adonis2(Jdist.mat ~ ATTRIBUTE_Site, data=pcatflex2, permutations = 999, method = "euclidean")
perm.pareto
write.csv(as.matrix(Jdist.mat), file = "MS1jaccdist_sedimentV2.csv", quote = F)

clust.res <- hclust(Jdist.mat)
plot(clust.res)

library(ape)

# Principal coordinate analysis and simple ordination plot
PCOAJ <- pcoa(Jdist.mat)
# plot the eigenvalues and interpret
barplot(PCOAJ$values$Relative_eig[1:10])
PCOAJ$values$Relative_eig[1:3]
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Plot your results
biplot(PCOAJ)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
# However, we could work around this problem like this:
bipcoa <- biplot.pcoa(PCOAJ, tflexJ)


PCOAaxesJ <- PCOAJ$vectors[,c(1,2)]
PCOAmetaJ <- metadata[, c(1,2,3)]

PCOAplotJ <- merge(PCOAaxesJ, PCOAmetaJ, by = "row.names")

# should look the same as the above one, just without the nice graphing stuff
ggplot(data = PCOAplotJ, aes(Axis.1, Axis.2)) +
  stat_ellipse(geom = "polygon",
               aes(fill = ATTRIBUTE_Site), 
               alpha = 0.25, level = 0.85) +
  geom_point(aes(color = ATTRIBUTE_Site), size = 6) +
  scale_colour_manual(values = staxcolors) +
  scale_fill_manual(values = samplecolors) +
  theme_bw()




################################################
################################################
### compare MS data to geographic distances
################################################
################################################
geodist <- as.matrix(read.table("../../16S_analysis/geodist.txt", header = T, row.names = 1))
library(reshape2)

geodistm <- melt(geodist)
colnames(geodistm) <- c("sample1", "sample2", "geodist")

MSpcoa <- melt(as.matrix(Jdist.mat))
colnames(MSpcoa) <- c("sample1", "sample2", "Jac")

MSpcoa$sample1 = gsub(".mzXML", "", MSpcoa$sample1)
MSpcoa$sample1 = gsub("MO_", "MO18_", MSpcoa$sample1)
MSpcoa$sample2 = gsub(".mzXML", "", MSpcoa$sample2)
MSpcoa$sample2 = gsub("MO_", "MO18_", MSpcoa$sample2)

MSpcoa$simDist <- 0 + MSpcoa$Jac
MSpcoa$Jac <- NULL
MSpcoa$distmet <- "JacMS"

jac16S <- as.matrix(read.table("../../16S_analysis/qiime2_nosite6/coremetrics/jaccard_dist/data/distance-matrix.tsv", header = T, row.names = 1))
jac16Sm <- melt(jac16S)
colnames(jac16Sm) <- c("sample1", "sample2", "Jac16S")
jac16Sm$simDist <- 0 + jac16Sm$Jac16S
jac16Sm$Jac16S <- NULL
jac16Sm$distmet <- "Jac16S"

totaldistmet <- rbind(MSpcoa, jac16Sm)

totaldist3 <- merge(geodistm, totaldistmet, by = c("sample1", "sample2"))
# remove samples with the same 
totaldist3 <- totaldist3[totaldist3$geodist != 0,]

jacMSstat <- totaldist3[totaldist3$distmet == "JacMS",]
jac16stat <- totaldist3[totaldist3$distmet == "Jac16S",]

summary(lm(jacMSstat$geodist ~ jacMSstat$simDist))$adj.r.squared 
summary(lm(jac16stat$geodist ~ jac16stat$simDist))$adj.r.squared 


pdf("commdistXgeodist.pdf", height = 6, width = 4)
ggplot(totaldist3, aes(x = geodist, y = simDist, color = distmet, factor = distmet)) +
  stat_smooth(method = "glm", se = F, color = "#B1B1B1", size = 2, alpha = 1, linetype = "solid") +
  geom_jitter(alpha = 0.1) +
  theme_bw() +
  scale_x_log10()
dev.off()

################################################
################################################
### other ways to perform PCA
################################################
################################################


library(ape)

PCOA <- pcoa(dist.mat)
# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
PCOA$values$Relative_eig[1:3]
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Plot your results
biplot(PCOA)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
# However, we could work around this problem like this:
bipcoa <- biplot.pcoa(PCOA, paretotflex)


PCOAaxes <- PCOA$vectors[,c(1,2)]
PCOAmeta = metadata[, c(1,2,3)]

PCOAplot <- merge(PCOAaxes, PCOAmeta, by = "row.names")

# should look the same as the above one, just without the nice graphing stuff
ggplot(data = PCOAplot, aes(Axis.1, Axis.2, colour = ATTRIBUTE_Site)) +
  geom_point(size = 6, stroke = 2) +
  scale_colour_manual(values = staxcolors) +
  theme_bw()

