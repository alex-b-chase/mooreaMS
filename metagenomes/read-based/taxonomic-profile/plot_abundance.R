setwd('/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/gene-centric-approach/taxonomic-profile/pplacer_output/pplacer_analysis/')

library(ggplot2)
library(reshape2) # for melt

abund <- read.table('phylaXsite_abundance_matrix.txt', header = T, sep = '\t')
abundm <- melt(abund, "sampleID")

### remove site 6 since it looks different in 16S data
site6 <- c("MO_167", "MO_170", "MO_180", "MO_184", "MO_188", "MO_192")
abundm2 <- abundm[ ! abundm$sampleID %in% site6, ]

mappingfile <- read.table('/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/mgenomeID.txt', header = T, sep = '\t', comment.char = "")

abundmtotal <- merge(abundm2, mappingfile, by = "sampleID")

abundmtotal$site = factor(abundmtotal$site, levels=c('site1', 'site2', 'site3', 'site4', 
                                                     'site5', 'site6', 'site7', 'site8'))

# read in phyla colors
colors <- read.table('phyla_colors.txt', header = T, sep ='\t', stringsAsFactors = T, comment.char = '')
# create new color pallette for ggplot
taxcolors <- setNames(as.character(colors$color), colors$variable)

meanabundm <- aggregate(value ~ variable + site + sampleID, abundmtotal, mean)

meanabundm$variable = factor(meanabundm$variable, levels = c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Firmicutes", 
                                                     "Planctomycetes", "Cyanobacteria", "Acidobacteria","Chloroflexi", "Other", "NA."))


# barplot of transplant and inoculum

pdf('phylaXMGID_nosite6.pdf', width = 14, height = 10)

ggplot(meanabundm, aes(y = value, x = sampleID, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') + facet_grid( ~ site, scales="free_x") +
  scale_fill_manual(values = taxcolors, name = "") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y = "% Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())

dev.off()


# PCO plot of abundances 
# read in the edgenum X site

# first need to do rarefaction
# Rarefying - normalizes read depth across all samples. 
# Allows for an equal comparison across different sample types, at the risk of excluding rarer taxa
# A "good" rarefaction depth should minimize sample loss while maximizing OTU richness.

library(vegan)
library(ggplot2)
library(EcolUtils)


# get quartile ranges for rarefaction
finalOTU <- read.table("../siteXedge_num_raw.txt", header = T)

finalOTU2 <- subset(finalOTU,!rownames(finalOTU) %in% site6)

transOTU <- rowSums(finalOTU2)
Q10 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.10)
Q15 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.15)
Q25 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.25)

barplot(sort(transOTU), ylim = c(0, max(transOTU)), 
        xlim = c(0, NROW(transOTU)), col = "Blue", ylab = "Read Depth", xlab = "Sample") 
abline(h = c(Q10, Q15, Q25), col = c("red", "pink", "yellow"))
plot.new()

rarecurve(finalOTU2, step = 100, cex = 0.1)
abline(v = c(Q10, Q15), col = c("red", "pink"))

rared_OTU <- as.data.frame((rrarefy.perm(finalOTU2, sample = Q10, n = 100, round.out = T)))
#This only keeps the samples that meet the rarefaction cutoff.
rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= Q10 - (Q10 * 0.1), colSums(rared_OTU) >= 1])

NMDS1 <- metaMDS(rared_OTU, distance = "bray", k = 2, trymax = 500, autotransform = F)

#This will make the first two columns as the x and y coordinates, so that they may be plotted easily.
coordinates <- data.frame(NMDS1$points[,1:2])

#Lets see how this looks like by actually plotting the data. 
#The plot can be viewed in the plots tab in the bottom right quadrant.
plot(x = coordinates$MDS1, y = coordinates$MDS2)

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
metadata <- mappingfile
rownames(metadata) <- metadata$sampleID
nmds_plus_metadata <- merge(coordinates, metadata, by = 0)

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

samplecolors <- setNames(as.character(nmds_plus_metadata$sitecolor), nmds_plus_metadata$site)

pdf("MDSedges_rarefied_averaged_nosite6.pdf", height = 8, width = 14)
ggplot(nmds_plus_metadata, aes(MDS1, MDS2, shape = reef)) +
  geom_point(aes(color = site), size = 12) +
  geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()

### now try with relative abundances
abund <- read.table('edgenumXsite_abundance_matrix.txt', header = T, sep = '\t')

library(vegan)
row.names(abund) <- abund$sampleID

### remove site 6 since it looks different in 16S data
site6 <- c("MO_167", "MO_170", "MO_180", "MO_184", "MO_188", "MO_192")
abund2 <- abund[ ! abund$sampleID %in% site6, ]

abund2$sampleID <- NULL

braycurt <- vegdist(abund2, method = "bray")
hClustering <- hclust(braycurt, method = 'complete')
plot(hClustering, hang = -1)

sol <- metaMDS(abund2, distance = "bray", k = 2, trymax = 500, autotransform = F)

row.names(mappingfile) <- mappingfile$sampleID

NMDS1 = data.frame(x = sol$point[, 1], y = sol$point[, 2])
NMDS2 = mappingfile

NMDS <- merge(NMDS1, NMDS2, by = "row.names")

row.names(NMDS) <- NMDS$Row.names
NMDS$Row.names <- NULL

pdf("MDSedges_nosite6_abund.pdf", height = 8, width = 14)

ggplot(NMDS, aes(x, y)) +
  geom_point(aes(color = site), size = 10, shape = 15) +
  #geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  #scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()

### permanova to test for site effects
metadata <- mappingfile
rownames(metadata) <- metadata$sampleID

tokeep <- rownames(NMDS1)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]

perm.log <- adonis2(braycurt ~ reef * site, data=metadata2, permutations = 999, method = "bray")
perm.log

