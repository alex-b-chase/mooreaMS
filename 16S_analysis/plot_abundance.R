setwd('/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/16S_analysis/plot_nosite6/')

library(ggplot2)
library(reshape2) # for melt

abund <- read.table('phylaXsite_abundance_matrix.txt', header = T, sep = '\t')

abundm <- melt(abund, "sampleID")

### remove site 6 since it looks different in 16S data
site6 <- c("MO18_161", "MO18_162", "MO18_163", "MO18_164", "MO18_165", "MO18_166", "MO18_167", 
           "MO18_168", "MO18_169", "MO18_170", "MO18_171", "MO18_172", "MO18_173", "MO18_174", 
           "MO18_175", "MO18_176", "MO18_177", "MO18_178", "MO18_179", "MO18_180", "MO18_181", 
           "MO18_182", "MO18_183", "MO18_184", "MO18_185", "MO18_186", "MO18_187", "MO18_188", 
           "MO18_189", "MO18_190", "MO18_191", "MO18_192")
abundm2 <- abundm[ ! abundm$sampleID %in% site6, ]

mappingfile <- read.table('16S_metadata.txt', header = T, sep = '\t', comment.char = "")

abundmtotal <- merge(abundm2, mappingfile, by = "sampleID")

abundmtotal$site = factor(abundmtotal$site, levels=c('site1', 'site2', 'site3', 'site4', 
                                                     'site5', 'site6', 'site7', 'site8'))

# read in phyla colors
colors <- read.table('phyla_colors.txt', header = T, sep ='\t', stringsAsFactors = T, comment.char = '')
# create new color pallette for ggplot
taxcolors <- setNames(as.character(colors$color), colors$variable)

meanabundm <- aggregate(value ~ variable + site, abundmtotal, mean)

meanabundm$variable = factor(meanabundm$variable, 
                             levels = c("Proteobacteria", "Desulfobacterota","Actinobacteria", "Bacteroidetes", "Firmicutes", 
                                        "Planctomycetes", "Cyanobacteria", "Acidobacteria","Chloroflexi", "Other", "Archaea", "NA."))


# barplot of transplant and inoculum

pdf('phylaXMGID.pdf', width = 6, height = 10)

ggplot(meanabundm, aes(y = value, x = site, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') +
  #+ facet_grid( ~ site, scales="free_x") +
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

### ran many of these diversity metrics on QIIME2
### read in those distance matrices and plot
library(vegan)
library(ggplot2)
library(EcolUtils)

unifracQ <- as.matrix(read.table("../qiime2_nosite6/coremetrics/weightedunifrac_dist/data/distance-matrix.tsv", header = T, row.names = 1))

NMDSQ <- metaMDS(unifracQ, distance = "bray", k = 3, trymax = 500, autotransform = T)

#This will make the first two columns as the x and y coordinates, so that they may be plotted easily.
coordinates <- data.frame(NMDSQ$points[,1:2])

#Lets see how this looks like by actually plotting the data. 
#The plot can be viewed in the plots tab in the bottom right quadrant.
# plot(x = coordinates$MDS1, y = coordinates$MDS2)

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
metadata <- mappingfile
rownames(metadata) <- metadata$sampleID

tokeep <- rownames(coordinates)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]
nmds_plus_metadata <- merge(coordinates, metadata, by = 0)

samplecolors <- setNames(as.character(nmds_plus_metadata$sitecolor), nmds_plus_metadata$site)


pdf("MDSedges_rarefied_unifrac.pdf", height = 10, width = 7)
#ggplot(nmds_plus_metadata, aes(MDS1, MDS2, shape = ecology)) +
ggplot(nmds_plus_metadata, aes(MDS1, MDS2)) +
  geom_point(aes(color = site), size = 6, shape = 15) +
  #geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  #scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()

### permanova to test for site effects
dist.mat <- vegdist(unifracQ, method = "bray")
perm.log <- adonis2(dist.mat ~ reef * site, data=metadata2, permutations = 999, method = "bray")
perm.log


# PCO plot of abundances by inoculum and transplant
# read in the edgenum X site

# first need to do rarefaction
# Rarefying - normalizes read depth across all samples. 
# Allows for an equal comparison across different sample types, at the risk of excluding rarer taxa
# A "good" rarefaction depth should minimize sample loss while maximizing OTU richness.

library(vegan)
library(ggplot2)
library(EcolUtils)


# get quartile ranges for rarefaction
finalOTU2 <- read.table("../qiime2_nosite6/feature-table-nosite6.tsv", sep = '\t', header = T, row.names = 1)
finalOTU3 <- as.data.frame(t(finalOTU2))

# site6 <- c("MO18_161", "MO18_162", "MO18_163", "MO18_164", "MO18_165", "MO18_166", "MO18_167", 
#            "MO18_168", "MO18_169", "MO18_170", "MO18_171", "MO18_172", "MO18_173", "MO18_174", 
#            "MO18_175", "MO18_176", "MO18_177", "MO18_178", "MO18_179", "MO18_180", "MO18_181", 
#            "MO18_182", "MO18_183", "MO18_184", "MO18_185", "MO18_186", "MO18_187", "MO18_188", 
#            "MO18_189", "MO18_190", "MO18_191", "MO18_192")
# finalOTU4 <- subset(finalOTU3,!rownames(finalOTU3) %in% site6)
finalOTU <- finalOTU3

# ### remove singletons from the ASV table (dunno why QIIME didn't work on this...)
# ### must be it removed ASV with a frequency of 1, not only in 1 sample....
# ### can remove more rarer taxa - doesn't affect beta diversity patterns anyways

transOTU <- rowSums(finalOTU)
Q01 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.001)
Q02 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.0025)
Q03 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.003)
Q05 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.005)
Q1 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.01)
Q5 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.05)
Q10 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.10)
Q15 <- quantile(transOTU[order(transOTU, decreasing = TRUE)], 0.15)

# barplot(sort(transOTU), ylim = c(0, max(transOTU)),
#         xlim = c(0, NROW(transOTU)), col = "Blue", ylab = "Read Depth", xlab = "Sample")
# abline(h = c(Q05, Q1, Q5, Q10, Q15), col = c("seagreen", "red", "pink", "yellow", "blue"))
# plot.new()

# pdf("rarefaction.pdf", height = 8, width = 8)
# rarecurve(finalOTU, step = 100, cex = 0.1)
# abline(v = c(Q01, Q02, Q03, Q05, Q1, Q5, Q10, Q15), col = c("green", "dodgerblue", "pink", "seagreen", "red", "pink", "yellow", "blue"))
# dev.off()

Q05
Q1
Q5

QINPUT <- Q1
## Q1 == 27757.75

rared_OTU <- as.data.frame((rrarefy.perm(finalOTU, sample = QINPUT, n = 100, round.out = T)))
## This only keeps the samples that meet the rarefaction cutoff.
## Also can remove rare taxa with <5 sequences, less likely to impact beta diversity, but cuts down computation
rared_OTU <- as.data.frame(rared_OTU[rowSums(rared_OTU) >= QINPUT - (QINPUT * 0.1), colSums(rared_OTU) >= 0])

NMDS1 <- metaMDS(rared_OTU, distance = "bray", k = 3, trymax = 500, autotransform = T)

#This will make the first two columns as the x and y coordinates, so that they may be plotted easily.
coordinates <- data.frame(NMDS1$points[,1:2])

#Lets see how this looks like by actually plotting the data. 
#The plot can be viewed in the plots tab in the bottom right quadrant.
# plot(x = coordinates$MDS1, y = coordinates$MDS2)

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
metadata <- mappingfile
rownames(metadata) <- metadata$sampleID

tokeep <- rownames(coordinates)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]
nmds_plus_metadata <- merge(coordinates, metadata, by = 0)

samplecolors <- setNames(as.character(nmds_plus_metadata$sitecolor), nmds_plus_metadata$site)


pdf("MDSedges_rarefied_bray.pdf", height = 10, width = 7)
#ggplot(nmds_plus_metadata, aes(MDS1, MDS2, shape = ecology)) +
ggplot(nmds_plus_metadata, aes(MDS1, MDS2)) +
  geom_point(aes(color = site), size = 6, shape = 15) +
  #geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  #scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()

### permanova to test for site effects
dist.mat <- vegdist(rared_OTU, method = "bray")
perm.log <- adonis2(dist.mat ~ reef * site, data=metadata2, permutations = 999, method = "bray")
perm.log

write.csv(as.matrix(dist.mat),'braycurt_dist.csv',quote=F);


### cannot compute these metrics with the current size of the dataframe
### the only way to get them, was to eliminate a lot of rare taxa, not ideal
### with that method, we get a betaMNTD << -2 == homogeneous selection

### we can sub-sample the dataframe to look at 8 samples per plot = 64 total samples
### do the numbers hold up?????

subsampleIDs <- c(
  "MO18_001", "MO18_003", "MO18_005", "MO18_007", "MO18_009", "MO18_011", "MO18_013", "MO18_015", "MO18_020", "MO18_024", "MO18_028", "MO18_032",
  "MO18_033", "MO18_035", "MO18_037", "MO18_039", "MO18_041", "MO18_043", "MO18_045", "MO18_047", "MO18_052", "MO18_056", "MO18_060", "MO18_064",
  "MO18_065", "MO18_067", "MO18_069", "MO18_071", "MO18_073", "MO18_075", "MO18_077", "MO18_079", "MO18_084", "MO18_088", "MO18_092", "MO18_096",
  "MO18_097", "MO18_099", "MO18_101", "MO18_103", "MO18_105", "MO18_107", "MO18_109", "MO18_111", "MO18_116", "MO18_120", "MO18_124", "MO18_128",
  "MO18_129", "MO18_131", "MO18_133", "MO18_135", "MO18_137", "MO18_139", "MO18_141", "MO18_143", "MO18_148", "MO18_152", "MO18_156", "MO18_160",
  "MO18_193", "MO18_195", "MO18_197", "MO18_199", "MO18_201", "MO18_203", "MO18_205", "MO18_207", "MO18_212", "MO18_216", "MO18_220", "MO18_224",
  "MO18_225", "MO18_227", "MO18_229", "MO18_231", "MO18_233", "MO18_235", "MO18_237", "MO18_239", "MO18_244", "MO18_248", "MO18_252", "MO18_256")
subOTU2 <- subset(finalOTU,rownames(finalOTU) %in% subsampleIDs)
## remove OTUs from other samples
subOTU <- subOTU2[, !colSums(subOTU2 == 0) >= 84, drop = FALSE]
## remove rare OTUs with frequency of <5
subOTU[subOTU == 1] <- NA
subOTU[subOTU == 2] <- NA
subOTU[subOTU == 3] <- NA
subOTU[subOTU == 4] <- NA
subOTU[subOTU == 5] <- NA
subOTU2 <- subOTU[ , colSums(is.na(subOTU)) == 0]

subRARE <- as.data.frame((rrarefy.perm(subOTU2, sample = QINPUT, n = 100, round.out = T)))
## This only keeps the samples that meet the rarefaction cutoff.
## Also can remove rare taxa with <5 sequences, less likely to impact beta diversity, but cuts down computation
subRARE <- as.data.frame(subRARE[rowSums(subRARE) >= QINPUT - (QINPUT * 0.1), colSums(subRARE) >= 0])

subNMDS <- metaMDS(subRARE, distance = "bray", k = 3, trymax = 500, autotransform = T)

#This will make the first two columns as the x and y coordinates, so that they may be plotted easily.
coordinates <- data.frame(subNMDS$points[,1:2])

#Lets see how this looks like by actually plotting the data. 
#The plot can be viewed in the plots tab in the bottom right quadrant.
# plot(x = coordinates$MDS1, y = coordinates$MDS2)

#To make a more sophisticated plot, we will merge the stress scores with the metadata.
metadata <- mappingfile
rownames(metadata) <- metadata$sampleID

tokeep <- rownames(coordinates)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]
nmds_plus_metadata <- merge(coordinates, metadata, by = 0)

samplecolors <- setNames(as.character(nmds_plus_metadata$sitecolor), nmds_plus_metadata$site)


pdf("MDSsub_rarefied_bray.pdf", height = 12, width = 14)
#ggplot(nmds_plus_metadata, aes(MDS1, MDS2, shape = ecology)) +
ggplot(nmds_plus_metadata, aes(MDS1, MDS2)) +
  geom_point(aes(color = site), size = 12) +
  geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  #scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()


##### Determination of Assembly Processes #####
### By Stacey J. Doherty ###
### last modified 8 October 2020 ###

## This script was used to identify assembly processes structuring bacterial communities using 16S rRNA data.
## This script is intended to accompany the publication Doherty et al. 2020. "The transition from stochastic to deterministic bacterial community assembly during permafrost thaw succession"
## The original code was published by Stegen et al. 2013. "Quantifying community assembly processes and identifying features that impose them." ISME J.
## Please cite the original paper if using this script! I have merely added notes for clarification. Orginial code can be found on Stegen GitHub.

Sys.setenv('R_MAX_VSIZE'=128000000000)
##### This section will generate the bNTI pairwise matrix for the dataset #####

## Change to the directory on your computer that contains the ASV table and associated phylogeny
## I used the rarefied ASV table and phylogeny. This one phylogeny will include all ASVs present in your dataset.
## Note that the 'slash' needs to be a forward slash like this /
## Load the picante library (Kembel et al., 2010)
## if not already installed, use install.packages('picante')
library(picante)

## Read in ASV table. It should be formatted as a .csv file. Rows are ASVs and columns are samples.
## This table was generated in qiime2, exported to a .biom file, then converted to a .csv file (code not shown).
## Note, Stegen et al. (2013) performed this analysis with OTU data, but for the sake of this code OTU = ASV.
# otu = read.table("feature-table-all.tsv",header=T,row.names=1);
# dim(otu); # this gives the dimensions
# otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## Read in the phylogeny. It should be formatted as a .nwk file. 
## Phylogeny was constructed in qiime2, exported as a .nwk format
phylo = read.tree("../qiime2_output/rooted_phylo/data/tree.nwk");
# phylo; # a summary of the phylogeny
#plot.phylo(phylo,typ="fan"); # a quick plot, this may not work well depending on computer specs, but it isn't necessary to run.

## Make sure the names on the phylogeny are ordered the same as the names in ASV table
match.phylo.otu = match.phylo.data(phylo, t(subRARE));
# str(match.phylo.otu);

## Calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
# dim(beta.mntd.weighted);
# beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'betaMNTDsub_weighted.csv',quote=F);

## plot the betaMNTD as MDS plot
betaMNTD_MDS <- metaMDS(beta.mntd.weighted, k = 3, trymax = 500)
betaMNTDcoord <- data.frame(betaMNTD_MDS$points[,1:2])

tokeep <- rownames(betaMNTDcoord)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]
betaMNTD_plus_meta <- merge(betaMNTDcoord, metadata, by = 0)

samplecolors <- setNames(as.character(betaMNTD_plus_meta$sitecolor), betaMNTD_plus_meta$site)


pdf("MDSsub_rarefied_betaMNTD.pdf", height = 12, width = 14)
#ggplot(nmds_plus_metadata, aes(MDS1, MDS2, shape = ecology)) +
ggplot(betaMNTD_plus_meta, aes(MDS1, MDS2)) +
  geom_point(aes(color = site), size = 12) +
  geom_text(aes(label = sampleID)) +
  theme_bw() +
  scale_colour_manual(values = samplecolors) +
  #scale_shape_manual(values = c(15,18)) +
  scale_fill_manual(values = samplecolors)

dev.off()

### permanova to test for site effects
dist.mat <- vegdist(subRARE, method = "bray")
perm.log <- adonis2(dist.mat ~ site, data=metadata2, permutations = 999, method = "bray")
perm.log


identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# Calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

pdf("weighted_bNTI_Histogram.pdf")
hist(weighted.bNTI)
dev.off()

## The weighted_bNTI.csv contains a pairwise matrix of bNTI values for the dataset.
## I opened the file in excel to assess and document assembly processes. 
## -2 < bNTI < 2 indicates stochastic process which can be further refined using values from RCbray in next section.
## bNTI < -2 indicates homogeneous selection
## bNTI > 2 indicates heterogeneous selection



##### This section will generate the RCbray pairwise matrix for the dataset #####

## This script is quite computationally heavy, I recommend running it on a super computer. 

## Change to the directory on your computer that contains the ASV table and associated phylogeny
## I used the rarefied ASV table and phylogeny. This one phylogeny will include all ASVs present in your dataset.
## Note that the 'slash' needs to be a forward slash like this /

## Read in ASV table 
## This table was generated in qiime2, exported to a .biom file, then converted to a .csv file (code not shown)
asv = t(subRARE)
dim(asv); # this gives the dimensions
asv[1:5,1:5]; # this gives a look at the first 5 rows and columns

spXsite1 = asv

# Load vegan package 
library(vegan)

## Only deviation from Stegen code: SJD edited this function to use vegan for distance matrix creation.
raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  
  ##By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). 
  ##Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. 
  ##The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number 
  ##of shared species to the calculation- this is highly recommended.  
  ## The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  
  ## Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  
  ## If ties are split (as we recommend) the dissimilarity (default) and similarity 
  ## (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) 
  ## or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
  ## If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. 
  ## The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  
  ## set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  
  ## Note that the choice of how many plots (rows) to include has a real impact on the metric,
  ## as species and their occurrence frequencies across the set of plots is used to determine gamma and 
  ## the frequency with which each species is drawn from the null model	
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null bray curtis
        null_bray_curtis[i] = vegdist(null.spXsite,method='bray');
        
      }; # end reps loop
      
      ## empirically observed bray curtis
      obs.bray = vegdist(spXsite[c(null.one,null.two),],method='bray');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }	
  
  return(results)
  
}; ## end function

raup_crick_abundance(spXsite1, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE)

## The RCbray pairwise matrix will print to the terminal (console). This code could be adapted to write the results to a .csv file.
## I copied the results into a .csv file. The file was opened in excel to further refine the stochastic assembly processes. 
## -2 < bNTI < 2 and RCbray < -0.95 indicates homogenizing dispersal
## -2 < bNTI < 2 and RCbray > 0.95 indicates dispersal limitation and drift
## -2 < bNTI < 2 and -0.95 < RCbray < 0.95 indicates drift



