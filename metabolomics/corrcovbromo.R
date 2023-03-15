setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/metabolomics")

brominated <- read.table("Feat772_452_dibrominated.txt", header = T, sep = '\t', row.names = 1)
brominated[brominated <= 1E3] <- NA
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
brominated[is.na(brominated)] <- my.min(brominated) / 10
brominated$sampleID <- rownames(brominated)


bgccov <- read.table("BGCcoverages.txt", header = T, sep = '\t', row.names = 1)

bgcvar <- bgccov[ , grepl( "var" , names( bgccov ) ) ]
bgcvar$contigID <- rownames(bgcvar)
bgccov2 <- bgccov[ , !grepl( "var" , names( bgccov ) ) ]
bgccov2$contigID <- rownames(bgccov2)

library(reshape2)
mbgcvar <- melt(bgcvar)
mbgcvar$sampleID <- lapply(strsplit(as.character(mbgcvar$variable), "\\."), "[", 2)
mbgcvar <- mbgcvar[ , c("sampleID", "contigID", "value")]
colnames(mbgcvar) <- c("sampleID", "contigID", "variance")

mbgccov <- melt(bgccov2)
mbgccov$sampleID <- lapply(strsplit(as.character(mbgccov$variable), "\\."), "[", 2)
mbgccov <- mbgccov[ , c("sampleID", "contigID", "value")]

totalcov <- merge(mbgccov, mbgcvar, by = c("contigID", "sampleID"))
totalcov$sampleID <- vapply(totalcov$sampleID, paste, collapse = ", ", character(1L))

totalbromo <- merge(totalcov, brominated, by = "sampleID", all = T)
totalbromo$logINT <- log(totalbromo$peak_intensity)

### get the p-values extracted from the linear models
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

## summarize the topn this by taking slopes (expect positive), R2, and p-value using above function
library(dplyr)
slopes <- totalbromo %>% 
  group_by(contigID) %>% 
  do({
    mod = lm(value ~ logINT, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2],
               R2 = summary(mod)$r.squared,
               pvalue = lmp(mod))
  })

library(ggplot2)

ggplot(totalbromo, aes(x = logINT, y = value, color = contigID)) +
  geom_point() +
  theme_bw() +
  theme(legend.position="none") +
  geom_smooth(method = "lm", se = F)

### check for the myxo MAG BGC
myxo <- totalbromo[totalbromo$contigID == "MO_082.c_4",]

fit <- lm(data=myxo, logINT ~ value)

summary(aov(fit))$r.squared

ggplot(myxo, aes(x = logINT, y = value, color = contigID)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_errorbar(aes(ymin=value-variance, ymax=value+variance), width=.2, position=position_dodge(0.05)) +
  theme_bw()
  
tophit <- totalbromo[totalbromo$contigID == "MO_077.c_3101",]

ggplot(tophit, aes(x = logINT, y = value, color = contigID)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

## find if any top hits are NRPS related
# read in GCF calling
gcfIDs <- read.table("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/bigscape/bigscape_out/network_files/2022-10-20_12-08-51_hybrids_glocal/mix/mix_clustering_c0.60.tsv", header = T, sep = '\t', comment.char = "")
colnames(gcfIDs) <- c("Id", "GCFnum")
gcfIDs$contigID <- sub("^([^.]*.[^.]*).*", "\\1", gcfIDs$Id)

finalatt <- merge(slopes, gcfIDs, by = "contigID")
r2xGCF <- aggregate(finalatt$R2, by = list(finalatt$GCFnum), mean)

write.table(finalatt, file = "topcorrelatedBGCs.csv", sep = ',', quote = F, row.names = F)

################################################
################################################
#### try it with MAGs now from avg contig cov
################################################
################################################

magcov <- read.table("MAGcoverages.txt", header = T, sep = "\t", row.names = 1)

magvar <- magcov[ , grepl( "var" , names( magcov ) ) ]
magvar$magID <- rownames(magvar)
magcov2 <- magcov[ , !grepl( "var" , names( magcov ) ) ]
magcov2$magID <- rownames(magcov2)

library(reshape2)
mmagvar <- melt(magvar)
mmagvar$sampleID <- lapply(strsplit(as.character(mmagvar$variable), "\\."), "[", 2)
mmagvar <- mmagvar[ , c("sampleID", "magID", "value")]
colnames(mmagvar) <- c("sampleID", "magID", "variance")

mmagcov <- melt(magcov2)
mmagcov$sampleID <- lapply(strsplit(as.character(mmagcov$variable), "\\."), "[", 2)
mmagcov <- mmagcov[ , c("sampleID", "magID", "value")]

totalcov <- merge(mmagcov, mmagvar, by = c("magID", "sampleID"))
totalcov$sampleID <- vapply(totalcov$sampleID, paste, collapse = ", ", character(1L))

totalbromo <- merge(totalcov, brominated, by = "sampleID", all = T)
totalbromo$logINT <- log(totalbromo$peak_intensity)

slopes2 <- totalbromo %>% 
  group_by(magID) %>% 
  do({
    mod = lm(value ~ logINT, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2],
               R2 = summary(mod)$r.squared,
               pvalue = lmp(mod))
  })

myxo <- totalbromo[totalbromo$magID == "MO_102.bin.8",]

fit <- lm(data=myxo, logINT ~ value)

summary(aov(lm(data=myxo, logINT ~ value)))

ggplot(myxo, aes(x = logINT, y = value, color = magID)) +
  geom_point() +
  geom_smooth(method = "lm") +
  #geom_errorbar(aes(ymin=value-variance, ymax=value+variance), width=.2, position=position_dodge(0.05)) +
  theme_bw()

################################################
################################################
#### try it with MAGs avg cov from BBMap
################################################
################################################

magcov <- read.table("../MAGs/assemblies/binned_MAGs/avgcov.txt", header = T, sep = "\t")
magcov2 <- na.omit(magcov)

#### need to calculate RPKM with genome size instead of gene length
#### RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )

## import genome size
magstats <- read.table("../MAGs/assemblies/binned_MAGs/total.magstats.txt", header = T, sep = "\t")
magstats <- magstats[, c(1,3)]
magstats$binID <- magstats$Sequence_ID

totalmag <- merge(magcov2, magstats, by = "binID")
totalmag$RPKM <- ((totalmag$mappedreads) / ( (totalmag$Genome_length / 1000) * (totalmag$totalreads / 1000000) ) )

ggplot(totalmag) +
  geom_point(aes(x = RPKM, y = magcov))

totalbromo <- merge(totalmag, brominated, by = "sampleID", all = T)
totalbromo$logINT <- log(totalbromo$peak_intensity)

slopes <- totalbromo %>% 
  group_by(binID) %>% 
  do({
    mod = lm(RPKM ~ logINT, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2],
               R2 = summary(mod)$r.squared,
               pvalue = lmp(mod))
  })

myxo <- totalbromo[totalbromo$binID == "MO_082.bin.26",]

fit <- lm(data=myxo, logINT ~ RPKM)
summary(aov(lm(data=myxo, logINT ~ RPKM)))

### get colors for graph
mappingfile <- read.table("../mgenomeID.txt", header = T, sep = '\t', comment.char = "")
myxo <- merge(myxo, mappingfile, by = "sampleID")

staxcolors <- setNames(as.character(myxo$sitecolor), myxo$site)

p1 <- ggplot(myxo, aes(x = RPKM, y = logINT)) +
  geom_point(aes(color = site), size = 4) +
  geom_smooth(method = "lm", color = "#D3D3D3", fill = "lightgray") +
  ylim(c(4.9, 10)) +
  theme_bw() +
  scale_color_manual(values = staxcolors) +
  theme(legend.position="none")


bromoXsite <- merge(brominated, mappingfile, by = "sampleID")
bromoXsite$logINT <- log(bromoXsite$peak_intensity)

model <- aov(peak_intensity ~ site, data=bromoXsite)
summary(model)
TukeyHSD(model, conf.level=.95)

p2 <- ggplot(bromoXsite, aes(x = site, y = logINT, fill = site)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom="point", shape=23, size=4) +
  ylim(c(4.9, 10)) +
  scale_x_discrete(limits=c("site6", "site1", "site2", "site3", "site4", "site5", "site7", "site8")) +
  theme_bw() + 
  scale_fill_manual(values = staxcolors) +
  theme(legend.position="none")


library(gridExtra)
pdf("dibromoXsite.pdf", height = 4, width = 6)
grid.arrange(p2, p1, nrow = 1)
dev.off()

