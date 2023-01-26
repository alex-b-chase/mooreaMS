setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/16S_analysis")

## compare the weighted unifrac method to comdistnt (picante package)
## unifrac is from QIIME2 output, comdistnt from "plot_abundance.R"

betaMNTD <- as.matrix(read.csv("plot_nosite6/betaMNTD_weighted_nosite6.csv", header = T, row.names = 1))
unifracQ <- as.matrix(read.table("qiime2_nosite6/coremetrics/weightedunifrac_dist/data/distance-matrix.tsv", header = T, row.names = 1))

library(reshape2)
betaMNTDm <- melt(betaMNTD)
colnames(betaMNTDm) <- c("sample1", "sample2", "betaMNTD")
unifracQm <- melt(unifracQ)
colnames(unifracQm) <- c("sample1", "sample2", "Wunifrac")

totaldist <- merge(unifracQm, betaMNTDm, by = c("sample1", "sample2"))

library(ggplot2)

summary(lm(totaldist$Wunifrac ~ totaldist$betaMNTD))$adj.r.squared 

pdf("weightedUnifracXbetaMNTD.pdf", height = 6, width = 10)
ggplot(totaldist, aes(x = Wunifrac, y = betaMNTD)) +
  stat_smooth(se=F,
              method = "lm", color = "#B1B1B1", size = 2, alpha = 1, linetype = "solid") +
  geom_jitter(alpha = 0.1) +
  theme_bw()
dev.off()


## similarly compare the two bray curtis calculations from QIIME2 and R
brayR <- as.matrix(read.csv("plot_nosite6/braycurt_dist.csv", header = T, row.names = 1))
brayQ <- as.matrix(read.table("qiime2_nosite6/coremetrics/braycurt_dist/data/distance-matrix.tsv", header = T, row.names = 1))

brayRm <- melt(brayR)
colnames(brayRm) <- c("sample1", "sample2", "brayR")
brayQm <- melt(brayQ)
colnames(brayQm) <- c("sample1", "sample2", "brayQ")

totaldist2 <- merge(brayRm, brayQm, by = c("sample1", "sample2"))

summary(lm(totaldist2$brayR ~ totaldist2$brayQ))$adj.r.squared 

pdf("brayRXbrayQ.pdf", height = 6, width = 10)
ggplot(totaldist2, aes(x = brayQ, y = brayR)) +
  stat_smooth(se=F,
              method = "lm", color = "#B1B1B1", size = 2, alpha = 1, linetype = "solid") +
  geom_jitter(alpha = 0.1) +
  theme_bw()
dev.off()


# Spatial Mantel correlogram analysis
## make geographical distance matrix
metadata <- read.table('plot_nosite6/16S_metadata.txt', header = T, sep = '\t', comment.char = "", row.names = 1)
tokeep <- rownames(brayR)
metadata2 <- metadata[which(rownames(metadata) %in% tokeep), ]

library(geosphere)
library(mpmcorrelogram)
Geodist <- as.dist(distm(cbind(metadata2$lon, metadata2$lat), fun=distGeo))

braydist <- as.dist(brayR)
bMNTDdist <- as.dist(betaMNTD)
unidist <- as.dist(unifracQ)

brayGeoMCA <- mpmcorrelogram(braydist, Geodist, method="spearman", permutations=999)
bMNTDGeoMCA <- mpmcorrelogram(bMNTDdist, Geodist, method="spearman", permutations=999)
uniGeoMCA <- mpmcorrelogram(unidist, Geodist, method="spearman", permutations=999)

# ## plot results
# pdf("GeoMCA.pdf", width=7, height=7)
# tipos <- brayGeoMCA$pval.Bonferroni < 0.05
# tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
# xval <- c()
# for(j in 1:(length(brayGeoMCA$breaks) - 1)) {
#   xval[j] <- (brayGeoMCA$breaks[j] + brayGeoMCA$breaks[j + 1]) / 2
# }
# par(mar=c(5, 4, 4, 4))
# plot(Geodist, braydist, pch=16, cex=1, type="p", xlim=c(0, tail(brayGeoMCA$breaks, n=1)), ylim=c(0,1), col=alpha("black", 0.3), ann=F, axes=F)
# axis(4)
# mtext("Bray-Curtis distance", side=4, line=2.8)
# par(new=T)
# plot(xval, brayGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(brayGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="BrayCurtisGeoMCA")
# abline(v=brayGeoMCA$breaks, lty=2)
# 
# tipos <- bMNTDGeoMCA$pval.Bonferroni < 0.05
# tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
# xval <- c()
# for(j in 1:(length(bMNTDGeoMCA$breaks) - 1)) {
#   xval[j] <- (bMNTDGeoMCA$breaks[j] + bMNTDGeoMCA$breaks[j + 1]) / 2
# }
# par(mar=c(5, 4, 4, 4))
# plot(Geodist, bMNTDdist, pch=16, cex=1, type="p", xlim=c(0, tail(bMNTDGeoMCA$breaks, n=1)), ylim=c(0,0.1), col=alpha("black", 0.3), ann=F, axes=F)
# axis(4)
# mtext("betaMNTD distance", side=4, line=2.8)
# par(new=T)
# plot(xval, bMNTDGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(bMNTDGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="betaMNTDGeoMCA")
# abline(v=bMNTDGeoMCA$breaks, lty=2)
# 
# tipos <- bMNTDGeoMCA$pval.Bonferroni < 0.05
# tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,15,22))
# xval <- c()
# for(j in 1:(length(bMNTDGeoMCA$breaks) - 1)) {
#   xval[j] <- (bMNTDGeoMCA$breaks[j] + bMNTDGeoMCA$breaks[j + 1]) / 2
# }
# par(mar=c(5, 4, 4, 4))
# plot(Geodist, bMNTDdist, pch=16, cex=1, type="p", xlim=c(0, tail(bMNTDGeoMCA$breaks, n=1)), ylim=c(0,0.1), col=alpha("black", 0.3), ann=F, axes=F)
# axis(4)
# mtext("betaMNTD distance", side=4, line=2.8)
# par(new=T)
# plot(xval, bMNTDGeoMCA$rM, pch=tipos, cex=1, type="b", xlim=c(0, tail(bMNTDGeoMCA$breaks, n=1)), ylim=c(-1,1), xlab="geographic distance", ylab="Mantel correlation", main="betaMNTDGeoMCA")
# abline(v=bMNTDGeoMCA$breaks, lty=2)
# 
# dev.off()
# 
# 
# p <- as.matrix(Geodist)
# write.csv(p, "geodist.csv", quote = F)
### calculate this and now import new distances adjusted for plot spacing

geodist <- as.matrix(read.table("geodist.txt", header = T, row.names = 1))
library(reshape2)

geodistm <- melt(geodist)
colnames(geodistm) <- c("sample1", "sample2", "geodist")

brayRm <- melt(brayR)
colnames(brayRm) <- c("sample1", "sample2", "bray")
brayRm$simDist <- 1 - brayRm$bray
brayRm$bray <- NULL
brayRm$distmet <- "Bray"

unifracQm <- melt(unifracQ)
colnames(unifracQm) <- c("sample1", "sample2", "unifrac")
unifracQm$simDist <- 1 - unifracQm$unifrac
unifracQm$unifrac <- NULL
unifracQm$distmet <- "Unifrac"

totaldistmet <- rbind(unifracQm, brayRm)

totaldist3 <- merge(geodistm, totaldistmet, by = c("sample1", "sample2"))

# remove samples with the same 
totaldist3 <- totaldist3[totaldist3$geodist != 0,]

summary(lm(totaldist3$geodist ~ totaldist3$simDist * totaldist3$distmet))

braystat <- totaldist3[totaldist3$distmet == "Bray",]
uniFstat <- totaldist3[totaldist3$distmet == "Unifrac",]

summary(lm(braystat$geodist ~ braystat$simDist))$adj.r.squared 
summary(lm(uniFstat$geodist ~ uniFstat$simDist))$adj.r.squared 

pdf("commdistXgeodist.pdf", height = 6, width = 4)
ggplot(totaldist3, aes(x = geodist, y = simDist, color = distmet, factor = distmet)) +
  stat_smooth(method = "glm", se = F, color = "#B1B1B1", size = 2, alpha = 1, linetype = "solid") +
  geom_jitter(alpha = 0.1) +
  theme_bw() +
  scale_x_log10()
dev.off()
