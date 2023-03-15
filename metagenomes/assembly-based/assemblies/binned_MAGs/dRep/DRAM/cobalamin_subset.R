rm(list=ls())
setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/binned_MAGs/dRep/binMAGanalysis/DRAM")

miscmeta <- read.table("misc.metabolism.tsv", header = T, sep = '\t', quote = "")
miscmeta2 <- miscmeta[, -c(1,2,3,5)]

rownames(miscmeta2) <- miscmeta2$header
library(reshape2)

miscmeta2m <- melt(miscmeta2)
aggmisc <- aggregate(value ~ header + variable, data = miscmeta2m, sum, na.rm = TRUE)
miscmeta3 <- dcast(aggmisc, header ~ variable, value.var = "value")
rownames(miscmeta3) <- miscmeta3$header
miscmeta3$header <- NULL
miscmeta3 <- as.matrix(miscmeta3)

row_names_df_to_remove <- c("Ribosome, archaea", "Ribosome, bacteria")
miscmeta4 <- miscmeta3[!(row.names(miscmeta3) %in% row_names_df_to_remove),]

heatmap.2(x=t(miscmeta3), scale="column", col="bluered")

library(gplots)

heatmap2()