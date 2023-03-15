setwd('/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/bigscape')

mixednet <- read.table('bigscape_out/network_files/2022-10-20_12-08-51_hybrids_glocal/mix/mix_c0.60.network', sep = '\t', header = T)

hist(mixednet$Raw.distance, breaks = 600)
hist(mixednet$DSS.index, breaks = 600)
hist(mixednet$Squared.similarity, breaks = 600)



finalletters <- mixednet[, c("Clustername.1", "Clustername.2", "Squared.similarity") ]
finalletters$type <- "Directed"
finalletters$BGC1 <-finalletters$Clustername.1
finalletters$BGC2 <- finalletters$Clustername.2

library(splitstackshape)
letterstemp <- cSplit(finalletters, "Clustername.1", ".")
letterstemp2 <- cSplit(letterstemp, "Clustername.2", ".")

finalletters2 <- letterstemp2[, c(3,4,2,1)]
# 
# # want to remove hits that only have MIBIG clusters in them
finalletters2$mibig1 <- NA
finalletters2$mibig2 <- NA
finalletters2$mibig1 <- ifelse(grepl("^MB", finalletters2$BGC1, ignore.case = T), "yes", "no")
finalletters2$mibig2 <- ifelse(grepl("^MB", finalletters2$BGC2, ignore.case = T), "yes", "no")
 
finalletters2$mibigcheck <- NA
finalletters2$mibigcheck <- paste(finalletters2$mibig1, finalletters2$mibig2, sep ="")
finalletters3 <- finalletters2[finalletters2$mibigcheck != "yesyes",]
 
finalletters3 <- finalletters3[, c(1,2,3,4)]
colnames(finalletters3) <- c("Source", "Target", "Type", "Weight")

write.table(finalletters3, file = "bgc_cluster_network.csv", row.names = F, sep = ',')
write.table(finalletters3, file = "bgc_cluster_network.tsv", row.names = F, sep = '\t', quote = F)

########################################################################
####  GET ATTRIBUTES FOR ALL NODES
########################################################################

### run parseBGCs.sh first to get input files
netatt <- read.table('network.nodes.uniq.txt', sep = '\t', header = F)
colnames(netatt) <- c("cluster1")

library(splitstackshape)
netatt$genome1 <- netatt$cluster1
netatt <- cSplit(netatt, "cluster1", ".")

colnames(netatt) <- c("BGC1", "genome1", "contig1", "region1")

colors <- read.table('../../../mgenomeID.txt', header =T, sep ='\t', stringsAsFactors = T, comment.char = '')
colors$genome1 <- colors$sampleID
colors <- colors[, c(7,2,6)]
library(dplyr)
colors <- distinct(colors)

colors2 <- merge(netatt, colors, by = "genome1", all = T)
colors2$region1 <- gsub( "region", "", as.character(colors2$region1))
colors2$contig1 <- gsub( "c_", "c", as.character(colors2$contig1))
colors2$cluster1 <- paste(colors2$contig1, colors2$region1, sep = "_")
colors2$cluster1 <- gsub( "region_", "", as.character(colors2$cluster1))

colors3 <- colors2[,c("BGC1", "site", "genome1" ,"cluster1", "sitecolor")]

colors4 <- colors3 %>% distinct

bgcatt <- read.table('mibig_BGC_innetwork.reduce10.txt', header = F, sep = '\t')
colnames(bgcatt) <- c("cluster1", "product")

colors6 <- merge(colors4, bgcatt, by = "cluster1", all = T)

# color known BGC by themselves to distinguish
colors6$clade2 <- ifelse(grepl("MB", colors6$BGC1, ignore.case = T), "knownBGC", 
                         paste(colors6$site))
colors6$size <- ifelse(grepl("BGC", colors6$BGC1, ignore.case = T), "5", 
                       "1")
colors6$cluster2 <- ifelse(grepl("MB", colors6$BGC1, ignore.case = T), paste(colors6$product), 
                           paste(colors6$genome1, colors6$cluster1, sep = '_'))

colors6 <- colors6[!is.na(colors6$BGC1),]

tmpatt <- colors6[, c(2,9,7,8)]
colnames(tmpatt) <- c("Id","Label","clade","size")

contigedge <- read.table("BGCinformation.txt", header = T, sep = '\t')
colnames(contigedge) <- c("Id", "BGCtype", "complete", "MAGhit", "MAG")
finalatt <- merge(tmpatt, contigedge, by = "Id", all = T)
finalatt <- as.matrix(finalatt)
finalatt[is.na(finalatt)] <- "NA."
finalatt <- as.data.frame(finalatt)

# get GCF info
gcfIDs <- read.table("bigscape_out/network_files/2022-10-20_12-08-51_hybrids_glocal/mix/mix_clustering_c0.60.tsv", header = T, sep = '\t', comment.char = "")
colnames(gcfIDs) <- c("Id", "GCFnum")
finalatt2 <- merge(finalatt, gcfIDs, by = "Id")

write.table(finalatt2, file = "bgc_cluster_attributes.csv", row.names = F, sep = ',')
write.table(finalatt2, file = "bgc_cluster_attributes.tsv", row.names = F, sep = '\t', quote = F)
