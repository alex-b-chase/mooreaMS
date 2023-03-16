setwd("/Volumes/JensenLabMGs/alex_alyssa-MooreaMGs/moorea2020/metagenome_processing/MAGs/assemblies/bigscape")

numbgcs <- read.table("numBGCsample.txt", header = F, sep = '\t')
colnames(numbgcs) <- c("numBGC", "sampleID")

metadata <- read.table("../../../mgenomeID.txt", header = T, sep = '\t', comment.char = "")

totaldata <- merge(numbgcs, metadata, by = "sampleID")

library(ggplot2)


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

summtot <- summarySE(totaldata, measurevar = "numBGC", groupvars = c("site", "sitecolor"))

staxcolors <- setNames(summtot$sitecolor, summtot$site)

pdf("bgccountsSite.pdf", height = 12, width = 8)
ggplot(summtot, aes(x = site, y = numBGC, fill = site)) +
  stat_summary(fun = "mean", geom = "bar", aes(colour="black")) +
  geom_errorbar(aes(ymin = numBGC - se, ymax = numBGC + se),
                width = 0,                    # Width of the error bars
                position = position_dodge(.9)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(y = "Number of BGCs", x = "") +
  scale_fill_manual(values = staxcolors, guide = F) +
  scale_color_manual(values = c("black"), guide = F) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
dev.off()
  