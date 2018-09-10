library(phytools)
library(ggplot2)
library(Rmisc)
library(dplyr)

bird <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Meliphagoids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
  bird$clade <- "Meliphagoidea"
    bm.bird <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.BM.Meliphagoids.AllNodes.All.Unique.Comparisons.RDS")
  
agam <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Agamids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
  agam$clade <- "Agamidae"
    bm.agam <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.BM.Agamids.AllNodes.All.Unique.Comparisons.RDS")
    trc.agam <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Agamids.AllNodes.All.Unique.Comparisons.RDS")
mars <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Marsupials.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
  mars$clade <- "Marsupials"
pygo <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
  pygo$clade <- "Pygopodoidea"
skink <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
  skink$clade <- "Skinks"

#########################

all.clades <- rbind(mars, agam, bird, pygo, skink)
all.one <- rbind(pbt[[1]], mbt[[1]], abt[[1]], bbt[[1]], sbt[[1]])


#short.d <-  filter(bird, dist_tree < 23 & dist_tree > 2)
short.d <- filter(agam, dist_tree < 23)
#short.d <- filter(short.d, dist_trait < 0.5)
lines <- (ggplot(short.d, aes(x=dist_tree, y=dist_trait, color=range_overlap))
          #+ geom_point(alpha=0.1)
          + geom_smooth(aes(fill=range_overlap))
          + scale_x_reverse()
          + ggtitle("Agamid Dragons Empirical:
Trait Distance between Overlapping and Non-overlapping Taxa")
          + theme_classic())
#### Plot the 
square <- (ggplot(short.d)
           + geom_density(aes(dist_tree, ..count.., fill=range_overlap, alpha=0.5), position="fill") # or dist_trait to see if the two modes differ
           + scale_x_reverse()
           + scale_y_reverse(position="right")
           + ggtitle("Agamid Dragons SIMULATED BM:
Trend in Mode of Geographic Speciation")
           + theme_classic())


#+ geom_vline(data=means, aes(xintercept=grp.mean, color=range_overlap),
#              linetype="dashed"))
multiplot(lines, square, layout=matrix(c(1,1,2), nrow=1, byrow=T))
#multiplot(lines, square, cols=2)
multiplot(bm.a, bm.a,
          a, a, 
          trc.a, trc.a,
          cols=3)

lines <- (ggplot(shorties, aes(x=dist_tree, y=dist_trait))
  #+ geom_point(alpha=0.1)
  + geom_smooth(aes(fill=clade, color=clade))
  + scale_x_reverse()
  + ggtitle("All Clades:
            Trait Distance between Overlapping and Non-overlapping Taxa")
  + theme_classic())
(ggplot(short.d)
           + geom_density(aes(dist_tree, ..count.., fill=clade, alpha=0.5), position="fill") # or dist_trait to see if the two modes differ
           + scale_x_reverse()
           #+ scale_y_reverse(position="right")
           + ggtitle("All Clades:
                     Trend in Mode of Geographic Speciation")
           + theme_classic())
shorties <- filter(pygo, dist_tree < 23)
dists <- (ggplot(shorties, aes(dist_tree, ..count.., fill=clade, color=clade))
  + geom_density(alpha=0.25)
  + scale_x_reverse()
  + theme_classic()
  + ggtitle("All Clades:
             Trend in Mode of Geographic Speciation")
  + scale_y_continuous(position="right"))

multiplot(lines, dists, layout=matrix(c(1,1,2), nrow=1, byrow=T))


###
bins <- seq(from=0, to=max(agam$dist_tree), by=0.25)
trend <- NULL;
for (k in 1:length(bins)) {
  period <- filter(agam, dist_tree > bins[[k]] & dist_tree < bins[[k]]+0.25)
  overlaps  <- filter(period, range_overlap==T)
  noverlaps <- filter(period, range_overlap==F)
  proportion <- (nrow(noverlaps))/(nrow(period))
  trend <- rbind(trend, proportion)
}
rownames(trend) <- bins
write.csv(trend, file="/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Allopatry.Trend.Agamids.csv",
          na="0")

############################
# this is going to be ugly, but if we read in the estimates from biogeobears for each clade, we can summarize the trends

agam.time <- read.delim("/Users/Ian/Google.Drive/R.Analyses/BioGeoBEARS/Data.EMPIRICAL/BGB.Pygopodoidea.UPDATED.Empirical.vj.with.CIs.txt", header=T)
    agam.time <- as.data.frame(agam.time$meanCI)
        agam.time$time <- (0:(nrow(agam.time)-1))/10
        
(ggplot(agam.time, aes(x=time, y=agam.time$meanCI))
  + geom_smooth())


test.stack$cent <- test.stack$cent*100
stack.expanded <- test.stack[rep(row.names(test.stack), test.stack$frequency), 2:3]
(ggplot(stack.expanded)
  + geom_density(aes(stack.expanded$time, ..count.., fill=range_overlap, alpha=0.5), position="fill"))

test.stack <- filter(pygo.stack, time<=20)

#