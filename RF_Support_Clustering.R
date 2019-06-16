library(phytools)
library(RPANDA)

enviro.data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Data/Andrae_S1.csv", header=T)
flux.data <-   read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Data/Aeolian_Flux.csv", header=T)
data(InfTemp)



library(deeptime); library(gridExtra)

pp <- ggplot(enviro.data, aes(Age)) +
  geom_ribbon(aes(ymin = C4_recon_lower, ymax = C4_recon_upper), fill = "pink") +
  geom_line(aes(y = C4_recon_mean), color="red") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0,80), expand = FALSE) 

qq <- gggeo_scale(pp, dat="epochs")


rr <- ggplot(flux.data, aes(Age)) +
  geom_ribbon(aes(ymin = A_Flux-35, ymax = A_Flux+35), fill = "light blue") +
  geom_line(aes(y = A_Flux), color="blue") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 13), ylim = c(0,150), expand = FALSE) 

ss <- gggeo_scale(rr, dat="epochs")

grid.arrange(qq, ss, nrow=4, ncol=2)

sampled.trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/MODEL110_Sampled_Run2_Macropodinae.trees")
hypsodonty.index <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/CrownHeight_Macropodinae_spMEANS.csv", header=T)


# Trim tree and data down to overlapping taxa (for SAMPLED TREES)
overlaps <-   intersect(sampled.trees[[1]]$tip.label, hypsodonty.index$Taxon)
tip.drops <-    setdiff(sampled.trees[[1]]$tip.label, overlaps)
sampled.trees <- lapply(sampled.trees, drop.tip, tip=tip.drops)
sampled.data <- filter(hypsodonty.index, Taxon %in% overlaps)
sampled.HI <- sampled.data[,2]; names(sampled.HI) <- sampled.data[,1]; geiger::name.check(sampled.trees[[1]], sampled.HI) 

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/MDS_Clustering_Source.R")
class(sampled.trees) <- "multiPhylo"
inputs <- unroot(sampled.trees)
initial <- topclustMDS(inputs, mdsdim=2, max.k=10, makeplot = T, criterion = "max", trdist = "score")

##### We want to be able to compare differences among trees to differences in model fit
##### so we'll use a few different metrics to try and get at what causes discrepancies

# Calculate the Pairwise Robinson-Foulds distances among all trees
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/MDS_Clustering_Source.R")
testo <- pairwise.RF(sampled.trees, measure="score")

# Or alternatively use the Quartet Dissimilarity method which is more sensitive
library(Quartet)
pairwise.TQ <- function(input.trees){
  pTQ <- as.data.frame(TQDist(input.trees)); colnames(pTQ) <- rownames(pTQ)
  TQtable <- NULL
  for(i in 1:nrow(pTQ)){
    for(j in i:ncol(pTQ)){
      TQtable <- rbind(TQtable, data.frame(tree1=i, tree2=j, TQdist=pTQ[i,j],
                                           AGEdist=abs(max(nodeHeights(input.trees[[i]])) - max(nodeHeights(input.trees[[j]])))))
    }
  }
  return(list(TQmatrix=pTQ, TQtable=TQtable))
}
chub <- pairwise.TQ(sampled.trees)

# From the model fitting data, extract differences in model preference among trees
all.aics <- read.RDS("/Users/Ian/Desktop/SAMPLED500_Model_Fitting_AICCs.RDS")
pairwise.AIC <- function(input.AIC, target.model, comparison.table){
  AICtable <- NULL
  #timer <- progress_estimated(length(tree.span))
  
  #for(i in 1:length(unique(input.AIC$tree))){
  #  for(j in 1:length(unique(input.AIC$tree))){
  #    rowi <- filter(input.AIC, tree==i & model==target.model)
  #    rowj <- filter(input.AIC, tree==j & model==target.model)
  #    AICtable <- rbind(AICtable, data.frame(tree1=i, tree2=j, 
  #                                           AICCWdiff=abs(rowi$aiccw - rowj$aiccw)))
  #  }
  #  print(paste("finished column", i))
  #}
  
  for(i in 1:nrow(comparison.table)){
    rowi <- filter(input.AIC, tree==comparison.table[i,"tree1"] & model==target.model)
    rowj <- filter(input.AIC, tree==comparison.table[i,"tree2"] & model==target.model)
    AICtable <- rbind(AICtable, data.frame(tree1=rowi$tree, tree2=rowj$tree, 
                                           AICCWdiff=abs(rowi$aiccw - rowj$aiccw)))
    print(paste("finished row", i))
    #print(timer$tick())
  }
  
  return(AICtable)
}
testa <- pairwise.AIC(all.aics, target.model="FLUXexp", chub$TQtable)

# Now combine the different information together into a single dataframe (remove comparisons of a tree with itself)
besto <- left_join(chub$TQtable, testa); besto <- filter(besto, !tree1==tree2)

# Plot AICcWeight difference as a result of Age difference
agedist <- ggplot(besto, aes(x=AGEdist, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#F21A00", shape=16) +
  geom_smooth(col="black") + theme_classic()

# Plot AICcWeight difference as a result of Quartet Dissimilarity or Robinson-Foulds distances
tqdist <- ggplot(besto, aes(x=TQdist, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#EBCC2A", shape=16) +
  geom_smooth(col="black") + theme_classic()

grid.arrange(agedist, tqdist, nrow=1)

# Plot AICcWeight difference as a result of Quartet Dissimilarity or Robinson-Foulds distances
ggplot(besto, aes(x=TQAGE, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#F21A00") +
  geom_smooth() + theme_classic()

filter(besto, TQdist==0)

grassexp <- filter(all.aics, model=="FLUXexp")
ggplot(grassexp, aes(x=age, y=aiccw)) +
  geom_point() +
  geom_smooth() + theme_classic()

# extract trees which had high support for the FLUXexp model
hi.AICS <- filter(all.aics, model=="FLUXexp" & aiccw >= 0.5)
hi.trees <- sampled.trees[hi.AICS$tree]; class(hi.trees) <- "multiPhylo"
hi.chub <- filter(chub$TQtable, tree1 %in% hi.AICS$tree & tree2 %in% hi.AICS$tree)
hi.testa <- pairwise.AIC(hi.AICS, target.model="FLUXexp", hi.chub)
hi.besto <- left_join(hi.chub, hi.testa); hi.besto <- filter(hi.besto, !tree1==tree2)


gvf <- filter(all.aics, model=="GRASSexp" | model=="FLUXexp")
fiftyAIC <- filter(gvf, aiccw>=0.50)

(ggplot(fiftyAIC, aes(x=age, colour=model))
                    #+ geom_density(alpha=0.75, adjust=0.5)
                    #+ geom_histogram(alpha=0.75)
                    + geom_freqpoly(bins=3)
                    + theme_classic()
                    + scale_fill_manual(values=wes_palette("Zissou1", type="continuous", 3))
                    + scale_x_reverse(lim=c(15,5)))


