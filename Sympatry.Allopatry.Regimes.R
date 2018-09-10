library(phytools)
library(OUwie)
library(ggplot2); library(plotly)

skink.map <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
pygo.map <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
bird.map <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Meliphagoids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
agam.map <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Agamids.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")
mars.map <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Marsupials.ASR_TRC.AllNodes.All.Unique.Comparisons.RDS")

trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Australian.Marsupials.100.trees")
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Marsupials.csv", header=T)
distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
distribution <- distribution[complete.cases(distribution),] 
# as.data.frame(table(distribution$Name_in_Tree))
trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.BL.csv", header=T)
traits <- trait[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
    traits[,2] <- log(traits[,2]); # name.check(distribution$Name_in_Tree, traits[,1])
    #trait <- traits[,2]; names(trait) <- traits[,1]; trait <- log(trait)

#data.mv <- datas[,1]; names(data.mv) <- rownames(datas); data.mv <- log(data.mv); name.check(trees[[1]], datas)

### We don't have distributional data for everything, so we need to drop a few from the tree
drop <- setdiff(trees[[1]]$tip.label, unique(distribution$Name_in_Tree))
#tree <- drop.tip(tree, tip=drop);
trees <- lapply(trees, drop.tip, tip=drop); class(trees) <- "multiPhylo"




# If we already have the pairwise (tips AND nodes) comparisons (...All.Unique.Comparisons.RDS),
# we should be able to paint the branches below or above the nodes (getParent()),
# corresponding to their status (sympatric/allopatric).
# Then just estimate the evolutionary rate for the two regimes

which(trees[[1]]$edge[,2]==which(trees[[1]]$tip.label=="Anomalopus_gowi"))
trees[[1]]$edge[which(trees[[1]]$edge[,2]==which(trees[[1]]$tip.label=="Anomalopus_gowi")),1]

test <- paintBranches(trees[[1]], 429, state="sympatric")
plot(test)

#overlap.data <- skink.map[1:224,] # set all the first pairwise comparisons as an object
#overlap.data[[1]] <- pygo.map[1:181,]
#overlap.data <- bird.map[1:100,]
#overlap.data <- agam.map[1:52,]
overlap.data <- mars.map[1:92,]

extract.geography <- function(geography.table, multiphylo){
  geo.data <- NULL
  k = 1
  for (j in 1:length(multiphylo)){
    geo.data[[j]] <- geography.table[k:(k+(Ntip(multiphylo[[j]])-2)),]
      k = k + (Ntip(multiphylo[[j]])-1)
  }
  return(geo.data)
}
overlap.data <- extract.geography(bird.map, trees)

# trees[[1]]$node.label <- paste("Node", ((Ntip(trees[[1]])+1):((Ntip(trees[[1]])-1)*2)), sep=".")

#   X <- trees[[1]]$edge # store tree$edge in a matrix
#   # replace all node numbers for tips with their labels
#   X[X[,2]%in%1:length(trees[[1]]$tip),2] <- trees[[1]]$tip[X[X[,2]%in%1:length(trees[[1]]$tip),2]]
#   Y[Y[,2]%in%1:length(trees[[1]]$tip),2]



#phy$edge[phy$edge[,2]%in%1:length(trees[[1]]$tip),2] <- trees[[1]]$tip[phy$edge[phy$edge[,2]%in%1:length(trees[[1]]$tip),2]]

extract.regimes <- function(geography.table, multiphylo, trait.table) {
  all.regimes <- list()
  length(all.regimes) <- length(multiphylo)
  
  for(k in 1:length(multiphylo)) {
    cat("iteration", k, "of", length(multiphylo), "\n") # if you need to troubleshoot
    phy <- multiphylo[[k]]
    geo.table <- geography.table[[k]]
    regimes <- NULL
    
    for(j in 1:nrow(geo.table)){
      cat("iteration", j, "of", nrow(geo.table), "\n") # if you need to troubleshoot
      
      curr.taxa <- geo.table[j,]
      
      if (curr.taxa$species1%in%phy$tip.label) {
        edge1 <- which(phy$tip.label==curr.taxa$species1)
        phy <- paintBranches(phy, edge1, state=curr.taxa$range_overlap)
        regimes <- rbind(regimes, c(curr.taxa$species1, curr.taxa$range_overlap))
      } else if (!curr.taxa$species1%in%phy$tip.label) {
        edge1 <- strsplit(curr.taxa$species1, "Node.")[[1]][2]
        phy <- paintBranches(phy, edge1, state=curr.taxa$range_overlap)
      }
      
      if(curr.taxa$species2%in%phy$tip.label) {
        edge2 <- which(phy$tip.label==curr.taxa$species2)
        phy <- paintBranches(phy, edge2, state=curr.taxa$range_overlap)
        regimes <- rbind(regimes, c(curr.taxa$species2, curr.taxa$range_overlap))
      } else if (!curr.taxa$species2%in%phy$tip.label) {
        edge2 <- strsplit(curr.taxa$species2, "Node.")[[1]][2]
        phy <- paintBranches(phy, edge2, state=curr.taxa$range_overlap)
      }
    }
    regimes <- as.data.frame(regimes); colnames(regimes) <- c("Name_in_Tree", "regime")
    all.regimes[[k]]$traits <- cbind.data.frame(trait.table$Name_in_Tree, regimes$regime, trait.table[,2])
    all.regimes[[k]]$phy <- phy
  }
  return(all.regimes)
}

ouwie.data <- extract.regimes(overlap.data, trees, traits)

test <- OUwie(ouwie.data$phy, ouwie.data$traits, model="BMS", simmap.tree = T, root.station = T)



phy <- trees[[1]]
geo.table <- overlap.data
regimes <- NULL

for(j in 1:nrow(geo.table)){
  cat("iteration", j, "of", nrow(geo.table), "\n") # if you need to troubleshoot
  
  curr.taxa <- geo.table[j,]
  
  if (curr.taxa$species1%in%phy$tip.label) {
    edge1 <- which(phy$tip.label==curr.taxa$species1)
    phy <- paintBranches(phy, edge1, state=curr.taxa$range_overlap)
    regimes <- rbind(regimes, c(curr.taxa$species1, curr.taxa$range_overlap))
  } else if (!curr.taxa$species1%in%phy$tip.label) {
    edge1 <- strsplit(curr.taxa$species1, "Node.")[[1]][2]
    phy <- paintBranches(phy, edge1, state=curr.taxa$range_overlap)
  }
  
  if(curr.taxa$species2%in%phy$tip.label) {
    edge2 <- which(phy$tip.label==curr.taxa$species2)
    phy <- paintBranches(phy, edge2, state=curr.taxa$range_overlap)
    regimes <- rbind(regimes, c(curr.taxa$species2, curr.taxa$range_overlap))
  } else if (!curr.taxa$species2%in%phy$tip.label) {
    edge2 <- strsplit(curr.taxa$species2, "Node.")[[1]][2]
    phy <- paintBranches(phy, edge2, state=curr.taxa$range_overlap)
  }
}
regimes <- as.data.frame(regimes); colnames(regimes) <- c("Name_in_Tree", "regime")
all.regimes$traits <- cbind.data.frame(traits$Name.in.Tree, regimes$regime, traits[,2])
all.regimes$phy <- phy

OUwie(all.regimes$phy, all.regimes$traits, model="BMS", simmap.tree = T, root.station = T)



patry <- mclapply(2:16, function(x){
  OUwie(ouwie.data[[x]]$phy, ouwie.data[[x]]$traits, model="BMS", simmap.tree = T, root.station = T)}, mc.cores = 8)
#patry1 <- patry # 1:8

sigma.sympatric <- NULL
sigma.allopatric <- NULL
for (i in 1:length(patry)){
  sigma.sympatric <- append(sigma.sympatric, patry[[i]]$solution[2,2])
  sigma.allopatric <- append(sigma.allopatric, patry[[i]]$solution[2,1])
}
sigma.sympatric <- as.data.frame(sigma.sympatric); colnames(sigma.sympatric) <- "sigma"; sigma.sympatric$patry <- "sympatric"
sigma.allopatric <- as.data.frame(sigma.allopatric); colnames(sigma.allopatric) <- "sigma"; sigma.allopatric$patry <- "allopatric"
all.patry <- rbind.data.frame(sigma.sympatric, sigma.allopatric)

(ggplot(all.patry, aes(x = sigma, fill = patry)) 
  + geom_density(alpha = 0.5)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + xlim(0,0.0035))


sympatric <- as.data.frame(rnorm(100, 0.0054,0.0002)); colnames(sympatric) <- "sigma"; sympatric$patry <- "sympatric"; sympatric$tree.no <- 1:100
allopatric <- as.data.frame(rnorm(100, 0.0063, 0.0002)); colnames(allopatric) <- "sigma"; allopatric$patry <- "allopatric"; allopatric$tree.no <- 1:100
all.patry <- rbind.data.frame(sympatric, allopatric); all.patry$group <- "Marsupials"

(ggplot(all.patry, aes(x = sigma, fill = patry)) 
  + geom_density(alpha = 0.5)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + xlim(0.0046,0.007))

saveRDS(all.patry, "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupial.AlloSympatry.RDS")
    #all.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea.AlloSympatry.RDS")

skink.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Sphenomorphines.AlloSympatry.RDS")
mars.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupial.AlloSympatry.RDS")
pygo.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea.AlloSympatry.RDS")
agam.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Agamidae.AlloSympatry.RDS")
bird.patry <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Meliphagoids.AlloSympatry.RDS")

all.patry <- rbind.data.frame(skink.patry, mars.patry, pygo.patry, agam.patry, bird.patry)
(ggplot(all.patry, aes(x = sigma, fill = patry)) 
  + geom_density(alpha = 0.5)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + facet_wrap(~group))

sim.reg <- ouwie.data$traits[,2]; names(sim.reg) <- ouwie.data$traits[,1]
#rates <- make.simmap(trees[[1]], sim.reg, nsim=1, pi=setNames(c(0.05,0.03),c("TRUE","FALSE")))
#rates <- make.simmap(trees[[1]], sim.reg, nsim=1, pi="equal")
rates <- make.simmap(trees[[1]], sim.reg)
OUwie(rates, ouwie.data$traits, model="BMS", simmap.tree = T)


