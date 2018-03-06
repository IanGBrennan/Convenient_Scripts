## Source code for data simulation functions:
  ## TRC
  ## SRC
  ## BM
  ## OU



library(phytools)
library(DataCombine)
library(sp)
library(adehabitatHR)
library(rgeos)
library(rworldmap); library(ggmap)
library(Rmisc)
library(diversitree)
library(mvtnorm)
library(dplyr)
library(phangorn)

#### Read in your distribution of trees, or a single tree
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")
#### Read in your empirical trait values
trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.logMASS.csv", header=F)
#### Read in your distribution data (if valuable, ours helps trim the tree/traits)
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Meliphagoids.csv", header=T)
  distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
    distribution <- distribution[complete.cases(distribution),] 
# as.data.frame(table(distribution$Name_in_Tree))

### We don't have distributional data for everything, so we need to drop a few from the tree
drop <- setdiff(trees[[1]]$tip.label, unique(distribution$Name_in_Tree))
trees <- lapply(trees, drop.tip, tip=drop); class(trees) <- "multiPhylo"
### and from the trait data
trait <- trait[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
traits <- trait[,2]; names(traits) <- trait[,1] #read in data file in RPANDA format
  ## or read in trait data you created previously
  ts <- readRDS("/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoid.Traits.RDS")
    traits <- ts[[1]][1:Ntip(trees[[1]])]


source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R");

## Edit the 'split.vcv' function to allow internal nodes
########################################################
split.vcv.internal <- function(phy, time) {
  mat <- vcvPhylo(phy)
  mat1 <- mat
  mat2 <- mat
  n <- nrow(mat)
  root <- max(mat)
  shift.from.root <- root - time
  
  make.mat.1 <- function(x, shift.from.root) {
    if(x == 0) {
      return(0)
    }
    if(x<shift.from.root) {
      return(x)
    } 
    if(x > shift.from.root){
      return(shift.from.root)
    }
  }
  
  make.mat.2 <- function(x, time, shift.from.root) {
    if(x == 0) {
      return(0)
    }
    
    if(x < shift.from.root) {
      return(0)
    } else{
      return(x-shift.from.root)
    }
  }
  
  mat1 <- matrix(sapply(mat1, make.mat.1, shift.from.root =  shift.from.root), nrow = n, ncol = n, byrow=T)
  
  diag1 <- diag(mat);
  diag.foo <- function(x) {
    if(x<time) {
      return(x)
    } 
    if(x > time){
    } 
  }
  mat2 <- matrix(sapply(mat2, make.mat.2, time =  time, shift.from.root = shift.from.root), nrow = n, ncol = n, byrow=T)
  
  rownames(mat1) <- rownames(mat2) <- rownames(mat)
  colnames(mat1) <- colnames(mat2) <- colnames(mat)
  
  return(list(mat1 = mat1, mat2 = mat2))
}
########################################################


# trees - multiphylo object
# traits - vector, each object has name of a taxon in tree
# model - model for generating the data (TRC, SRC, BM, OU)
# parameters - sig2, optional: alpha, sig2, etc
# shift.time - time to split the VCV matrix in time-shifting models



backwards.sim <- function(trees, traits, model = c("BM", "OU", "TRC", "SRC"),
                          sigma1, sigma2, alpha, shift.time.min, shift.time.max, iterations) {
                          #parameters = c(sigma1, sigma2, alpha, shift.time),
                          #iterations) {
  sim.traits <- NULL
  
  for (rr in 1:iterations) {
    traits.geiger <- NULL; # make intermediate trait lists for each tree in geiger and ouwie data format
    phy <- trees[[rr]] # designating the target tree
    traitz <- NULL; #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
    cat("iteration", rr, "of", iterations, "\n") #keep track of what tree/loop# we're on
    
    if (model == "BM") {
      m <- vcvPhylo(phy)
      m.rev.rel.rad <- m * sigma1
    } else if (model == "OU") {
      m <- vcvPhylo(phy)
      m <- ouMatrix(m, alpha)
      m.rev.rel.rad <- m * sigma1
      #m <- m * sigma1
      #m.rev.rel.rad <- ouMatrix(m, alpha)
    } else if (model == "TRC") {
      shift.time <- sample((seq(from=shift.time.min, to=shift.time.max, by=1)), size=1)
      m <- split.matrices <- split.vcv.internal(phy, shift.time) # divide the vcv matrix at a static time
      m[[1]] <- m[[1]] * sigma1 # adjust BM (old era) vcv according to a rate scalar (usually = 1)
      m[[2]] <- ouMatrix(m[[2]], alpha) # transform the second half of the vcv matrix according to your alpha
      m[[2]] <- m[[2]] * sigma2 # adjust OU (new era) vcv according to a rate scalar (faster or slower)
      m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
    } else if (model == "SRC") {
      shift.time <- sample((seq(from=shift.time.min, to=shift.time.max, by=1)), size=1)
      m <- split.matrices <- split.vcv.internal(phy, shift.time) # divide the vcv matrix at a static time
      m[[1]] <- m[[1]] * sigma1 # adjust BM (old era) vcv according to a rate scalar (usually = 1)
      m[[2]] <- ouMatrix(m[[2]], alpha) # transform the second half of the vcv matrix according to your alpha
      m[[2]] <- m[[2]] * sigma2 # adjust OU (new era) vcv according to a rate scalar (faster or slower)
      m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
    }
    
    # name the internal nodes of the tree
    rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")
    colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")
    
    # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
    traitz <- traits
    node.numbers <- as.data.frame(((Ntip(phy)+1):((Ntip(phy)*2)-1)))
    node.numbers[,2] <- node.numbers[,1]
    node.numbers[,1] <- sapply(node.numbers[,1], function(x) paste("Node", x, sep="."))
    colnames(node.numbers) <- c("name", "number")
    tip.numbers <- as.data.frame(phy$tip.label)
    tip.numbers[,2] <- rownames(tip.numbers)
    colnames(tip.numbers) <- c("name", "number")
    node.numbers <- rbind(node.numbers, tip.numbers)
    
    all.matches <- NULL
    for(p in 1:length(unique(rownames(m.rev.rel.rad)))) {
      taxon <- unique(rownames(m.rev.rel.rad))[p] # choose a current taxon
      
      if (p>(Ntip(phy))) {
        sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
        pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
        pair.matrix[,1] <- taxon # who they're being compared against (taxon)
        if (sisters<=Ntip(phy)) {
          pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters]
        } else {
          pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters-1] # add all the comparisons
        }
      } else {
        sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
        pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
        pair.matrix[,1] <- taxon # who they're being compared against (taxon)
        if (sisters<Ntip(phy)) {
          pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters]
        } else {
          pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters-1] # add all the comparisons
        }
      }
      all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
    }
    
    single.matches <- NULL
    skippers <- NULL
    for(i in 1:nrow(all.matches)){
      curr.guy <- unique(all.matches[i,1])
      if(!curr.guy %in% skippers){
        x <- rbind(all.matches[which(all.matches[,1] %in% curr.guy),], all.matches[which(all.matches[,2] %in% curr.guy),])
        skippers <- rbind(skippers,x[which(!x[,1] %in% curr.guy),1])
        single.matches <- rbind(single.matches,x[1,])
      }
    }
    single.matches <- as.data.frame(single.matches)
    
    tip.tip  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & single.matches[,2] %in% phy$tip.label)
    tip.tip <- as.matrix(tip.tip)
    tip.node  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & !(single.matches[,2] %in% phy$tip.label))
    node.node <- filter(single.matches, !(single.matches[,1] %in% phy$tip.label) & !(single.matches[,2] %in% phy$tip.label))
    
    tips.nodes <- rbind(tip.node, node.node)
    for (uu in 1:nrow(tips.nodes)) {
      tips.nodes[uu, "max_depth"] <- (max(nodeHeights(phy))-nodeheight(phy, getParent(phy, filter(node.numbers, name==tips.nodes[uu,2])[,2])))
    }
    tips.nodes <- tips.nodes[order(tips.nodes$max_depth),]
    tips.nodes <- as.matrix(tips.nodes)
    
    for (jj in 1:nrow(tip.tip)) {
      sister.tips <- tip.tip[jj,]
      parent.mean <- NULL; parent.estimate <- NULL
      
      parent.mean <- mean(c(traits[sister.tips[1]], traits[sister.tips[2]]))
      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
      parent.node <- getParent(phy, which(rownames(m.rev.rel.rad)==sister.tips[1]))
      names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1] 
      
      traitz <- append(traitz, parent.estimate)
    }
    
    for (jj in 1:nrow(tips.nodes)) {
      sister.tips <- tips.nodes[jj,]
      parent.mean <- NULL; parent.estimate <- NULL
      
      if (sister.tips[1] %in% phy$tip.label && !(sister.tips[2] %in% phy$tip.label)) {
        parent.mean <- mean(c(traitz[sister.tips[1]], traitz[sister.tips[2]]))
        parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
        parent.node <- getParent(phy, which(rownames(m.rev.rel.rad)==sister.tips[1]))
        names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1] 
      }
      if (!(sister.tips[1] %in% phy$tip.label) && !(sister.tips[2] %in% phy$tip.label)) {
        parent.mean <- mean(c(traitz[sister.tips[1]], traitz[sister.tips[2]]))
        parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
        parent.node <- which(rownames(m.rev.rel.rad)==sister.tips[1])
        if (parent.node == Ntip(phy)+1) {
          names(parent.estimate) <- paste("Node", parent.node, sep=".")
        } else {
          names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1]
        }
      }
      
      if (is.na(parent.estimate)) {
        redo.list <- rbind(redo.list, tips.nodes[jj,])
      }
      
      traitz <- append(traitz, parent.estimate)
    }
    sim.traits[[rr]] <- traitz
  }
  return(sim.traits)
}

simulated.traits <- backwards.sim(trees, traits, model="TRC", sigma1=0.01, sigma2=0.2, 
                                  shift.time.min=5, shift.time.max=9, alpha=1,  iterations=100)

saveRDS(simulated.traits, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoids.Traits.RDS")

### If you'd like to visualize the simulated traits, plot below
sims <- simulated.traits[[1]]
s <- names(sims)[(Ntip(trees[[1]])+1) : length(sims)]
s1 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[2]))
names(sims)[(Ntip(trees[[1]])+1) : length(sims)] <- s1
phenogram(trees[[1]], sims, fsize=0.5, lwd=0.5, spread.labels=T)
#































sim.traits <- NULL
for (z in 1:1) {
  traits.geiger <- NULL; # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- trees[[z]] # designating the target tree
  traitz <- NULL; #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  
  # Option A (comment out top when simulating different shift times, comment out bottom when simulating different alphas)
  m <- split.matrices <- split.vcv.internal(phy, sim.shifts[z]) # divide the vcv matrix at a static time
  #m <- split.matrices <- split.vcv(phy, sim.shifts[i]) # or at differing times
  
  # Option B (adjust to change simulating model, comment out both to get the SRC model)
  m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
  m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
  
  # Option C (comment out top when simulating different times, comment out bottom when simulating different alphas)
  m[[2]] <- ouMatrix(m[[2]], alpha=diff.alpha[z]) # transform the second half of the vcv matrix according to your alpha
  #m[[2]] <- ouMatrix(m[[2]], alpha=alpha) # transform the second half of the vcv matrix according to your alpha
  
  m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
  rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")
  colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")
  
  # OR, do it like the 'ecological release' model
  # m <- lapply(m, function(x) x*0.1)
  # m.rev.rel <- m[[1]] + m[[2]]
  
  # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
  traitz <- traits
  node.numbers <- as.data.frame(((Ntip(phy)+1):((Ntip(phy)*2)-1)))
  node.numbers[,2] <- node.numbers[,1]
  node.numbers[,1] <- sapply(node.numbers[,1], function(x) paste("Node", x, sep="."))
  colnames(node.numbers) <- c("name", "number")
  tip.numbers <- as.data.frame(phy$tip.label)
  tip.numbers[,2] <- rownames(tip.numbers)
  colnames(tip.numbers) <- c("name", "number")
  node.numbers <- rbind(node.numbers, tip.numbers)
  
  all.matches <- NULL
  for(p in 1:length(unique(rownames(m.rev.rel.rad)))) {
    taxon <- unique(rownames(m.rev.rel.rad))[p] # choose a current taxon
    
    if (p>(Ntip(phy))) {
      sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<=Ntip(phy)) {
        pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters]
      } else {
        pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters-1] # add all the comparisons
      }
    } else {
      sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<Ntip(phy)) {
        pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters]
      } else {
        pair.matrix[,2] <- rownames(m.rev.rel.rad)[sisters-1] # add all the comparisons
      }
    }
    all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
  }
  
  single.matches <- NULL
  skippers <- NULL
  for(i in 1:nrow(all.matches)){
    curr.guy <- unique(all.matches[i,1])
    if(!curr.guy %in% skippers){
      x <- rbind(all.matches[which(all.matches[,1] %in% curr.guy),], all.matches[which(all.matches[,2] %in% curr.guy),])
      skippers <- rbind(skippers,x[which(!x[,1] %in% curr.guy),1])
      single.matches <- rbind(single.matches,x[1,])
    }
  }
  single.matches <- as.data.frame(single.matches)
  
  tip.tip  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & single.matches[,2] %in% phy$tip.label)
  tip.tip <- as.matrix(tip.tip)
  tip.node  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & !(single.matches[,2] %in% phy$tip.label))
  node.node <- filter(single.matches, !(single.matches[,1] %in% phy$tip.label) & !(single.matches[,2] %in% phy$tip.label))
  
  tips.nodes <- rbind(tip.node, node.node)
  for (uu in 1:nrow(tips.nodes)) {
    tips.nodes[uu, "max_depth"] <- (max(nodeHeights(phy))-nodeheight(phy, getParent(phy, filter(node.numbers, name==tips.nodes[uu,2])[,2])))
  }
  tips.nodes <- tips.nodes[order(tips.nodes$max_depth),]
  tips.nodes <- as.matrix(tips.nodes)
  
  for (jj in 1:nrow(tip.tip)) {
    sister.tips <- tip.tip[jj,]
    parent.mean <- NULL; parent.estimate <- NULL
    
    parent.mean <- mean(c(traits[sister.tips[1]], traits[sister.tips[2]]))
    parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
    parent.node <- getParent(phy, which(rownames(m.rev.rel.rad)==sister.tips[1]))
    names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1] 
    
    traitz <- append(traitz, parent.estimate)
  }
  
  for (jj in 1:nrow(tips.nodes)) {
    sister.tips <- tips.nodes[jj,]
    parent.mean <- NULL; parent.estimate <- NULL
    
    if (sister.tips[1] %in% phy$tip.label && !(sister.tips[2] %in% phy$tip.label)) {
      parent.mean <- mean(c(traitz[sister.tips[1]], traitz[sister.tips[2]]))
      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
      parent.node <- getParent(phy, which(rownames(m.rev.rel.rad)==sister.tips[1]))
      names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1] 
    }
    if (!(sister.tips[1] %in% phy$tip.label) && !(sister.tips[2] %in% phy$tip.label)) {
      parent.mean <- mean(c(traitz[sister.tips[1]], traitz[sister.tips[2]]))
      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
      parent.node <- which(rownames(m.rev.rel.rad)==sister.tips[1])
      if (parent.node == Ntip(phy)+1) {
        names(parent.estimate) <- paste("Node", parent.node, sep=".")
      } else {
        names(parent.estimate) <- rownames(m.rev.rel.rad)[parent.node-1]
      }
    }
    
    if (is.na(parent.estimate)) {
      redo.list <- rbind(redo.list, tips.nodes[jj,])
    }
    
    traitz <- append(traitz, parent.estimate)
  }
  
  sim.traits[[z]] <- traitz  
  saveRDS(sim.traits, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/ASR_TRC.Pygopodoids.Traits.RDS")
}

## BELOW THIS IS JUNK
#################################################




rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", rownames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")
colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])] <- paste("Node", colnames(m.rev.rel.rad)[(length(phy$tip.label)+1):length(m.rev.rel.rad[,1])], sep = ".")

outlist <- list()
for (oo in 1:nrow(m.rev.rel.rad)) {
  for (uu in 1:ncol(m.rev.rel.rad)) {
    #value <- setNames(rmvnorm(n=1, mean=mean(traits[m.rev.rel.rad[oo]], traits[m.rev.rel.rad[[uu]]]), sigma=m.rev.rel.rad[oo,uu]), rownames(m.rev.rel.rad))
    value <- mvrnorm(n=1, mu=mean(traits[rownames(m.rev.rel.rad)[oo]], traits[colnames(m.rev.rel.rad)[uu]]), Sigma=m.rev.rel.rad[oo,uu])
    node.num1 <- which(trees[[1]]$tip.label==rownames(m.rev.rel.rad)[oo])
    node.num2 <- which(trees[[1]]$tip.label==colnames(m.rev.rel.rad)[uu])
    if (node.num1 == node.num2) {
      name.it <- rownames(m.rev.rel.rad)[oo]
    } else {
      parent.node <- getMRCA(trees[[1]], c(node.num1, node.num2))
      name.it <- rownames(m.rev.rel.rad)[parent.node]
    }
    names(value) <- name.it # HOW DO I NAME THE OUTPUT VALUE???!!
    outlist <- append(outlist, value)
  }
}




traits.list <- NULL
for (oo in 1:length(trees)) {
  cat("iteration", oo, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  all.traits <- traits
  phy <- trees[[oo]]
  #phy$node.label <- paste("Node", ((Ntip(phy)+1):((Ntip(phy)-1)*2)), sep=".")
  #plot(phy, cex=0.3, show.node.label=T)
  node.numbers <- as.data.frame(((Ntip(phy)+1):((Ntip(phy)*2)-1)))
  node.numbers[,2] <- node.numbers[,1]
  node.numbers[,1] <- sapply(node.numbers[,1], function(x) paste("Node", x, sep="."))
  colnames(node.numbers) <- c("name", "number")
  
  all.matches <- NULL
  for(p in 1:length(unique(rownames(dist.mat)))) {
    taxon <- unique(rownames(dist.mat))[p] # choose a current taxon
    
    if (p>(Ntip(phy))) {
      sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<=Ntip(phy)) {
        pair.matrix[,2] <- rownames(dist.mat)[sisters]
      } else {
        pair.matrix[,2] <- rownames(dist.mat)[sisters-1] # add all the comparisons
      }
    } else {
      sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<Ntip(phy)) {
        pair.matrix[,2] <- rownames(dist.mat)[sisters]
      } else {
        pair.matrix[,2] <- rownames(dist.mat)[sisters-1] # add all the comparisons
      }
    }
    all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
  }
  
  single.matches <- NULL
  skippers <- NULL
  for(i in 1:nrow(all.matches)){
    curr.guy <- unique(all.matches[i,1])
    if(!curr.guy %in% skippers){
      x <- rbind(all.matches[which(all.matches[,1] %in% curr.guy),], all.matches[which(all.matches[,2] %in% curr.guy),])
      skippers <- rbind(skippers,x[which(!x[,1] %in% curr.guy),1])
      single.matches <- rbind(single.matches,x[1,])
    }
  }
  single.matches <- as.data.frame(single.matches)
  #max_depth <- NULL
  #for (pp in 1:nrow(single.matches)) {
  #  single.matches[pp, "max_depth"] <- (max(nodeHeights(phy))-nodeheight(phy, which(rownames(dist.mat)==single.matches[pp,2])+1))
  #}
  #single.matches[,"max_depth"] <- lapply(single.matches, )
  
  tip.tip  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & single.matches[,2] %in% phy$tip.label)
  tip.tip <- as.matrix(tip.tip)
  tip.node  <- filter(single.matches, single.matches[,1] %in% phy$tip.label & !(single.matches[,2] %in% phy$tip.label))
  node.node <- filter(single.matches, !(single.matches[,1] %in% phy$tip.label) & !(single.matches[,2] %in% phy$tip.label))
  
  tips.nodes <- rbind(tip.node, node.node)
  for (uu in 1:nrow(tips.nodes)) {
    #parent.node <- # NEED TO SORT OUT THE AGE OF THE PARENT NODE
    #tips.nodes[uu, "max_depth"] <- (max(nodeHeights(phy))-nodeheight(phy, which(rownames(dist.mat)==tips.nodes[uu,2])+1))
    tips.nodes[uu, "max_depth"] <- (max(nodeHeights(phy))-nodeheight(phy, getParent(phy, filter(node.numbers, name==tips.nodes[uu,2])[,2])))
  }
  tips.nodes <- tips.nodes[order(tips.nodes$max_depth),]
  tips.nodes <- as.matrix(tips.nodes)
  
  
  for (jj in 1:nrow(tip.tip)) {
    sister.tips <- tip.tip[jj,]
    parent.mean <- NULL; parent.estimate <- NULL
    
    parent.mean <- mean(c(traits[sister.tips[1]], traits[sister.tips[2]]))
    parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
    parent.node <- getParent(phy, which(rownames(dist.mat)==sister.tips[1]))
    names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
    
    #intermediate <- append(intermediate, parent.estimate)
    all.traits <- append(all.traits, parent.estimate)
  }
  
  for (jj in 1:nrow(tips.nodes)) {
    sister.tips <- tips.nodes[jj,]
    parent.mean <- NULL; parent.estimate <- NULL
    
    if (sister.tips[1] %in% phy$tip.label && !(sister.tips[2] %in% phy$tip.label)) {
      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
      parent.node <- getParent(phy, which(rownames(dist.mat)==sister.tips[1]))
      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
    }
    if (!(sister.tips[1] %in% phy$tip.label) && !(sister.tips[2] %in% phy$tip.label)) {
      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
      parent.node <- which(rownames(dist.mat)==sister.tips[1])
      #parent.node <- getParent(phy, (which(rownames(dist.mat)==sister.tips[1])+1))
      if (parent.node == Ntip(phy)+1) {
        names(parent.estimate) <- paste("Node", parent.node, sep=".")
      } else {
        names(parent.estimate) <- rownames(dist.mat)[parent.node-1]
      }
    }
    
    if (is.na(parent.estimate)) {
      redo.list <- rbind(redo.list, tip.node[jj,])
    }
    
    #intermediate <- append(intermediate, parent.estimate)
    all.traits <- append(all.traits, parent.estimate)
  }
  traits.list[[oo]] <- all.traits
  
  
  
  
  #  redo.list <- NULL
  #  for (jj in 1:nrow(tip.node)) {
  #    sister.tips <- tip.node[jj,]
  #      #sis1 <- sister.tips[,1]; sis2 <- sister.tips[,2]
  #    parent.mean <- NULL; parent.estimate <- NULL
  #    
  #      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
  #      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
  #      parent.node <- getParent(phy, which(rownames(dist.mat)==sister.tips[1]))
  #      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
  #      
  #      if (is.na(parent.estimate)) {
  #        redo.list <- rbind(redo.list, tip.node[jj,])
  #      }
  #    
  #      #intermediate <- append(intermediate, parent.estimate)
  #      all.traits <- append(all.traits, parent.estimate)
  #  }
  #  
  #  for (jj in 1:nrow(node.node)) {
  #    sister.tips <- node.node[jj,]
  #    parent.mean <- NULL; parent.estimate <- NULL
  #
  #      parent.mean <- mean(c(all.traits[sister.tips[,1]], all.traits[sister.tips[,2]]))
  #      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[,1], sister.tips[,2]])
  #      parent.node <- which(rownames(dist.mat)==sister.tips[,1])
  #      #parent.node <- getParent(phy, (which(rownames(dist.mat)==sister.tips[1])+1))
  #      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
  #    
  #      intermediate <- append(intermediate, parent.estimate)
  #      #all.traits <- append(all.traits, parent.estimate)
  #  }
  #  
  #
  #
  #
  #
  #
  #  redo.list <- NULL
  #  for (jj in 1:nrow(single.matches)) {
  #    sister.tips <- single.matches[jj,]
  #    parent.mean <- NULL; parent.estimate <- NULL
  #    if (sister.tips[1] %in% phy$tip.label && !(sister.tips[2] %in% phy$tip.label)) {
  #      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
  #      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
  #      parent.node <- getParent(phy, which(rownames(dist.mat)==sister.tips[1]))
  #      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
  #      if (is.na(parent.estimate)) {
  #        redo.list <- rbind(redo.list, single.matches[jj,])
  #      }
  #    }
  #    all.traits <- append(all.traits, parent.estimate)
  #  }
  #  all.traits <- all.traits[!is.na(all.traits)]
  #  
  #  for (jj in 1:nrow(redo.list)) {
  #    sister.tips <- redo.list[jj,]
  #    parent.mean <- NULL; parent.estimate <- NULL
  #    if (sister.tips[1] %in% phy$tip.label && !(sister.tips[2] %in% phy$tip.label)) {
  #      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
  #      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
  #      parent.node <- getParent(phy, which(rownames(dist.mat)==sister.tips[1]))
  #      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
  #      if (is.na(parent.estimate)) {
  #        redo.list <- rbind(redo.list, single.matches[jj,])
  #      }
  #    }
  #    all.traits <- append(all.traits, parent.estimate)
  #  }
  #  all.traits <- all.traits[!is.na(all.traits)]
  #
  #  for (jj in nrow(single.matches):1) {
  #    sister.tips <- single.matches[jj,]
  #    parent.mean <- NULL; parent.estimate <- NULL
  #    if (!(sister.tips[1] %in% phy$tip.label) && !(sister.tips[2] %in% phy$tip.label)) {
  #      parent.mean <- mean(c(all.traits[sister.tips[1]], all.traits[sister.tips[2]]))
  #      parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[sister.tips[1], sister.tips[2]])
  #      parent.node <- which(rownames(dist.mat)==sister.tips[1])
  #      #parent.node <- getParent(phy, (which(rownames(dist.mat)==sister.tips[1])+1))
  #      names(parent.estimate) <- rownames(dist.mat)[parent.node-1] 
  #    }
  #    all.traits <- append(all.traits, parent.estimate)
  #  }
  #  all.traits <- all.traits[!is.na(all.traits)]
  #  
  #  intermediate <- intermediate[!is.na(intermediate)]
  #  redo.list <- names(intermediate[is.na(intermediate)])
  #  setdiff(rownames(m.rev.rel.rad), names(c(traits, intermediate)))
  #
  #
  #
  #
  #  which(rownames(dist.mat)=="Node.91")
  #
  #  
  #  parent <- getParent(phy, oo)
  #  daughters <- Descendants(phy, parent, type="tips")
  #  daughters <- Descendants(phy, parent, type="children")
  #  
  #  daughter.list <- NULL
  #  for (kk in 1:length(daughters[[1]])) {
  #    daughter.list <- append(daughter.list, traits[daughters[[1]]][kk])
  #  }
  #  parent.mean <- mean(c(daughter.list))
  #  parent.estimate <- mvrnorm(n=1, mu=parent.mean, Sigma=m.rev.rel.rad[parent])
  #  names(parent.estimate) <- parent
  #  all.traits <- append(all.traits, parent.estimate)
}







dist.mat <- max(nodeHeights(phy)) - vcvPhylo(phy)
rownames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])] <- paste("Node", rownames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])], sep = ".")
colnames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])] <- paste("Node", colnames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])], sep = ".")

all.matches <- NULL
#phy <- trees[[1]]
for(p in 1:length(unique(rownames(dist.mat)))) {
  taxon <- unique(rownames(dist.mat))[p] # choose a current taxon
  
  if (p>(Ntip(phy))) {
    sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
    pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
    pair.matrix[,1] <- taxon # who they're being compared against (taxon)
    if (sisters<=Ntip(phy)) {
      pair.matrix[,2] <- rownames(dist.mat)[sisters]
    } else {
      pair.matrix[,2] <- rownames(dist.mat)[sisters-1] # add all the comparisons
    }
  } else {
    sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
    pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
    pair.matrix[,1] <- taxon # who they're being compared against (taxon)
    if (sisters<Ntip(phy)) {
      pair.matrix[,2] <- rownames(dist.mat)[sisters]
    } else {
      pair.matrix[,2] <- rownames(dist.mat)[sisters-1] # add all the comparisons
    }
  }
  all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
}
#


