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
library(ggplot2)

# read in your tree, trait data, and rase output data
tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.tre")
trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.MlogBL.csv", header=F); 
    colnames(trait) <- c("Name_in_Tree", "trait")
dist.data <- readRDS("/Users/Ian/Desktop/MioceneAustralia/Marsupials.RASE.RDS")

# determine the taxon overlap among your data sets
#keep <- intersect(tree$tip.label, unique(dist.data$Tip_DistData$Name_in_Tree))
keep <- Reduce(intersect, list(tree$tip.label, unique(dist.data$Tip_DistData$Name_in_Tree), trait$Name_in_Tree))

# trim down the tree, traits, and dist.data (as necessary)
tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
trait <- filter(trait, Name_in_Tree %in% keep)
    trait.data <- trait[,2]; names(trait.data) <- trait[,1] #read in data file in RPANDA format

#### Start by building pairwise trait and tree distance matrices, and a table
build.TTR <- function(phy, traits) {
  
  ## Build a distance matrix for all taxa in the tree
  dist.mat <- max(nodeHeights(phy)) - vcvPhylo(phy)
  
  rownames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])] <- paste0("n", rownames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])])
  colnames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])] <- paste0("n", colnames(dist.mat)[(length(phy$tip.label)+1):length(dist.mat[,1])])
  
  
  ## Build a distance matrix for the trait of interest for all taxa in the tree
  # IF YOU'RE USING EMPIRICAL DATA (use this, not the below)
  ANC <- fastAnc(phy, traits, var=F, CI=F)
  all.traits <- append(traits, ANC)
  
  # Create a matrix of pairwise distances
  diff.mat <- matrix(NA,length(all.traits), length(all.traits))
  colnames(diff.mat) <- names(all.traits); rownames(diff.mat) <- names(all.traits)
  diff.mat <- as.data.frame(diff.mat)
  for (j in 1:length(rownames(diff.mat))) {
    for (k in 1:length(colnames(diff.mat))) {
      diff.mat[j,k] <- abs(all.traits[j]-all.traits[k])
    }
  }
  ## IF USING SIMULATED DATA COMMENT OUT THE RENAMING BELOW!
  rownames(diff.mat)[(length(phy$tip.label)+1):length(diff.mat[,1])] <- paste0("n", rownames(diff.mat)[(length(phy$tip.label)+1):length(diff.mat[,1])])
  colnames(diff.mat)[(length(phy$tip.label)+1):length(diff.mat[,1])] <- paste0("n", colnames(diff.mat)[(length(phy$tip.label)+1):length(diff.mat[,1])])
  
  
  ## Now build a data frame to hold the Tree/Trait/Range distance data
  # First, determine tips in common:
  #rownames(dist.mat)[1:length(tree$tip.label)] <- tree$tip.label
  #rownames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])] <- paste0("n", rownames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])])
  #colnames(dist.mat)[1:length(tree$tip.label)] <- tree$tip.label
  #colnames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])] <- paste0("n", colnames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])])
  
  inboth <- intersect(rownames(dist.mat), rownames(diff.mat))
  
  # Get unique combinations of this set:
  unique.combo <-  combn(inboth, m = 2)
  
  # make vectors to hold results
  dist_trait <- rep(NA, ncol(unique.combo)) # the trait distances
  dist_tree <-  rep(NA, ncol(unique.combo)) # the tree distances
  
  TTR <- data.frame(species1 = unique.combo[1,], 
                    species2 = unique.combo[2,], 
                    dist_trait, dist_tree, stringsAsFactors=F)
  
  # now fill it with the pairwise comparisons
  for (ii in 1:nrow(TTR)){
    TTR$dist_trait[ii]  <- diff.mat[TTR$species1[ii], TTR$species2[ii]]
    TTR$dist_tree[ii]   <- dist.mat[TTR$species1[ii], TTR$species2[ii]]    
  }
  #TTR$dist_tree <- TTR$dist_tree/2
  return(list(TTR = TTR, distance.matrix = dist.mat, difference.matrix = diff.mat))
}

ttrr <- build.TTR(tree, trait.data)


#### Pull out pairwise distances between ALL sister TIPS and NODES
sister.pairs <- function(phy, TTR.object) {
  all.matches <- NULL
  for(p in 1:length(unique(rownames(TTR.object$distance.matrix)))) {
    taxon <- unique(rownames(TTR.object$distance.matrix))[p] # choose a current taxon
    
    if (p>(Ntip(phy))) {
      sisters <- getSisters(phy, node=(p+1), mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<=Ntip(phy)) {
        pair.matrix[,2] <- rownames(TTR.object$distance.matrix)[sisters]
      } else {
        pair.matrix[,2] <- rownames(TTR.object$distance.matrix)[sisters-1] # add all the comparisons
      }
    } else {
      sisters <- getSisters(phy, node=p, mode="number") # get the sister node/tips to this taxon
      pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
      pair.matrix[,1] <- taxon # who they're being compared against (taxon)
      if (sisters<Ntip(phy)) {
        pair.matrix[,2] <- rownames(TTR.object$distance.matrix)[sisters]
      } else {
        pair.matrix[,2] <- rownames(TTR.object$distance.matrix)[sisters-1] # add all the comparisons
      }
    }
    all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
  }
  return(all.matches)
}

sissies <- sister.pairs(tree, ttrr)


#### Extract matching pairs from our Tree/Trait/Range data frame
get.matches <- function(TTR.object, Sister.object) {
  matches <- NULL
  all.matches <- Sister.object
  for (i in 1:length(all.matches[,1])) {
    current  <- filter(TTR.object$TTR, species1==all.matches[i,1] | species2==all.matches[i,1])
    current.matches <- subset(all.matches, 
                              all.matches[,1]==all.matches[i,1] | all.matches[,2]==all.matches[i,1])
    for (zz in 1:length(current.matches[,1])) {
      match1 <- subset(current, current[,1]==current.matches[zz,1] & current[,2]==current.matches[zz,2])
      match2 <- subset(current, current[,1]==current.matches[zz,2] & current[,2]==current.matches[zz,1])
      matches <- rbind(matches, match1)
      matches <- rbind(matches, match2)
    }
  }
  #all.sister.pairs <- rbind(all.sister.pairs, matches)
  matches <- unique(matches[,1:length(matches)])
  return(matches)
}

single.matches <- get.matches(ttrr, sissies)


#### Determine geographic overlap of all pairs of sisters
overlaps <- function(get.matches.object, combined.polygons) {
  current <- get.matches.object
  for (t in 1:length(current[,1])) {
    overlap <- gOverlaps(combined.polygons[[current[t,1]]], combined.polygons[[current[t,2]]])
    current[t,5] <- overlap
  }
  return(current)
}

comb.poly <- c(dist.data$ConvexHulls, dist.data$Tip_ConvexHulls) # combine the node and tip distributions
full.TTR <- overlaps(single.matches, comb.poly) # actually do the pairwise comparisons of sisters
colnames(full.TTR) <- c("species1", "species2", "dist_trait", "dist_tree", "range_overlap")

# save the output for plotting again later
saveRDS(full.TTR, "/Users/Ian/Desktop/MioceneAustralia/TTR_Table_Marsupials.RDS")

# if you want to plot the results from all groups together
a.TTR <- readRDS("/Users/Ian/Desktop/MioceneAustralia/TTR_Table_Agamids.RDS")
m.TTR <- readRDS("/Users/Ian/Desktop/MioceneAustralia/TTR_Table_Marsupials.RDS")
p.TTR <- readRDS( "/Users/Ian/Desktop/MioceneAustralia/TTR_Table_Pygopodoidea.RDS")
b.TTR <- readRDS( "/Users/Ian/Desktop/MioceneAustralia/TTR_Table_Meliphagoidea.RDS")

all.TTR <- rbind(a.TTR, rbind(m.TTR, rbind(p.TTR, b.TTR)))
all.TTR <- rbind(m.TTR, rbind(p.TTR, b.TTR))

mid.miocene <- filter(all.TTR, dist_tree <= 23)

lines <- (ggplot(mid.miocene, aes(x=dist_tree, y=dist_trait, color=range_overlap))
          #+ geom_point(alpha=0.1)
          + geom_smooth(aes(fill=range_overlap))
          + scale_x_reverse()
          + ggtitle("Marsupial Mammals:
                    Trait Distance between Overlapping and Non-overlapping Taxa")
          + theme_classic())
square <- (ggplot(mid.miocene)
           + geom_density(aes(dist_tree, ..count.., fill=range_overlap, alpha=0.5), position="fill") # or dist_trait to see if the two modes differ
           + scale_x_reverse()
           + scale_y_reverse(position="right")
           + theme_classic())
square + expand_limits(y=0.8,0.5)

multiplot(lines, square, layout=matrix(c(1,1,2), nrow=1, byrow=T))

test <- filter(mid.miocene, range_overlap=="FALSE")

ggplot(mid.miocene, aes(x=dist_tree, ..count..)
       + geom_point(aes(fill=range_overlap)))
(ggplot(mid.miocene)
  + geom_density(aes(dist_tree, ..count.., fill=range_overlap, alpha=0.5)))
