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

tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.tre")
  trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")
    trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.logMASS.csv", header=F)
#range <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Pygopodoidea.Range.Test.csv", header=T, row.names=1)
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Meliphagoids.csv", header=T)
  distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
    distribution <- distribution[complete.cases(distribution),] 
# as.data.frame(table(distribution$Name_in_Tree))

### We don't have distributional data for everything, so we need to drop a few from the tree
  drop <- setdiff(trees[[1]]$tip.label, unique(distribution$Name_in_Tree))
  tree <- drop.tip(tree, tip=drop);
    trees <- lapply(trees, drop.tip, tip=drop); class(trees) <- "multiPhylo"
### and from the trait data
  trait <- trait[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
    traits <- trait[,2]; names(traits) <- trait[,1] #read in data file in RPANDA format
  #sim.traits <- readRDS("/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/ASR_TRC.Skinks.Traits.RDS")
  sim.traits <- readRDS("/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoids.Traits.RDS")
    


all.TTR <- NULL   
for (pp in 1:length(trees)) {
  cat("iteration", pp, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  tree <- trees[[pp]]
  traits <- sim.traits[[pp]] # ONLY IF YOU'RE USING SIMULATED DATA!
  
  ## Build a distance matrix for all taxa in the tree
  #dist.mat <- cophenetic.phylo(tree) # only works the tips
  #dist.mat <- dist.nodes(tree) # includes internal and terminal nodes
  dist.mat <- max(nodeHeights(tree)) - vcvPhylo(tree)
       # for (uu in 1:nrow(dist.mat)) {
       #   for (jj in 1:ncol(dist.mat)) {
       #     if (uu > Ntip(tree) && jj > Ntip(tree)) {
       #       dist.mat[[uu,jj]] <- 0
       #     }
       #   }
       # }
       # 
       # node_distances <- mrca(tree, full=T)
       # for (uu in 1:nrow(node_distances)) {
       #   for (jj in 1:ncol(node_distances)){
       #     if (uu < Ntip(tree)+1 && jj < Ntip(tree)+1){
       #       node_distances[[uu,jj]] <- 0
       #     } else {
       #       node_distances[[uu,jj]] <- max(nodeHeights(tree)) - nodeheight(tree, node_distances[[uu,jj]])
       #     }
       #   }
       # }
       # dist.mat <- node_distances + dist.mat
  
  
  ## Build a distance matrix for the trait of interest for all taxa in the tree
      # IF YOU'RE USING EMPIRICAL DATA (use this, not the below)
      #ANC <- fastAnc(tree, traits, var=F, CI=F)
      #all.traits <- append(traits, ANC)
      
      # IF YOU'RE USING SIMULATED DATA (use this, not the above)
      all.traits <- traits
  
  diff.mat <- matrix(NA,length(all.traits), length(all.traits))
  colnames(diff.mat) <- names(all.traits); rownames(diff.mat) <- names(all.traits)
  diff.mat <- as.data.frame(diff.mat)
  for (j in 1:length(rownames(diff.mat))) {
    for (k in 1:length(colnames(diff.mat))) {
      diff.mat[j,k] <- abs(all.traits[j]-all.traits[k])
    }
  }
  ## IF USING SIMULATED DATA COMMENT OUT THE RENAMING BELOW!
  #rownames(diff.mat)[(length(tree$tip.label)+1):length(diff.mat[,1])] <- paste("Node", rownames(diff.mat)[(length(tree$tip.label)+1):length(diff.mat[,1])], sep = ".")
  #colnames(diff.mat)[(length(tree$tip.label)+1):length(diff.mat[,1])] <- paste("Node", colnames(diff.mat)[(length(tree$tip.label)+1):length(diff.mat[,1])], sep = ".")
  
  
  ## Now build a data frame to hold the Tree/Trait/Range distance data
  # First, determine tips in common:
  #rownames(dist.mat)[1:length(tree$tip.label)] <- tree$tip.label
    rownames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])] <- paste("Node", rownames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])], sep = ".")
  #colnames(dist.mat)[1:length(tree$tip.label)] <- tree$tip.label
    colnames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])] <- paste("Node", colnames(dist.mat)[(length(tree$tip.label)+1):length(dist.mat[,1])], sep = ".")
  
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
  all.TTR[[pp]] <- TTR
}
#saveRDS(all.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoids.Matrices.RDS")

    ### If you're using simulated data, use the steps below here
    #test <- readRDS("/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Agamids.Matrices.RDS")
    #for (kk in 1:length(all.TTR)) {
    #  all.TTR[[kk]]$range_overlap <- test[[kk]]$range_overlap
    #}

#### If you've already calculated the necessary information
#for (rr in 1:length(all.TTR)) {
#  cat("iteration", rr, "of", length(all.TTR), "\n") #keep track of what tree/loop# we're on
#  TTR <- all.TTR[[rr]]
#  for(yy in 1:length(TTR[,1])){
#    if (!(TTR[yy,1] %in% tree$tip.label) | !(TTR[yy,2] %in% tree$tip.label)) {
#      taxon1 <- which(rownames(dist.mat)==TTR[yy,1])
#      if (taxon1 > length(tree$tip.label)){
#        taxon1 <- min(getDescendants(tree, taxon1))
#      }
#      taxon2 <- which(rownames(dist.mat)==TTR[yy,2])
#      if (taxon2 > length(tree$tip.label)){
#        taxon2 <- min(getDescendants(tree, taxon2))
#      }
#      TTR[yy,4] <- max(nodeHeights(tree))-(nodeheight(tree, node=getMRCA(tree,c(taxon1,taxon2))))
#    }
#  }
#  all.TTR[[rr]] <- TTR
#}


## Next we need to determine (pairwise) if taxa overlap in their ranges
  ## this step only needs to be done once! (won't change with changes to the tree)
    ## this version works by assuming an ancestral node has the distribution of its daughters, combined.
all.taxa <- unique(distribution$Name_in_Tree)
all.sp <- list()
all.hull <- list()
#all.poly <- list()
sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")

for (p in 1:length(unique(distribution$Name_in_Tree))) {
  cat("iteration", p, "of", length(unique(distribution$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
  current.taxon <- all.taxa[[p]]
  current.data <- subset(distribution, distribution$Name_in_Tree == current.taxon)
  if (nrow(current.data) > 1000) {
    current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
  } 
  points <- current.data[,c(2,3)]
  all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
  all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)

  # only plot the maps below if you need to check/clean the data
    #fortified.data <- fortify(distribution.hull)
    #pdf(paste("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Skinks_Maps/", current.taxon, ".pdf", sep=""))
    #print(ggmap(sq_map) 
    #  + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
    #  + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
    #dev.off()
}
names(all.sp) <- unique(distribution$Name_in_Tree)
names(all.hull) <- unique(distribution$Name_in_Tree)

# make a loop that creates distributions at nodes by combining daughter distributions
#for (kk in length(rownames(dist.mat)):(length(tree$tip.label)+1)){
#  descendantz <- Descendants(tree, kk, type="children")
#  combined.range <- gUnion(all.hull[[descendantz[1]]], all.hull[[descendantz[2]]])
#  all.hull[[kk]] <- combined.range
#}
polygon.list <- NULL
for (ii in 1:length(trees)) {
  cat("iteration", ii, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  for (kk in (length(rownames(dist.mat))+1):(length(tree$tip.label)+2)){
    descendantz <- Descendants(trees[[ii]], kk, type="children")
    combined.range <- gUnion(all.hull[[descendantz[1]]], all.hull[[descendantz[2]]])
    all.hull[[kk]] <- combined.range
  }
  polygon.list[[ii]] <- all.hull
}
for (zz in 1:length(polygon.list)){
  #names(polygon.list[[zz]])[length(polygon.list[[zz]]: (length(trees[[1]]$tip.label)+1)] <- rownames(dist.mat)[(length(trees[[1]]$tip.label)+1):length(rownames(dist.mat))]
  names(polygon.list[[zz]])[(length(polygon.list[[zz]])-1) : (Ntip(trees[[1]])+1)] <- rownames(dist.mat)[(Ntip(trees[[1]])+1):(length(rownames(dist.mat)))]
  names(polygon.list[[zz]])[length(polygon.list[[zz]])] <- paste("Node", (Ntip(trees[[1]])+1), sep=".")
  #names(polygon.list[[zz]])[(length(Ntip(trees[[1]]))):(length(polygon.list[[zz]]))] <- rownames(dist.mat)[length(rownames(dist.mat))]:rownames(dist.mat)[Ntip(phy)+1]
  #current.poly <- polygon.list[[zz]]
  #for (pp in length(current.poly):(Ntip(phy)+1)) {
  #  names(current.poly[[pp]]) <- rownames(dist.mat)[Ntip(phy)+pp]
  #}
}

#saveRDS(all.sp,   file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.SpatialPoints.RDS")
#saveRDS(all.hull, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.Polygons.RDS")
saveRDS(polygon.list, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.AllPolygons.RDS")


#range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
#arrange <- data.frame(species1 = unique.combo[1,], 
#                      species2 = unique.combo[2,], 
#                      range_overlap, stringsAsFactors=F)
#for (t in 1:length(arrange[,1])) {
#  overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
#  arrange[t,3] <- overlap
#}
# Get unique combinations of this set:
inboth <- intersect(rownames(dist.mat), rownames(diff.mat)) 
unique.combo <-  combn(inboth, m = 2) # check ncol(unique.combo) matches length(all.TTR[[1]][,1])

### This will determine pairwise overlap for ALL tips and nodes
  ### because of this, it takes FOREVER, so instead, I've implemented this step later
#arrange.list <- NULL
#for (ii in 1:length(trees)) {
#  cat("iteration", ii, "of", length(trees), "\n") #keep track of what tree/loop# we're on
#  range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
#  arrange <- data.frame(species1 = unique.combo[1,], 
#                        species2 = unique.combo[2,], 
#                        range_overlap, stringsAsFactors=F)
#  for (t in 1:length(arrange[,1])) {
#    overlap <- gOverlaps(polygon.list[[ii]][[arrange[t,1]]], polygon.list[[ii]][[arrange[t,2]]])
#    arrange[t,3] <- overlap
#  }
#  arrange.list[[ii]] <- arrange
#}


#all.TTR <- lapply(all.TTR, function(x) {x[,5] <- NULL; x}) # if you need to drop a column from the 'all.TTR' list
#for (jj in 1:length(all.TTR)){
#  all.TTR[[jj]] <- cbind(all.TTR[[jj]], arrange.list[[jj]]$range_overlap)
#  colnames(all.TTR[[jj]]) <- c("species1", "species2", "dist_trait", "dist_tree", "range_overlap")
#}
#all.TTR <- lapply(all.TTR, cbind, arrange$range_overlap)
#all.TTR <- lapply(all.TTR, setNames, c("species1", "species2", "dist_trait", "dist_tree", "range_overlap"))
saveRDS(all.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Skinks.ASR_TRC.Matrices.RDS")


################################################################################
## If you already have an 'all.TTR' object with range_overlap of all tips+nodes
### and all you need is to recalculate the dist_tree, then start HERE
################################################################################

## Pull out all sister species for comparison
## in a situation like Hesperoedura, this includes comparison to just the sister node:
#### This version gets pairwise distances between only sister TIPS (and a tip+node if necessary)
#all.out <- NULL
#for (y in 1:length(trees)) {
#  all.matches <- NULL
#  phy <- trees[[y]]
#  for(p in 1:length(phy$tip.label)) {
#    taxon <- phy$tip.label[p] # choose a current taxon
#    sisters <- getSisters(phy, node=taxon, mode="number") # get the sister node/tips to this taxon
#    
#    #all.desc <- getDescendants(phy, sisters) # get all the descendants of that sister node
#    #tip.desc <- subset(all.desc, all.desc <= length(phy$tip.label)) # keep only the tips (terminal nodes)
#    
#    pair.matrix <- matrix(NA, ncol=2, nrow=1) # make an empty matrix for the pairwise comparisons
#    pair.matrix[,1] <- taxon # who they're being compared against (taxon)
#    pair.matrix[,2] <- rownames(dist.mat)[sisters] # add all the comparisons
#    
#    #for (q in 1:length(tip.desc)) {
#    #  pair.matrix[q,2] <- phy$tip.label[tip.desc[q]] # add all the comparisons
#    #  #all.matches <- rbind(all.matches, pair.matrix)
#    #}
#    all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
#  }
#  all.out[[y]] <- all.matches
#}
#all.matches <- unique(all.matches[,1:2]) # drop any duplicates, which seem to happen
#all.matches <- lapply(all.out, unique)
## Pull out all sister species for comparison
## in a situation like Hesperoedura, this includes comparison to just the sister node:
#all.matches <- NULL # empty object

#### This alternative version gets pairwise distances between ALL sister TIPS and NODES
all.out <- NULL
for (y in 1:length(trees)) {
  all.matches <- NULL
  phy <- trees[[y]]
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
  all.out[[y]] <- all.matches
}


## Now we need to get the matching data from our Tree/Trait/Range data frame
#all.sister.pairs <- NULL
all.sister.bytree <- NULL
for (tt in 1:length(all.TTR)) {
  matches <- NULL
  cat("iteration", tt, "of", length(all.TTR), "\n") #keep track of what tree/loop# we're on
  #sister.pairs <- NULL
  all.matches <- all.out[[tt]]
  #for (i in 1:length(unique(all.TTR[[tt]]$species1))) {
  for (i in 1:length(all.matches[,1])) {
      current  <- filter(all.TTR[[tt]], species1==all.matches[i,1] | species2==all.matches[i,1])
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
  all.sister.bytree[[tt]] <- matches
  #all.sister.pairs <- rbind(all.sister.pairs, matches)
}
#all.sister.bytree <- lapply(all.sister.bytree, function(x) {x[,5] <- NULL; x}) # if you need to drop a column from the 'all.TTR' list

### Determine geographic overlap of all pairs of sisters
for (ii in 1:length(trees)) {
  cat("iteration", ii, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  current <- all.sister.bytree[[ii]]
  #range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
  #arrange <- data.frame(species1 = unique.combo[1,], 
  #                      species2 = unique.combo[2,], 
  #                      range_overlap, stringsAsFactors=F)
  for (t in 1:length(current[,1])) {
    overlap <- gOverlaps(polygon.list[[ii]][[current[t,1]]], polygon.list[[ii]][[current[t,2]]])
    current[t,5] <- overlap
  }
  all.sister.bytree[[ii]] <- current
}
all.sister.bytree <- lapply(all.sister.bytree, setNames, c("species1", "species2", "dist_trait", "dist_tree", "range_overlap"))

#singles <- unique(all.sister.pairs[,1:5]) # drop any duplicates, which seem to happen
bytree.final <- NULL
all.sister.pairs <- NULL
for (qq in 1:length(all.sister.bytree)) {
  matches <- NULL
  cat("iteration", qq, "of", length(all.sister.bytree), "\n") #keep track of what tree/loop# we're on
  #sister.pairs <- NULL
  final.frame <- all.sister.bytree[[qq]]
  #tree <- trees[[qq]]
  #node_distances <- mrca(tree, full=T)
  
 # for (yy in 1:length(final.frame[,1])) {
 #   if (!(final.frame[yy,1] %in% tree$tip.label) | !(final.frame[yy,2] %in% tree$tip.label)) {
 #     taxon1 <- which(rownames(dist.mat)==final.frame[yy,1])
 #     taxon2 <- which(rownames(dist.mat)==final.frame[yy,2])
#
 #     final.frame[yy,4] <- max(nodeHeights(tree))-(nodeheight(tree, node_distances[[taxon1, taxon2]]))
 #   }
    # node_distances <- mrca(tree, full=T)
    # for (uu in 1:nrow(node_distances)) {
    #   for (jj in 1:ncol(node_distances)){
    #     if (uu < Ntip(tree)+1 && jj < Ntip(tree)+1){
    #       node_distances[[uu,jj]] <- 0
    #     } else {
    #       node_distances[[uu,jj]] <- max(nodeHeights(tree)) - nodeheight(tree, node_distances[[uu,jj]])
    #     }
    #   }
    # }
  #}
  bytree.final[[qq]] <- final.frame
  all.sister.pairs <- rbind(all.sister.pairs, final.frame)
}

#saveRDS(all.sister.bytree, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.Comparisons.ByTree.RDS")
#saveRDS(all.sister.pairs,  file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.All.Unique.Comparisons.RDS")
saveRDS(bytree.final,      file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoids.AllNodes.Comparisons.ByTree.RDS")
saveRDS(all.sister.pairs,  file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoids.AllNodes.All.Unique.Comparisons.RDS")

#split <- filter(short, range_overlap==F) # if you want to just plot one of the groups
#### Plot the trends through time
  lines <- (ggplot(all.sister.pairs, aes(x=dist_tree, y=dist_trait, color=range_overlap))
    #+ geom_point(alpha=0.1)
    + geom_smooth(aes(fill=range_overlap))
    + scale_x_reverse()
    + ggtitle("Marsupial Mammals:
  Trait Distance between Overlapping and Non-overlapping Taxa")
    + theme_classic())

  #ggMarginal(t, type="density", fill="transparent")
#### Plot the 
  square <- (ggplot(all.sister.pairs)
    + geom_density(aes(dist_tree, ..count.., fill=range_overlap, alpha=0.5), position="fill") # or dist_trait to see if the two modes differ
    + scale_x_reverse()
    + scale_y_reverse(position="right")
    + theme_classic())
    #+ geom_vline(data=means, aes(xintercept=grp.mean, color=range_overlap),
    #              linetype="dashed"))
  multiplot(lines, square, layout=matrix(c(1,1,2), nrow=1, byrow=T))
### Plot the frequency of each type of event through time
###(ggplot(all.sister.pairs, aes(dist_tree))
###  + geom_area(aes(y= ..count.., fill=range_overlap, 
###                  group=range_overlap, alpha=0.5), stat="bin", binwidth=1) # or dist_trait to see if the two modes differ
###  + scale_x_reverse()
###  + theme_classic())


short.d <- filter(bytree.final[[3]], dist_tree < 23)
#short.d <- filter(short.d, dist_trait < 0.4)
lines <- (ggplot(short.d, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  #+ geom_point(alpha=0.1)
  + geom_smooth(aes(fill=range_overlap))
  #+ geom_smooth()
  + scale_x_reverse()
  + ggtitle("Agamid Dragons - BM Simulated Data:
Trait Distance between Overlapping and Non-overlapping Taxa")
  + theme_classic())
#### Plot the 
square <- (ggplot(short.d)
  + geom_density(aes(dist_tree, ..count.., fill=range_overlap, alpha=0.5), position="fill") # or dist_trait to see if the two modes differ
  + scale_x_reverse()
  + scale_y_reverse()
  + theme_classic()
  + ggtitle())
#+ geom_vline(data=means, aes(xintercept=grp.mean, color=range_overlap),
#              linetype="dashed"))
multiplot(lines, square, layout=matrix(c(1,1,2), nrow=1, byrow=T))
#multiplot(lines, square, cols=2)


pygo <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Pygopodoids.AllNodes.All.Unique.Comparisons.RDS")
  pygo$clade <- "Pygopodoidea"
mars <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Marsupials.AllNodes.All.Unique.Comparisons.RDS")
  mars$clade <- "Marsupials"
agam <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/New.Agamids.AllNodes.All.Unique.Comparisons.RDS")
  agam$clade <- "Agamidae"
bird <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/NEW.Meliphagoids.AllNodes.All.Unique.Comparisons.RDS")
  bird$clade <- "Meliphagoidea"
skink <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/New.Skinks.AllNodes.All.Unique.Comparisons.RDS")
  skink$clade <- "Skinks"
all.clades <- rbind(pygo, mars, agam, bird, skink)

#### Above line FINISHED ################################ Below line NOT FINISHED

all.clades <- rbind(agam, pygo, bird, skink, mars)
  saveRDS(all.clades, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/AllClades.AllTree.Comparisons.RDS")
    all.1 <- rbind(agam.1, bird.1)
      all.1 <- rbind(all.1, mars.1)
        all.1 <- rbind(all.1, skink.1)
          all.1 <- rbind(all.1, pygo.1)
multiplot(agam.plot, pygo.plot, bird.plot, mars.plot, skink.plot, ncol=1)
multiplot(agam.plot, pygo.plot, bird.plot, mars.plot, skink.plot,
          layout=matrix(c(1,2,3,4,5), nrow=2, byrow=F))


short <- filter(all.clades, dist_tree <= 23)
#split <- filter(short, range_overlap==F) # if you want to just plot one of the groups
split <- filter(short, dist_trait<0.4)
(ggplot(short , aes(x=dist_tree, y=dist_trait, color=range_overlap)) # alternatively color=clade
              + geom_point(alpha=0.1) # aes(shape=clade)
              + geom_smooth(aes(fill=range_overlap, 
                                linetype=range_overlap, color=range_overlap)) #alternatively fill=clade or fill=range_overlap
              + scale_x_reverse()
              + ggtitle("All Clades:
Trait Distance between Overlapping and Non-overlapping Taxa")
              + theme_classic()) # remember you can set the x axis limits with + xlim(40,0)



short.dff <- subset(dff, dff$dist_tree <= 23)

(ggplot(all.sister.pairs, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.25)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse()
  + ggtitle("Pygopodoids Trait Distance between Overlapping and Non-overlapping Taxa"))
(ggplot(fuller, aes(x=dist_trait, fill=range_overlap))
  + geom_density(alpha=0.5))


# if you want to plot a series of species within a given genus:
test <- distribution[which(distribution$Name_in_Tree %in% all.taxa[175:178]),]
sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
(ggmap(sq_map) + geom_point(data = test, mapping = aes(x = Longitude, y = Latitude, color=Name_in_Tree)))


pygo <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.Matrices.RDS")

test <- dff
test <- rbind(test, agam)



## Set your parameters relevant to your empirical parameter estimates
######################################################
# alpha
#alpha = 2 # either a static value
diff.alpha <- seq(from=1.5, to=4.5, by=0.01) # or set it as a vector of sampled values
diff.alpha <- sample(diff.alpha, size=100, replace=T) # or set it as a vector of sampled values
# pre/post shift sigma
preshift.sigma  = 0.01
postshift.sigma = 0.2
# and shift time
sim.shifts <- seq(from=5, to=9, by=1)
sim.shifts <- sample(sim.shifts, size=100, replace=T)
######################################################
bm.pars <- 0.001 # set the diffusion parameter of the BM process
ou.pars <- c(0.01, sample(diff.alpha, 1), 1) # set the diffusion parameter, the alpha, and the optimum

##########################################################################################################
## LOOP 1:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either the SRC ('Single-rate-constraint', BM-OU, 1 sig2) or the 
### TRC ('Two-rate-constraint', BM-OU, 2 sig2) model, then comparatively fit a set of standard models 
### (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a different model use the second loop, to simulate under a Brownian
### motion or Ornstein-Uhlenbeck process, using Diversitree (although you could also use Geiger).

### This loop will simulate data onto the tree with fossil tips, fit a series of models, then
### drop the extinct tips and associated data, refit the same models, and provide a summary for each
##########################################################################################################
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
    
  # Option C (comment out top when simulating different times, comment out bottom when simulating different alphas)
  m[[2]] <- ouMatrix(m[[2]], alpha=diff.alpha[z]) # transform the second half of the vcv matrix according to your alpha
  #m[[2]] <- ouMatrix(m[[2]], alpha=alpha) # transform the second half of the vcv matrix according to your alpha
  
  # Option B (adjust to change simulating model, comment out both to get the SRC model)
  m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
  m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
    
    
  m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
    
  # OR, do it like the 'ecological release' model
  # m <- lapply(m, function(x) x*0.1)
  # m.rev.rel <- m[[1]] + m[[2]]
    
  # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
  traitz <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), 
                                    sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))
  sim.traits[[z]] <- traitz  
  saveRDS(sim.traits, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SIMULATED.TRC.Meliphagoid.Traits.RDS")
}



##########################################################################################################
## LOOP 2:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either Brownian Motion (BM) or Ornstein-Uhlenbeck (OU) models then comparatively fit a set 
### of standard models (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a mode-variable model, use the first loop.
##########################################################################################################
bm.pars <- 0.004 # set the diffusion parameter of the BM process
ou.pars <- c(0.5, sample(diff.alpha, 1), 1) # set the diffusion parameter, the alpha, and the optimum
######################################################
## Set empty object to hold your results
sim.traits.geiger <- list(); sim.traits.ouwie <- list() # make trait lists for all trees in geiger and ouwie data format
save.sim.traits <- NULL
sim.traits <- NULL
extant.data <- list()
# make sure to the outputs depending on if you're simulating different times, or alphas
# this means changing the 'm' and the 'm[[2]]' objects below!
num.sims <- 1 # designate the number of simulated data sets you like to create per tree
### if trees are sensitive to shift dates: (this is for the Plio-Pleistocene trees)
######################################################

for (z in 1:1) {
  traits.geiger <- NULL; traits.ouwie <- NULL # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  
  simulated.traits <- NULL
  simulated.traits <- as.data.frame(sim.character(phy, model="bm", bm.pars, x0=2))
  dataz <- simulated.traits[,1]
  names(dataz) <- rownames(simulated.traits)
  out.data <- dataz[1:length(phy$tip.label)]
  
  sim.traits[[z]] <- out.data
  saveRDS(sim.traits, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Simulated.BM.Meliphagoid.Traits.RDS")
}





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
test <- split.vcv(trees[[1]], 10)
nrow(test$mat1)
test <- split.vcv.internal(trees[[1]], 10)
nrow(test$mat1)
#













sim.TTR <- NULL   
for (pp in 1:length(trees)) {
  cat("iteration", pp, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  tree <- trees[[pp]]
  
  ## Build a distance matrix for all taxa in the tree
  dist.mat <- cophenetic.phylo(tree)
  
  ## Build a distance matrix for the trait of interest for all taxa in the tree
  diff.mat <- matrix(NA,length(sim.traits[[pp]]), length(sim.traits[[pp]]))
  colnames(diff.mat) <- names(sim.traits[[pp]]); rownames(diff.mat) <- names(sim.traits[[pp]])
  diff.mat <- as.data.frame(diff.mat)
  for (j in 1:length(tree$tip.label)) {
    for (k in 1:length(tree$tip.label)) {
      diff.mat[j,k] <- abs(sim.traits[[pp]][j]-sim.traits[[pp]][k])
    }
  }
  
  ## Now build a data frame to hold the Tree/Trait/Range distance data
  # First, determine tips in common:
  inboth <- intersect(tree$tip.label, names(sim.traits[[pp]]))
  
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
  TTR$dist_tree <- TTR$dist_tree/2
  ## Save each/all of the Distance/Data Matrices externally
  sim.TTR[[pp]] <- TTR
}
saveRDS(sim.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SRC.Simulated.Agamids.Matrices.RDS")


#all.TTR <- lapply(all.TTR, function(x) {x["range_overlap"] <- NULL; x}) # if you need to drop a column from the 'all.TTR' list
sim.TTR <- lapply(sim.TTR, cbind, arrange$range_overlap)
sim.TTR <- lapply(sim.TTR, setNames, c("species1", "species2", "dist_trait", "dist_tree", "range_overlap"))
saveRDS(sim.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SRC.Simulated.Agamids.Matrices.RDS")


## Pull out all sister species for comparison
## in a situation like Hesperoedura, this includes comparison to all:
## Oedura, Amalosia, Nebulifera (make sense?)
all.matches <- NULL # empty object
for(p in 1:length(tree$tip.label)) {
  taxon <- tree$tip.label[p] # choose a current taxon
  sisters <- getSisters(tree, node=taxon, mode="number") # get the sister node/tips to this taxon
  all.desc <- getDescendants(tree, sisters) # get all the descendants of that sister node
  tip.desc <- subset(all.desc, all.desc <= length(tree$tip.label)) # keep only the tips (terminal nodes)
  
  pair.matrix <- matrix(NA, ncol=2, nrow=(length(tip.desc))) # make an empty matrix for the pairwise comparisons
  pair.matrix[,1] <- taxon # who they're being compared against (taxon)
  for (q in 1:length(tip.desc)) {
    pair.matrix[q,2] <- tree$tip.label[tip.desc[q]] # add all the comparisons
  }
  all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
}
#all.matches <- unique(all.matches[,1:2]) # drop any duplicates, which seem to happen

## Now we need to get the matching data from our Tree/Trait/Range data frame
all.sister.pairs <- NULL
for (tt in 1:length(sim.TTR)) {
  cat("iteration", tt, "of", length(sim.TTR), "\n") #keep track of what tree/loop# we're on
  sister.pairs <- NULL
  for (i in 1:length(unique(sim.TTR[[tt]]$species1))) {
    current <- filter(sim.TTR[[tt]], species1==all.matches[i,1])
    matches <- subset(current, current$species2==all.matches[i,2])
    sister.pairs <- rbind(sister.pairs, matches)
  }
  all.sister.pairs <- rbind(all.sister.pairs, sister.pairs)
}
saveRDS(all.sister.pairs, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/SRC.Simulated.Agamids.AllTree.Comparisons.RDS")

(ggplot(all.sister.pairs, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.1)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse()
  + ggtitle("Meliphagoid Birds:
            Trait Distance between Overlapping and Non-overlapping Taxa")
  + theme_classic())

short.d <- filter(all.sister.pairs, dist_tree <= 23)
sim.agam.plot <- (ggplot(short.d, aes(x=dist_tree, y=dist_trait, color=range_overlap))
              #+ geom_point(alpha=0.1)
              + geom_smooth(aes(fill=range_overlap))
              + scale_x_reverse()
              + ggtitle("Agamid Dragons Simulated Under SRC:
                        Trait Distance between Overlapping and Non-overlapping Taxa")
              + theme_classic())



densityplot(agam$dist_tree)

sub23 <- filter(all.clades, dist_tree <= 23)

(ggplot(sub23)
  + geom_density(aes(dist_tree, fill=range_overlap, alpha=0.5))
  + theme_classic())





overs <- filter(mars.1, dist_tree < 15 & range_overlap==T & dist_trait > 0.1)


#


































##### Below this line is the previous attempt to determine overlap using biomes, DISREGARD
##########################################################################################
range_overlap <- rep(NA, ncol(ucomb))
arrange <- data.frame(species1 = ucomb[1,], species2 = ucomb[2,], range_overlap, stringsAsFactors=F)
for (t in 1:length(arrange[,1])) {
  taxon1 <- subset(range, rownames(range)==arrange[t,1])
  t1 <- c(as.character(taxon1[[1]]), 
          as.character(taxon1[[2]]), 
          as.character(taxon1[[3]]), 
          as.character(taxon1[[4]]),
          as.character(taxon1[[5]]))
  taxon2 <- subset(range, rownames(range)==arrange[t,2])
  t2 <- c(as.character(taxon2[[1]]), 
          as.character(taxon2[[2]]), 
          as.character(taxon2[[3]]), 
          as.character(taxon2[[4]]),
          as.character(taxon2[[5]]))
  
  overlap <- intersect(t1, t2)
  overlap <- overlap[!is.na(overlap)]
  if (length(overlap > 0)) {
    arrange[t,3] <- "yes"
  } else {
    arrange[t,3] <- "no"
  }
}

dff <- cbind(dff, arrange$range_overlap)
colnames(dff) <- c("species1", "species2", "dist_trait", "dist_tree", "range_overlap")
short.dff <- subset(dff, dff$dist_tree <= 20)
saveRDS(dff, file="/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Pygopodoidea.Tree.Trait.Range.Matrices.RDS")


(ggplot(dff, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point()
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse())

(ggplot(short.dff, aes(x=dist_trait, fill=range_overlap))
  + geom_density(alpha=0.5))


saveRDS(all.taxa)
all.taxa <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Agamids.Tree.Trait.Range.Matrices.RDS")
all.taxa <- rbind.data.frame(all.taxa, dff)
shorty <- subset(all.taxa, all.taxa$dist_tree <= 20)
(ggplot(all.taxa, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.5)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse())


distributions <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Agamid_PointData/Agamids_PointData.csv", header=T)
distribution <- read.csv("/Users/Ian/Downloads/Clean_Agamids/Clean_Agamids.csv", header=T)
test <- unique(distribution$Species...matched)
tests <- unique(distributions$Genus_species)
setdiff(tests, test)
