
# function 'bind.at.depth' to choose where a new tip goes into the phylogeny
  # might be useful for 

# extinction is an additionally interesting derivative of missing sampling. 
# in studies that aren't looking at temporal patterns in trait change, you could
# simply add tips to an existing phylogeny ('bind.tip' or 'bind.at.depth').
# however, an extinct taxa is sampled in a given time period, then may not 
# be sampled in the subsequent period, changing the dynamics. 

library(phytools)
library(phangorn)
library(geiger)
library(mvtnorm)
library(MuMIn)
library(BioGeoBEARS)


#######################################################
# Read in all the data you'll use
######################################################
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")
  # or your empirical trees:
    agam <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.tre") # 100 tips
    mars <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.tre") # 133 tips
    bird <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.tre") # 149 tips
    pygo <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Pygopodoidea.tre") # 189 tips
    skink<- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Sphenomorphines.tre") # 240 tips

data  <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.logMASS.csv", 
                  row.names = 1, header=F) #read in data file in GEIGER format

#######################################################
#### Bind at Depth' function adds tips at desired age
#######################################################
# if the input 'length = NULL', then tip/edge will be added to keep tree ultrametric
bind.at.depth<-function(tree, tip.label, depth){
  H <- nodeHeights(tree)
  age.from.root<-max(H)-depth # how far is the shift from the root
  ii <- intersect(which(H[,1]<age.from.root),which(H[,2]>age.from.root)) # which edges span the desired placement depth
  edges <- tree$edge[ii,2] # what are the edge numbers, based on the above suitable matches?
  depths <- H[ii,2]-age.from.root
  jj <- sample(1:length(edges),1) # choose one of the suitable edges to paste the new tip to
  #mm <- sample(0.1:sample(1:depths[[jj]], 1), 1)
  bind.tip(tree, tip.label, where=edges[jj], position=depths[jj], 
           edge.length=(depths[jj]/(sample(1:depths[jj],1))))
  # by making the new branch length = depth/1:depth, this allows the
  # new branch to reach the present (depth/1) or be quite short (depth/depth)
}




#########################################################################
#### 'sim.fossil.trees' function adds tips to a base tree at desired age
#########################################################################
# phy = you input empirical phylogeny
# time.frame = vector that matches 'seq': c(young/lower.bound, old/upper.bound, interval)
# num.taxa = how many tips to add to each tree
# num.trees = number of trees to simulate
old.fossil.trees <- function (phy, time.frame, num.taxa, num.trees) {
  species<-as.character(c(1:num.taxa)) #identify how many new tips to include, package it up as "species"
  output.trees <<- NULL
  for (m in 1:num.trees){
    time.seq <- seq(time.frame[1], time.frame[2], time.frame[3]) # pick a time range
    # now a loop to add the desired number of taxa to the tree within the time period designated
    fossil.tree <- phy
    for(i in 1:length(species)) {
      cat(" tree", m ,"tip", i, "...")
      fossil.tree <- bind.at.depth(fossil.tree, species[i], sample(time.seq, 1, replace=T))
      output.tree <<- fossil.tree
    }      
    output.trees[[m]] <<- output.tree
  }
  class(output.trees) <<- "multiPhylo"
}

sim.fossil.trees <- function (phy, time.frame, num.taxa, num.trees) {
  species<-as.character(c(1:num.taxa)) #identify how many new tips to include, package it up as "species"
  output.trees <<- NULL
  iterations = 1
  while (iterations <= num.trees) {
    tryCatch({
      time.seq <- seq(time.frame[1], time.frame[2], time.frame[3]) # pick a time range
      # now a loop to add the desired number of taxa to the tree within the time period designated
      fossil.tree <- phy
      taxa.its <- 1
      while (taxa.its <= num.taxa) { 
        tryCatch({
          #cat(" tree", taxa.its ,"tip", i, "...")
          fossil.tree <- bind.at.depth(fossil.tree, species[taxa.its], sample(time.seq, 1, replace=T))
          output.tree <<- fossil.tree
          taxa.its = taxa.its+1
          }, error=function(e){})
        }
    }, error=function(e){})
    output.trees[[iterations]] <<- output.tree
    iterations = iterations+1
  }
  class(output.trees) <<- "multiPhylo"
}

# test out the function to see if it works:
sim.fossil.trees(phy=agam, time.frame=c(0.1,max(nodeHeights(agam)),0.1), num.taxa=100, num.trees=2) # output is called 'output.trees'
plot(output.trees)

####################################################################################
#### 'sim.fossil.data' function simulates trait data for the new tips in your trees
####################################################################################
# phy = your input phylogeny, or list of phylogenies as a multiPhylo object from 'sim.fossil.trees' (output.trees)
# num.taxa = the number of taxa you added to each phylogeny (same as in the fxn above)
# trait.data = trait data for the base tree (previous to adding the new tips), in geiger format (rownames=tips, col1=trait)
sim.fossil.data <- function (phy, num.taxa, trait.data) {
  species.list <- as.character(c(1:num.taxa)) #identify how many new tips to include, package it up as "species"
  fossil.data <<- NULL # empty object to hold all the simulated data (will be a list of length = length(phy))
  for (z in 1:length(phy)) {
    input.data <- NULL # empty object to hold data simulated for the current tree (data frame with length = ntips in tree + num.taxa)
    input.tree <- phy[[z]] # current phylogeny for simulating data
    new.species.data <- NULL # empty object to hold raw simulated data before moving it into 'input.data'
    #node.numbers <- as.data.frame(phy$tip.label); node.numbers[,"node.no"] <- rownames(node.numbers)
    edge.node.info <- cbind.data.frame(input.tree$edge, nodeHeights(input.tree))
    colnames(edge.node.info) <- c("descendant", "node", "start", "finish")
    
    for (t in 1:num.taxa) {
        sibs.list <- NULL # empty object to hold the list of siblings for our current tip
        sis.node <- getSisters(input.tree, paste(species.list[t]), mode="label") # find what the sister taxon or node is to your new tip
        outer.data <<- NULL
        if (length(sis.node$tips[1])==1 && (sis.node$tips[1] %in% species.list)) {
          parent.node <- findMRCA(input.tree, c(paste(species.list[t]), paste(sis.node$tips[1])), type="node")
          outer.data <<- parent.node
          sis.node <- Descendants(input.tree, parent.node, type="tips")
          for (i in 1:length(sibs[[1]])) {
            tip.num <- sibs[[1]][i]
            sibling <- input.tree$tip.label[tip.num]
            sibs.list <- append(sibs.list, sibling)
          } 
        } else if (is.null(sis.node$tips)) { 
          outer.data <<- sis.node$ndoes[1]
          sibs <- Descendants(input.tree, sis.node$nodes[1], type="tips") # from phangorn, get only the tips which descend from the node
          for (i in 1:length(sibs[[1]])) {
            tip.num <- sibs[[1]][i]
            sibling <- input.tree$tip.label[tip.num]
            sibs.list <- append(sibs.list, sibling)
          } 
        } else {
          sibs.list <- sis.node$tips[1]
          outer.data <<- sis.node$tips[1]
        }
        # now we can simulate data based on our real body sizes, start with a loop
        ## that pulls out empirical size data for all the sibling taxa
        size.dist <- NULL
        for (q in 1:length(sibs.list)) {
          size <- trait.data[sibs.list[q],] 
          size.dist <- append(size.dist, size)
        }
        size.dist <- as.numeric(size.dist[!is.na(size.dist)])
        
        # then get the range of body sizes from those sibling taxa
        ## and add uncertainty by expanding the range on each end (diff of max and min)
        ### again this if/else is set up to handle both groups of taxa (else) or singletons (if)
        min.size <- min(size.dist); max.size <- max(size.dist); difference <- max.size-min.size
        
        n <- length(input.tree$tip.label)
        terminal.edge.lengths <- setNames(input.tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=input.tree$edge[,2])],input.tree$tip.label)
        new.edge.length <- terminal.edge.lengths[species.list[t]] # edge length of the new tip
        target.node <- subset(edge.node.info, edge.node.info[,2] == outer.data)
        sister.branch.length <- max(nodeHeights(input.tree)) - target.node[,"finish"]
        print(sister.branch.length)
        
        if(min.size==max.size) {
          minie <- min.size - (min.size/10)
          maxie <- max.size + (max.size/10)
          size.range <- seq(minie, maxie, 0.05)
        } else {
          size.range <- seq(min.size-difference, max.size+difference, 0.1)
        }
        sim.size <- sample(size.range, 1) # finally pull a single estimate from the range (our simulated size data!) 
        output.data <- as.data.frame(t(c(paste(species.list[t]), sim.size)))
        new.species.data <- rbind(new.species.data, output.data)
      }
    new.species.data[,2] <- as.numeric(as.character(new.species.data[,2])); new.species.data[,1] <- as.character(new.species.data[,1])
    new.data <- subset(new.species.data, select = -c(V1))
    rownames(new.data) <- new.species.data$V1
    input.data <- rbind(trait.data, new.data)
    
    fossil.data[[z]] <<- input.data
  }
}

# test out the function to see if it works:
sim.fossil.data(phy=output.trees, num.taxa=20, trait.data = data) # output is called 'fossil.data'

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD

fitContinuous(output.trees[[1]], fossil.data[[1]], model="EB")
fitContinuous_paleo(output.trees[[1]], fossil.data[[1]], model="TRC", shift.time=7)
fitContinuous_paleo(output.trees[[1]], fossil.data[[1]], model="BM1.OU.BM1", shift.time=6, shift.time2=3)

fitContinuous_paleo(trees[[1]], data, model="TRC.time")
fitContinuous_paleo(trees[[1]], data, model="TRC", shift.time=5.26)
fitContinuous_paleo(trees[[1]], data, model="BM1.OU.BM1.time")

calc_AICc_vals(64.2227, 4, nTips(trees[[1]]))

#
























new.species.data <- NULL
for (t in 1:length(species)) {
  sibs.list <- NULL
  sis.node <- getSisters(output.trees[[1]], paste(species[t]), mode="label") # find what the sister taxon or node is to your new tip
  
  if (length(sis.node$tips[1])==1 && (sis.node$tips[1] %in% species)) {
    parent.node <- findMRCA(output.trees[[1]], c(paste(species[t]), paste(sis.node$tips[1])), type="node")
    sis.node <- Descendants(output.trees[[1]], parent.node, type="tips")
    for (i in 1:length(sibs[[1]])) {
      tip.num <- sibs[[1]][i]
      sibling <- output.trees[[1]]$tip.label[tip.num]
      sibs.list <- append(sibs.list, sibling)
    } 
  } else if (is.null(sis.node$tips)) { 
     sibs <- Descendants(output.trees[[1]], sis.node$nodes[1], type="tips") # from phangorn, get only the tips which descend from the node
      for (i in 1:length(sibs[[1]])) {
        tip.num <- sibs[[1]][i]
        sibling <- output.trees[[1]]$tip.label[tip.num]
        sibs.list <- append(sibs.list, sibling)
      } 
   } else {
    sibs.list <- sis.node$tips[1]
   }
  # now we can simulate data based on our real body sizes, start with a loop
  ## that pulls out empirical size data for all the sibling taxa
  size.dist <- NULL
  for (z in 1:length(sibs.list)) {
    size <- data[sibs.list[z],] 
    size.dist <- append(size.dist, size)
  }
  size.dist <- as.numeric(size.dist[!is.na(size.dist)])
  

  # then get the range of body sizes from those sibling taxa
  ## and add uncertainty by expanding the range on each end (diff of max and min)
  ### again this if/else is set up to handle both groups of taxa (else) or singletons (if)
  min.size <- min(size.dist); max.size <- max(size.dist); difference <- max.size-min.size
  if(min.size==max.size) {
    minie <- min.size - (min.size/10)
    maxie <- max.size + (max.size/10)
    size.range <- seq(minie, maxie, 0.05)
  } else {
    size.range <- seq(min.size-difference, max.size+difference, 0.1)
  }
  sim.size <- sample(size.range, 1) # finally pull a single estimate from the range (our simulated size data!) 
  output.data <- as.data.frame(t(c(paste(species[t]), sim.size)))
  new.species.data <- rbind(new.species.data, output.data)
}
new.species.data[,2] <- as.numeric(as.character(new.species.data[,2])); new.species.data[,1] <- as.character(new.species.data[,1])
new.data <- subset(new.species.data, select = -c(V1))
rownames(new.data) <- new.species.data$V1
data <- rbind(data, new.data)


# check to make sure the new tree and data match
name.check(tree, data); #check to make sure the tips match the data labels
# test it out:
fitContinuous(tree, data, model="BM")
fitContinuous_paleo(tree, data, model="TRC", shift.time=15)






sis.node <- getSisters(new.tree, "z", mode="label") # find what the sister taxon or node is to your new tip
sibs.list <- NULL
# this if/else pulls together a list of sibling taxa: either from groups (if) or singletons (else)
if(is.null(sis.node$tips)) {
  sibs <- Descendants(new.tree, sis.node$nodes[1], type="tips") # from phangorn, get only the tips which descend from the node
  for (i in 1:length(sibs[[1]])) {
    tip.num <- sibs[[1]][i]
    sibling <- new.tree$tip.label[tip.num]
    sibs.list <- append(sibs.list, sibling)
  } 
} else {
  sibs.list <- sis.node$tips[1]
}
# now we can simulate data based on our real body sizes, start with a loop
## that pulls out empirical size data for all the sibling taxa
size.dist <- NULL
for (z in 1:length(sibs.list)) {
  size <- data[sibs.list[z],] 
  size.dist <- append(size.dist, size)
}
# then get the range of body sizes from those sibling taxa
## and add uncertainty by expanding the range on each end (diff of max and min)
### again this if/else is set up to handle both groups of taxa (else) or singletons (if)
min.size <- min(size.dist); max.size <- max(size.dist); difference <- max.size-min.size
if(min.size==max.size) {
  minie <- min.size - (min.size/10)
  maxie <- max.size + (max.size/10)
  size.range <- seq(minie, maxie, 0.05)
} else {
  size.range <- seq(min.size-difference, max.size+difference, 0.1)
}
sim.size <- sample(size.range, 1) # finally pull a single estimate from the range (our simulated size data!) 






#### Now the Business, adding the tips
tree<-read.tree("best.mammalia.tre") #load your tree
species<-c(1:1500) #identify how many new tips to include, package it up as "species"
for(i in 1:length(species)) tree<-bind.at.depth(tree, "taxon",2) #loop that shit, it should loop it adding a new "taxon" 146 times at age 2my
plot.phylo(tree, show.tip.label=FALSE)
write.tree(tree, file="chiroptera.young.tre" )
