# source for function 'bind.at.depth' to choose where a new tip goes into the phylogeny
# source for function 'sim.fossil.trees' to simulate multiple trees with extinct tips


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

alt.bind.at.depth<-function(tree, tip.label, depth, min.depth, max.depth){
  tip.depth <- sample(seq(min.depth, max.depth,0.1),1)
  H <- nodeHeights(tree)
  age.from.root<-max(H)-depth # how far is the shift from the root
  ii <- intersect(which(H[,1]<age.from.root),which(H[,2]>age.from.root)) # which edges span the desired placement depth
  edges <- tree$edge[ii,2] # what are the edge numbers, based on the above suitable matches?
  depths <- H[ii,2]-age.from.root # what are the depths above the root, of each edge?
  jj <- sample(1:length(edges),1) # choose one of the suitable edges to paste the new tip to
  #mm <- sample(0.1:sample(1:depths[[jj]], 1), 1)
  #bind.tip(tree, tip.label, where=edges[jj], position=depths[jj], 
  #         edge.length=(depths[jj]/(sample(1:depths[jj],1))))
  out.tree <- bind.tip(tree, tip.label, where=edges[jj], position=depths[jj], 
           edge.length=0.01)
  new.tip <- which(out.tree$tip.label==tip.label)
  rownames(out.tree$edge) <- c(1:length(out.tree$edge[,1])) # give the tree edge frame rownames
  target.almost <- subset(out.tree$edge, out.tree$edge[,2]==new.tip) # pull out the ancestor and descendant nodes of the target edge
  target.edge <- as.numeric(rownames(target.almost)) # get the number of the target edge
  out.tree$edge.length[[target.edge]] <- out.tree$edge.length[[target.edge]]+((max(nodeHeights(out.tree))-nodeheight(out.tree, new.tip))-tip.depth) # add the desired length to the branch, per shift (here, 0.01)
  return(out.tree)
  # by making the new branch length = depth/1:depth, this allows the
  # new branch to reach the present (depth/1) or be quite short (depth/depth)
}
#get the depth of the "where", then get the difference to your target depth


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
  return(output.trees)
}

alt.sim.fossil.trees <- function (phy, time.frame, die, num.taxa, num.trees) {
  species<-as.character(c(1:num.taxa)) #identify how many new tips to include, package it up as "species"
  output.trees <<- NULL
  iterations = 1
  while (iterations <= num.trees) {
    tryCatch({
      time.seq <- seq(time.frame[1], time.frame[2], time.frame[3]) # pick a time range
      die.seq <- seq(die[1], die[2], die[3])
      # now a loop to add the desired number of taxa to the tree within the time period designated
      fossil.tree <- phy
      taxa.its <- 1
      while (taxa.its <= num.taxa) { 
        tryCatch({
          #cat(" tree", taxa.its ,"tip", i, "...")
          fossil.tree <- alt.bind.at.depth(fossil.tree, species[taxa.its], sample(time.seq, 1, replace=T), min(die.seq), max(die.seq))
          output.tree <<- fossil.tree
          taxa.its = taxa.its+1
        }, error=function(e){})
      }
    }, error=function(e){})
    output.tree$edge.length[output.tree$edge.length<0]<-0.1 # in case any branch lengths are negative
    output.trees[[iterations]] <<- output.tree
    iterations = iterations+1
  }
  class(output.trees) <<- "multiPhylo"
  return(output.trees)
}



