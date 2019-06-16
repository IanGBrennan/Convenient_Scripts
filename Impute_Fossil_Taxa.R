library(phytools)

#######################################################
# Attach extinct tips to your empirical phylogeny
######################################################
base.tree <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Trees/Macro_AgesRanges_CON.tre")
## simulate trees under the "stochastic extinction" concept
#sim.extinct <- sim.fossil.trees(phy=base.tree, time.frame=c(0.1,max(nodeHeights(base.tree)),0.1), 
#                                num.taxa=(length(base.tree$tip.label))/2, num.trees=3) # output is called 'output.trees'

fossils <- data.frame(taxa=c("sp.A", "sp.B", "sp.C"), Min_age=c(5,10,8), Max_age=c(13,15,12),
                          clade_sample1=c("Dorcopsis_hageni", "Thylogale_thetis", "Wallabia_bicolor"),
                          clade_sample2=c("Dorcopsoides_fossilis", "Dendrolagus_dorianus", "Macropus_rufus"))


sim.extinct <- alt.sim.fossil.trees(phy=base.tree, time.frame=c(10,max(nodeHeights(base.tree)),0.1),
                                    die=c(5,10,0.1), 
                                    num.taxa=(length(base.tree$tip.label))*.50, num.trees=10) # output is called 'output.trees'

node<-124
tt<-splitTree(input.tree,split=list(node=node,
                              bp=input.tree$edge.length[which(input.tree$edge[,2]==node)]))
#tt[[2]]<-add.random(tt[[2]],tips="tip to add")
H <- nodeHeights(tt[[2]])
age.from.root <- max(H)-depth
ii
new.tree<-paste.tree(tt[[1]],tt[[2]])
plotTree(new.tree)

# I REALIZED I COULD ADD A AN ADDITIONAL COLUMN FOR THE OLDEST FOSSIL, SO THAT IT HAS TO SPLIT FROM THE TREE BEFORE AT LEAST THEN!
impute.fossils <- function(fossil.list, input.tree){
  for(i in 1:nrow(fossil.list)){
    # split the tree to isolate the clade of interest (target.clade[[2]])
    node <- getMRCA(input.tree, c(as.character(fossil.list[i,"MRCA1"]), as.character(fossil.list[i,"MRCA2"])))
    target.clade <- splitTree(input.tree, split=list(node=node, bp=input.tree$edge.length[which(input.tree$edge[,2]==node)])) # this makes the cut at the crown of the clade of interest
    other.clade <- splitTree(input.tree, split=list(node=node, bp=0)) # makes the cut up at the start of the root edge

    # identify the node heights for that clade
    H <- nodeHeights(target.clade[[2]])
    I <- nodeHeights(other.clade[[2]], root.edge=TRUE)

    # now make some stipulations about the maximum and minimum age of the fossil you're adding
    if (fossil.list[i,"Max_age"] >= max(H)){fossil.list[i,"Max_age"] = max(I)-0.1}
    #if (fossil.list[i,"Min_age"] <= min(H)){fossil.list[i,"Min_age"] = min(H)+0.01}
    if (fossil.list[i,"Min_age"] >= max(H)){fossil.list[i,"Min_age"] = max(H)-0.1}
    
    # set the depth of the tip (extinction time) of the new taxon
    tip.depth <- runif(1, fossil.list[i,"Min_age"], fossil.list[i,"Max_age"])

    # set the depth that the new taxon will join the tree
    depth <- runif(1, tip.depth+0.01, max(I)-0.01)
    
    # if the fossil speciation age (depth) is less than the maximum age of the target clade, add it in there
    if (depth < max(H)){
      age.from.root<-max(H)-depth # how far is the shift from the root
      ii <- intersect(which(H[,1]<age.from.root),which(H[,2]>age.from.root)) # which edges span the desired placement depth
      edges <- target.clade[[2]]$edge[ii,2] # what are the edge numbers, based on the above suitable matches?
      depths <- H[ii,2]-age.from.root # what are the depths above the root, of each edge?
      jj <- sample(1:length(edges),1) # choose one of the suitable edges to paste the new tip to
      target.clade[[2]] <- bind.tip(target.clade[[2]], as.character(fossil.list[i,"Taxon"]), where=edges[jj], 
                                    position=depths[jj], edge.length=depth-tip.depth)
      input.tree <- paste.tree(target.clade[[1]], target.clade[[2]])
    } else if (depth > max(H)){
      above.crown <- depth - max(H)
      target.clade[[1]] <- bind.tip(target.clade[[1]], as.character(fossil.list[i,"Taxon"]), 
                                    where=which(target.clade[[1]]$tip.label=="NA"),
                                    position=above.crown, edge.length=depth - tip.depth)
      input.tree <- paste.tree(target.clade[[1]], target.clade[[2]])
    } # if the fossil speciation age (depth) is greater than the maximum age of the target clade, add it somewhere to the root edge (stem)
    #print(as.character(fossil.list[i,"Taxon"]))
    #plot(input.tree); i<-i+1
  }
  return(input.tree)
}

testo <- impute.fossils(fossils, sorted.trees[[1]]); plot(testo)

trees <- sorted.trees[1:10]


#trees <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Trees/Macro_AgesRanges_NEWICK.trees")
#trees <- trees[990:1000]
#


impute.fossils.multiphylo <- function(fossil.list, input.trees){
  new.trees <- list()
  iterations = 1
  while (iterations <= length(input.trees)) {
    tryCatch({
      new.trees[[iterations]] <- impute.fossils(fossil.list = fossils, input.tree = input.trees[[iterations]])
    }, error=function(e){})
    iterations = iterations + 1
  }
  return(new.trees)
}


for (j in 1:length(trees)){
  new.trees[[j]] <- impute.fossils(fossil.list = fossils, input.tree = trees[[j]])
  print(j)
}

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



