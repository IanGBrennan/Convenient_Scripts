library(ape)
library(phytools)
library(strap)
library(devtools)


tree<-read.tree("pygopodoidea.tre")
branching.times(tree) #find the age of the deepest node
tree$root.time<-56.7 #input the crown age of the tree (deepest node)
geoscalePhylo(tree=ladderize(tree, right=FALSE), label.offset=0, cex.age=0.4, cex.ts=0.4, cex.tip=0.4)

tree<-read.tree("pygopodoidea.tre")
sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
tree$root.time<-56.7
geoscalePhylo(tree=ladderize(tree, right=FALSE), label.offset=0, cex.age=0.4, cex.ts=0.4, cex.tip=0.2)

