library(phytools)

#######################################################
#### Create the "Bind at Depth" function for use later
#######################################################
#this is the source code, and must be run beforehand
bind.at.depth<-function(tree, tip.label, depth){
  H<-nodeHeights(tree)
  h<-max(H)-depth
  ii<-intersect(which(H[,1]<h),which(H[,2]>h))
  edges<-tree$edge[ii,2]
  depths<-H[ii,2]-h
  jj<-sample(1:length(edges),1)
  tree1<-bind.tip(tree, tip.label, where=edges[jj],position=depths[jj])
  bind.tip(tree, tip.label, where=edges[jj],position=depths[jj])
}

#### Now the Business, adding the tips
tree<-read.tree("best.mammalia.tre") #load your tree
species<-c(1:1500) #identify how many new tips to include, package it up as "species"
for(i in 1:length(species)) tree<-bind.at.depth(tree, "taxon",2) #loop that shit, it should loop it adding a new "taxon" 146 times at age 2my
plot.phylo(tree, show.tip.label=FALSE)
write.tree(tree, file="chiroptera.young.tre" )


#####################################################################
#### Plot Distribution of Speciation Events (not all branching events) of the modified, younger phylogeny
#####################################################################
tree<-read.tree("best.mammalia.tre") 
n<-length(tree$tip.label)
t<-setNames(tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)
hist<-hist(t, breaks=80, xlab="Branching Times", col="lightpink", ylim=c(0,500), xlim=c(40,0)) #breaks determines the # of bins distributed across the whole xlim=
multiplier<-hist$counts/hist$density
mydensity<-density(t) #pull the speciation frequencies out
mydensity$y<-mydensity$y*multiplier[1]
lines(mydensity) #plot the smoothed-out line of best fit across our histogram
abline(v=mean(t), col="blue", lwd=2) #add a line for the mean
abline(v=median(t), col="red", lwd=2) #add a line for the median
median(t)

