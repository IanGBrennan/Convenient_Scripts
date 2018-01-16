library(geiger)
library(ape)
library(BAMMtools)
library(phytools)


#####################################################################
#### Plot Distribution of Speciation Events (not all branching events)
#####################################################################
mar<-read.tree("oz.marsupials.tre")
n<-length(mar$tip.label)
ms<-setNames(mar$edge.length[sapply(1:n,function(x,y)   which(y==x),y=mar$edge[,2])],mar$tip.label)
hist<-hist(ms, breaks=40, xlab="Branching Times", col="lightpink", ylim=c(0,30), xlim=c(40,0)) #breaks determines the # of bins distributed across the whole xlim=
multiplier<-hist$counts/hist$density
mydensity<-density(ms) #pull the speciation frequencies out
mydensity$y<-mydensity$y*multiplier[1]
lines(mydensity) #plot the smoothed-out line of best fit across our histogram
