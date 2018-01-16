#if you don't already have it, install PhyTools:
install.packages('phytools')
library(phytools) #open up the package and let's get going

#######################################
# CoPhylo Plots for Visual Comparison #
#######################################

#first things first, make sure to set your working directory wherever your trees are
#you can read in your trees either as nexus or newick strings,
#just make sure that you do the same for both, or it jacks up the script

astral<-read.nexus("astral.consensus.tre") #read in your ASTRAL phylogeny
raxml<-read.nexus("raxml.concatenated.tre") #read in your concatenated RAxML phylogeny
#if 'read.nex' doesn't work, try out 'read.tree'
co<-cophylo(astral, raxml, rotate=TRUE) #create an object ('co') of the trees facing one another, with branches rotated to reduce crossed lines
plot(co, fsize=0.3, pts=FALSE) #plot the object, control font size and tipballs ('pts'), among other things here

#boom, that was easy. export it so you can fiddle in illustrator or whathaveyou.


#if you want to just look at the two trees plotted separately:
par(mfrow=c(1,2)) #set up a plot with two figures
#go back to 1 figure at a time with par(mfrow=c(1,1))
plot(astral) #plot first tree
plot(raxml) #plot second tree
