library(phytools)

##############################################
#### Pull a Set of Taxa out of a Larger Tree #####
##############################################
tree<-read.tree("Gekkota.families.tre")
plot(tree, cex=0.3); tiplabels(cex=0.2, frame="circle") #grab the tip numbers (easier than using names)
species<-tree$tip.label[c(205, 244, 102, 60, 52, 28, 22:20, 18, 15, 5:3)] #vector our requested tips
pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, species));
plot(pruned.tree)
