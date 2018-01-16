library(phylometrics)

tree<-read.tree("pygopodoidea.tre")
d<-read.csv("treestat.arid.csv")
arid<-d$arid #call your data column, data must be "0" (not focal trait), or "1" (focal trait)


treestat(tree, state=arid, func=tars, traitevol="TBM", a=1000, alternative="greater", simplify=TRUE)
