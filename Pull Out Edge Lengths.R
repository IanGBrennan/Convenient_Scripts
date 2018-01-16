library(phytools)
library(geiger)
library(ape)
tree<-read.tree("Pygopodoidea.tre")
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)

write.csv(ee, file="Pygopodoidea.edgelenths.csv")
plot(tree, cex=0.3); edgelabels(cex=0.2, frame="circle")
plot(tree, cex=0.3); tiplabels(cex=0.2, frame="circle")
plot(tree, cex=0.3); nodelabels(cex=0.2, frame="circle")
plot(tree, cex=0.3); edgelabels(tree$edge.length, cex=0.3)


cc<-tree$edge.length[c(303:308, 299,298,296,295,291,290)] #remove 5 units from terminal edges 303-308
cf<-(tree$edge.length[c(303:308, 299,298,296,295,291,290)]-5) #remove 5 units from terminal edges 303-308

ch<-tree$edge.length[140]
plot(tree); edge

length(tree$tip.lable[145])

ke<-(ee["Saltuarius.salebrosus"]-5)
minfive<-(ee-5)
plot(minfive)


gammaStat(tree)
pygogamma<-2*(1-pnorm(abs(gammaStat(tree))))
pygogamma1<-1-pnorm(abs(gammaStat(tree)))

pygo<-read.tree("pygopodidae.tre")
gammaStat(tree)
pygamma<-2*(1-pnorm(abs(gammaStat(pygo))))
pygamma1<-1-pnorm(abs(gammaStat(pygo)))


require(geiger)
geo<-get(data(geospiza))
lphy<-rescale(geo$phy, "lambda", 0.5)
lphy <- rescale(geo$phy, "lambda", 0.5)
lphy <- geiger::rescale(geo$phy, "lambda", 0.9)

tree<-read.tree("oz.myo.tre")
times<-branching.times(tree)
short<-geiger::rescale(tree, "depth", 105)
write.tree(short, file="oz.myo.rescaled")

rescale<-function(tree,scale){
  tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
}
branching.times(rescale)
tree$edge.length<-tree$edge.length/max(nodeHeights(tree)[,2])*scale
