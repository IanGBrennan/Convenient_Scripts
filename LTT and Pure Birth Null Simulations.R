install.packages('TESS')
??TESS
??phytools
library(help=phytools)
library(phytools)
library(laser)
library(BAMMtools)
library(ape)
library(geiger)
tree<-read.tree('Clipped.Pygopodoidea.tre')
plot(tree)
ltt(tree, log=TRUE)
ltt(tree)
plot(ltt(tree, log=TRUE), xaxis="flipped")
tree$tip.label

#####################################################################################
#assess the Gamma statistics of the tree (tree shape, branching times, distribution of BTimes)
#####################################################################################
gammaStat(tree)
2*(1 - pnorm(abs(gammaStat(tree))))
#or
gamStat(branching.times(tree), return.list=TRUE) #fastest way

#compare distribution of your tree against differing models of evolution (compare aic values to establish best model fit)
pureBirth(branching.times(tree)) #under a pure birth model
DDX(branching.times(tree)) #density dependent: speciation rate decreases as the number of species in the clade, with exponential decay
DDL(branching.times(tree)) #density dependent: speciation rate decreases as the number of species in the clade, logarithmic decay
fitSPVAR(branching.times(tree), init=c(.3, .2, .1)) # an exponentially declining speciation rate through time and constant extinction
# or just compare them all by running em at the same time
fitdAICrc(branching.times(tree), modelset = c("pureBirth", "bd", "DDL", "DDX"), ints = NULL)



#fastest way to extract a clade
carpygo<-extract.clade(tree, node=171)
write.tree(carpygo, file="Carphodactylidae.Pygopodidae.tre")
diplo<-extract.clade(tree, node=240)
write.tree(diplo, file="Diplodactylidae.tre")

ltt(carpygo, log=TRUE)
plot(carpygo, font=2, cex=0.3); nodelabels(cex=0.3, frame="circle", bg="pink")
carpygo$tip.label

plot(carpygo)
ltt(diplo, log=TRUE)
diplo$tip.label

171-carpygo
240-diplodactylidae

pbtree<-pbtree(b=0.0938, d=0.00547, n=220, t=57, nsim=100, type='continuous', extant.only=TRUE)
pbtree<-pbtree(b=1, d=0, t=57, nsim=10, type='discrete', extant.only=TRUE, ape=FALSE)

#####################################################################################
# Create and Plot the Null Pure Birth model simulations, including a 95% confidence interval
#####################################################################################
pbtree<-pbtree(n=155, scale=56.7, nsim=1000, extant.only=TRUE) #simulated for all Pygopodoidea
ltt(pbtree)
ltt95(pbtree, log.lineages=TRUE)
plot(ltt95(pbtree, log=TRUE, mode="mean"), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped

pbtree<-pbtree(n=70, scale=51.4, nsim=5000, extant.only=TRUE) #simulated for Carpho+Pygo
ltt(pbtree)
ltt95(pbtree, log.lineages=TRUE)
plot(ltt95(pbtree, log=TRUE), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped

pbtree<-pbtree(n=99, scale=45.7, nsim=5000, extant.only=TRUE) #simulated for Diplodactylidae
ltt(pbtree)
ltt95(pbtree, log.lineages=TRUE)
plot(ltt95(pbtree, log=TRUE), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped

#####################################################################################
# Create and Plot the Null Birth-death model simulations (specified), under a constant BD model, including a 95% confidence interval
#####################################################################################
bdtree<-pbtree(b=0.1, d=0.0001, n=155, scale=57.3, nsim=1000, extant.only=TRUE)
#ltt(bdtree)
#bdnullCI<-ltt95(bdtree, log=TRUE)
plot(ltt95(bdtree, log=TRUE), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped
par(new=T)
ltt(tree, log=TRUE)

birthdeath(tree) #quickly get a rough estimate for speciation/extinction rates of your tree

mltt.plot(tree, bdnullCI, log="y")
mltt.plot(tree)


#### Model PB trees based on diversification rates from crown family groups #####
pb.t.r<-pbtree(b=0.07, d=0, scale=57, nsim=100, type='continuous', extant.only=TRUE)
t.r.95<-plot(ltt95(pb.t.r, log=TRUE), xaxis="flipped")

kj<-pbtree(b=0.15, d=0, t=57, nsim=10, type='continuous', extant.only=TRUE)
mltt.plot(kj, log="y")

?bd.tree

pygo.pb<-pbtree(b=0.1522, d=0, t=57, scale=57, nsim=100, type='continuous', extant.only=TRUE)
mltt.plot(pygo.pb, log="y")
pygo.95<-plot(ltt95(pygo.pb, log=TRUE), xaxis="flipped")

carpho.pb<-pbtree(b=0.096, d=0, t=57, nsim=1000, type='continuous', extant.only=TRUE)
carpho.95<-plot(ltt95(carpho.pb, log=TRUE), xaxis="flipped")

diplo.pb<-pbtree(b=0.1, d=0, t=57, nsim=1000, type='continuous', extant.only=TRUE)
diplo.95<-plot(ltt95(diplo.pb, log=TRUE), xaxis="flipped")


#### Diversity Dependent Check ####
library(DDD)
tree<-read.tree("Pygopodoidea.TESS.tre")
y<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
k<-2*(length(y))
res<-dd_ML(brts=y, initparsopt = c(0.1, 0.003, k, idparsopt=(1:4), cond=0))
ll<-res[,4]
AIC(-570.79, k=2)
