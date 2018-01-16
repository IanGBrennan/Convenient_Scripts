library(ape)
library(phytools)

pygo<-read.tree("Pygopodoidea.tre")
spheno<-read.tree("Sphenomorphine.tre")
plot(g, cex=0.2)
write.tree(pygo, file="Full.tre")
branching.times(full)
plot(full, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="pink"); tiplabels(cex=0.2, frame="circle", bg="yellow")
plot(g, cex=0.2); edgelabels(cex=0.2, frame="circle", bg="pink")
gre<-reroot(g, 545, position=0)
plot(gre, cex=0.2)



##### This can not possibly be the easiest way to do this, but it's something
# we're going to combine trees by first changing a node in the target phylogeny to match the crown divergence of the arrow phylogeny
# then we're going to check the branching times, to determine the difference between where we inserted the new tree,
# and where it needs to be to make the tree ultrametric
# finally, we'll have to remove the taxa that were changed in the target tree, because they're no longer ultrametric
plot(full, cex=0.2)
combined<-bind.tree(full, spheno, where=546, position=0) #11.596020657
plot(combined, cex=0.2)
plot(combined, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="yellow");
branching.times(combined)
is.ultrametric(combined) #even though it's now placed properly, it's not yet UM because of the displaced taxa from above
is.ultrametric(spheno)
comb<-drop.tip(combined, c('Sphaerodactylus.roosevelti', 'Sphaerodactylus.torrei')) #remove the shitheads
is.ultrametric(comb) # check to make sure it's UM now
write.tree(comb, file="Pygo.Spheno.UM") #save the tree externally

siz<-read.tree("Pygo.Spheno.UM.tre")
plot(siz, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="yellow");tiplabels(cex=0.2, frame="circle", bg="yellow")
elapids<-read.tree("Oz.Elapidae.tre")
com<-bind.tree(siz, elapids, where=336, position=15.95389677)
plot(com, cex=0.2);nodelabels(cex=0.2, frame="circle", bg="yellow")
write.tree(com, file="Pygo.Spheno.Elap.tre")

branching.times(com)
py.sp.el<-drop.tip(com, c(158:246)) #remove Ctenotus
plot(py.sp.el, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="yellow");

branching.times(py.sp.el)

py.sp.el<-drop.tip(com, c(158:246))
write.tree(py.sp.el, file="Pygo.Spheno.Elap.tre")

pyspel<-read.tree("Pygo.Spheno.Elap.tre")
plot(pyspel, cex=0.3); tiplabels(cex=0.1, frame="circle", bg="yellow"); nodelabels(cex=0.1, frame="circle", bg="yellow");
getMRCA(pyspel, c("Teratoscincus.roborowskii", "Teratoscincus.scincus"))
noerem<-drop.tip(pyspel, c(108:115))
write.tree(noerem, "Pygo.Spheno.Elap.tre")

plot(pyspel, cex=0.3);tiplabels(cex=0.1, frame="circle", bg="yellow"); nodelabels(cex=0.1, frame="circle", bg="yellow");
is.ultrametric(pyspel)
is.binary.tree(pyspel)

#### Add a tip to the tree, bind your tree, and erase the added tip
getMRCA(pyspel, c("Sphenodon.punctatus", "Gallus.gallus"))
tip<-"A2"
sister<-"Demansia_psammophis"
tree<-bind.tip(pyspel,tip,where=262,position=0.5)

plot(tree, cex=0.2); tiplabels(cex=0.1, frame="circle", bg="yellow"); nodelabels(cex=0.1, frame="circle", bg="yellow")
pyspel<-drop.tip(pyspel, c("Teratoscincus.roborowskii", "Teratoscincus.scincus"))
is.ultrametric(pyspel)
is.binary.tree(pyspel)

python<-read.tree("Oz.Pythonidae.UM.tre")
branching.times(tree)
com<-bind.tree(tree, python, where=99, position=16.186894)
py.sp.el.py<-drop.tip(com, "A2")
plot(py.sp.el.py, cex=0.2);nodelabels(cex=0.2, frame="circle", bg="yellow")
write.tree(py.sp.el.py, file="Pygo.Spheno.Elap.Pyth.tre")

branching.times(com)
is.binary.tree(py.sp.el.py)
is.ultrametric(py.sp.el.py)

limbless<-read.tree("Australian.Limb.Reduced.tre")
getMRCA(limbless, c("Anilios_diversus", "Anilios_pilbarensis"))
blindies<-extract.clade(limbless, 7162)
plot(blindies)
write.tree(blindies, file="OZ.Blindies.tre")

##### Combining Ultrametric Trees, the better way
add<-read.tree("Pygo.Spheno.Elap.Pyth.tre") #load your traget tree
plot(add, cex=0.2);nodelabels(cex=0.2, frame="circle", bg="yellow");tiplabels(cex=0.2, frame="circle", bg="yellow")
tip<-"A2" #designate name for a new tip
sister<-c(50:113) # identify where you want to place the tip, this can be a node, a sister taxon (tip), or here, a group of tips
tree<-bind.tip(add,tip,where=277,position=32.687719) #bind the tip in place, give it the same age as the crown divergence of the add-in tree
plot(tree, cex=0.2); tiplabels(cex=0.1, frame="circle", bg="yellow"); nodelabels(cex=0.1, frame="circle", bg="yellow")
com<-bind.tree(tree, blindies, where=114, position=32.687719) #bind the added tree to the new tip, designate the crown age as the position
py.sp.el.py.an<-drop.tip(com, "A2") #get rid of the added tip, to make it binary again
plot(py.sp.el.py.an, cex=0.2);nodelabels(cex=0.2, frame="circle", bg="yellow")
py.sp.el.py.an<-drop.tip(py.sp.el.py.an, "Sphaerodactylus.glaucus")
write.tree(py.sp.el.py.an, file="Pygo.Spheno.Elap.Pyth.Anil.tre") # write the tree and have a look!
oz.limbless<-extract.clade(py.sp.el.py.an, 247)
write.tree(oz.limbless, file="Oz.Limbless.tre")
oz.limbless$tip.label
limbless<-read.tree("Oz.Limbless.tre")
limbless$tip.label
is.ultrametric(py.sp.el.py.an) #check to make sure it's ultrametric
is.binary.tree(py.sp.el.py.an) #check to make sure it's binary/branching

#### if you work from the top down, it will keep the same labels, as they count from the bottom up!
base<-read.tree("Base.Tree.tre")
pygo<-read.tree("pygopodidae.tre")
plot(base, cex=0.2);tiplabels(cex=2, frame="circle", bg="yellow")
com<-bind.tree(base, pygo, where=5, position=24.804503)
com<-drop.tip(com, "pygo")
plot(com, cex=0.2);tiplabels(cex=0.5, frame="circle", bg="yellow")
spheno<-read.tree("Oz.limbless.Spheno.tre")
coms<-bind.tree(com, spheno, where=4, position=24.2349953257)
plot(coms, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
com<-drop.tip(coms, "spheno")
blind<-read.tree("OZ.Blindies.tre")
com<-bind.tree(com, blind, where=3, position=32.687719)
elapid<-read.tree("Oz.Elapidae.tre")
com<-bind.tree(com, elapid, where=2, position=24.23461777)
python<-read.tree("Oz.Pythonidae.UM.tre")
com<-bind.tree(com, python, where=1, position=16.186894)
plot(com)
is.ultrametric(com)
is.binary.tree(com)
com<-drop.tip(com, c("elapid", "python", "aniliid"))
write.tree(com, file="Oz.Limbless.tre")

tree<-read.tree("Oz.Limbless.tre")
plot(tree, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
trial<-drop.tip(tree, c(17:64))
plot(trial, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
com<-bind.tree(trial, elap, where=16, position=24.23461777)
plot(com, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
final<-drop.tip(com, 16)
write.tree(final, file="New.Oz.Limbless.tre")

elap<-read.tree("Elapidae.tre")
plot(elap, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
elap$tip.label
elap<-drop.tip(elap, c("Aspidomorphus_muelleri","Aspidomorphus_lineaticollis","Aspidomorphus_schlegeli"))
write.tree(elap, file="Oz.Elapidae.tre")
x<-read.tree("Oz.Elapidae.tre")
x$tip.label

tree<-read.tree("New.Oz.Limbless.tre")
tree$tip.label
plot(tree, cex=0.2);nodelabels(cex=0.3, frame="circle", bg="yellow")
plot(tree, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
plot(tree, cex=0.2);edgelabels(cex=0.3, frame="circle", bg="yellow")
new<-bind.tip(tree, tip.label="Pseudonaja.guttata", where=289, position=1)
plot(new, cex=0.2);tiplabels(cex=0.3, frame="circle", bg="yellow")
write.tree(new, file="New.Oz.Limbless.tre")



