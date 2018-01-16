library(phytools) #uses ape as a dependency
library(ggtree)

#Plotting Trees with Support Values as Node Labels

tree <- read.tree("Sub.all.loci.tre") # read in your tree
tree <- midpoint.root(tree) # reroot it if necessary
tree <- ladderize(tree, right=T) # ladderize it for presentation
tree$node.label # check to make sure your node labels are there
plot(tree, cex=0.3); nodelabels(tree$node.label,node=2:tree$Nnode+Ntip(tree),
           adj=c(1,-0.2),frame="none") # plot it with node labels, adjust size as necessary


# if you want to remove a clade, but keep the node labels
plot(tree, cex=0.3); nodelabels(cex=0.3, frame="circle") # plot the tree with node numbers
oreophilus.group <- extract.clade(tree, node=947) # extract the clade of interest 
# alternatively you can go the setdiff route
write.tree(oreophilus.group, file="P.oreophilus.XXX.tre") # don't 'write.nexus' or it'll lose the node values


tree <- read.tree("/Users/Ian/Google.Drive/Villanova Herp Work/Assorted Geckos - Blaesodactylus, Lygodactylus, Dixonius & Others/Pachydactylus and Chondrodactylus/SubSaharan.ND2.RAG1.PDC.April.20.tre")
tree <- read.tree("P.oreophilus.XXX.tre")
drip <- c("Pachydactylus.maraisi.JV.1856", "Pachydactylus.otaviensis.MCZ.F.38512",
          "Pachydactylus.scutatus.CAS.254826.Iona", "Pachydactylus.scutatus.AMB.4041",
          "Pachydactylus.parascutatus.b", "Pachydactylus.parascutatus.AMB.7633", "Pachydactylus.parascutatus.a",
          "Pachydactylus.caraculicus.AG.50", "Pachydactylus.caraculicus.JET158.ARG321.SW.ANG", "Pachydactylus.caraculicus.JET157.ARG320.SW.ANG",
          "Pachydactylus.caraculicus.AMB10347.Dolondolo", "Pachydactylus.caraculicus.AG96.Namibe",
          "Pachydactylus.caraculicus.AMB10622.Virei", "Pachydactylus.caraculicus.AMB10626",
          "Pachydactylus.caraculicus.MCZ.R.185767", )
