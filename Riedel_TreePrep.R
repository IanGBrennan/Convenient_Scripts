
pygo.tree <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Dating Speciation and Diversity/Trees/Pygopodoidea.total.100mil.1+2.tre")
taxon.list <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Riedel_TaxonList.csv", header=T)
keepers <- taxon.list$Name_in_Tree; length(keepers)
new.tree <- drop.tip(pygo.tree, setdiff(pygo.tree$tip.label, keepers)); Ntip(new.tree)
write.tree(new.tree, file="/Users/Ian/Google.Drive/ANU Herp Work/Riedel_TaxonList.tre")
