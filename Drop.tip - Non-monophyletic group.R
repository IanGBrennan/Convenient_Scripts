library(phytools)

####################################################
# Keep a designated set of non-monophyletic tips
####################################################

tree<-read.tree("pygopodoidea.tre")
samples<-tree$tip.label #set the tip labels as a list

#set the non-mono tips you want to keep
targets<-c("Saltuarius.salebrosus",
           "Saltuarius.cornutus.B",
           "Saltuarius.wyberba",
           "Saltuarius.moritzi",
           "Saltuarius.swaini", 
           "Delma.australis", 
           "Lialis.burtonis")

save<-setdiff(samples, targets) #create a list of everything except your targets (setdiff)
missing<-drop.tip(tree, save) #drop all the other crap


