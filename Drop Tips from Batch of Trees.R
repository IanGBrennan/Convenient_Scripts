library(ape)
library(phytools)


trees<-read.nexus("mil.trees") #read in your 1000 trees
class(trees)<-"multiPhylo" #make the trees a "multiPhylo" object, necessary
new.trees<-lapply(trees, drop.tip, tip=c("Gallus.gallus")) #create a function that applies the drop.tip bit to all the trees, separate by commas
class(new.trees)<-"multiPhylo" #set your new trees as a multiphylo object
write.nexus(new.trees, file="new.trees.nex") #save them externally
plot(new.trees[[1]]) #plot any single tree, designate the number within [[x]]



#################################################################
# What if you want to make a batch match a single tree you've already got?
#################################################################

input.trees <- read.tree("Pygopodoidea.1.burned_random100.tre") #input your multiphylo
samples <- input.trees[[1]]$tip.label #get the names of all samples in the multiphylo trees

single <- read.tree("Pygopodoidea.tre") #read in the single tree you want the multiphylo to match
targets <- as.list(single$tip.label) #get the names of all the samples you want to keep


remove <- setdiff(samples, targets) #create a list of everything except your targets (setdiff)
new.trees<-lapply(input.trees, drop.tip, tip=remove) #drop all the other crap
new.trees[[1]] # check it out!

class(new.trees) <- "multiPhylo" # change the class of the list to a multiPhylo object
write(new.trees, file="XXX.trees")
