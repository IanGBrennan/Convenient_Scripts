library(PhySortR)
library(phylotools)

# Just a little TREE preparation:
#################################
# PhySortR doesn't like using multiphylo objects, and it only works with *newick* trees
# so we'll either read a multiphylo object, or read in each tree, and
# write out an individual newick tree per locus
#################################
# if you have your trees in a multiphylo object and just need to write them out to a file:
trees = read.nexus("Data.TREES/PB.Meliphagides.100.trees") #pick out our set of posterior trees
#trees <- lapply(trees, drop.tip, "A_Physignathus_cocincinus") #if you have to drop tips from the set of trees
#class(trees)<-"multiPhylo"
# Run a loop that writes each one to an individual file in a named folder
dir.create(PhySortR/Elapidae) #create folder for the trees
for (i in 1:length(trees)){
  name <- paste("/Users/Ian/Desktop/PhySortR/Elapidae/.",i,".tre", sep="")
  write.tree(trees[[i]], file=name)
} #you should now have a folder with 100 tree separate tree files

# if you have your trees, but they're nexus, it's a little more complicated
setwd("/Users/Ian/Desktop") # set your working directory
dir.create("PhySortR"); dir.create("PhySortR/Elapidae") # create a folder for the new trees
# then pull in all the files, and read-in as nexus, write-out as newick
path = "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/BEAST/Third Pass - 25 loci - /Results/Top 100 Loci" # make sure to designate you path appropriately!
target.path = "/Users/Ian/Desktop/PhySortR/Elapidae"
file.names <- dir(path, pattern ="_con.tre")
for (m in 1:length(file.names)) {
  tree.name <- paste(path, "/", file.names[m], sep="")
  locus.tree <- read.nexus(tree.name)
  name.file <- paste(target.path, "/", file.names[m], ".newick", sep="")
  write.tree(locus.tree, file = name.file)
}


########################################################################
#### Now we're past the annoying bits, we can focus on sorting our trees
########################################################################
# locate the directory you want to pull your trees from
#elapid.gene.trees <- "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/RAxML_conSeqs/RAxML_Seqs/"
elapid.s.beast.trees <- "/Users/Ian/Desktop/PhySortR/Elapidae/"

sortTrees(target.groups = ("Acanthophis, Pseudechis"),
          in.dir = elapid.s.beast.trees,
          #out.dir = out.trees,
          mode = "c",
          clades.sorted = "E",
          extension = "con.tre.newick",
          clade.exclusivity = 0.95)





