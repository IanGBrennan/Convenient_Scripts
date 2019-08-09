library(PhySortR)

# Just a little TREE preparation:
#################################
# PhySortR doesn't like using multiphylo objects, and it only works with *newick* trees
# so we'll either read a multiphylo object, or read in each tree, and
# write out an individual newick tree per locus
#################################
# if you have your trees in a multiphylo object and just need to write them out to a file:
trees = read.nexus("PATH_to_DATA/SOME.trees") #pick out our set of posterior trees
#trees <- lapply(trees, drop.tip, "A_Physignathus_cocincinus") #if you have to drop tips from the set of trees
#class(trees)<-"multiPhylo"
# Run a loop that writes each one to an individual file in a named folder
dir.create(PhySortR/Elapidae) #create folder for the trees
for (i in 1:length(trees)){
  name <- paste("/PAHT_TO_DIRECTORY/.",i,".tre", sep="")
  write.tree(trees[[i]], file=name)
} #you should now have a folder with 100 tree separate tree files

# if you have your trees, but they're nexus, it's a little more complicated
setwd("/Users/Ian/Desktop") # set your working directory
dir.create("PhySortR"); dir.create("PhySortR/Elapidae") # create a folder for the new trees
# then pull in all the files, and read-in as nexus, write-out as newick
path = "/PATH_TO_DIRECTORY" # make sure to designate you path appropriately!
target.path = "/PATH_TO_DIRECTORY/DIR"
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
target.trees <- "/PATH_DIR/Filtering_Trees"
sortTrees(target.groups = ("Psammosaurus,Polydaedalus,Soterosaurus,Euprepiosaurus,Varanae,Odatria"),
          in.dir = target.trees,
          #out.dir = out.trees,
          mode = "c",
          clades.sorted = "E",
          extension = ".tre",
          clade.exclusivity = 0.95)

# Let's add the cluster information to our per locus summary
library(tidyverse); library(dplyr)

# read in your summary info from AMAS
locus_summary <- read.table("/PATH_to_FILE/PerLocus_summary.txt", sep="\t", header=T)
#lsum <- as.tibble(locus_summary)

# write the table as an RDS file for R
saveRDS(locus_summary, "/PATH_to_FILE/PerLocus_summary_CLUSTER.RDS")

# sort the matrix by whatever columns
locus_summary <- arrange(locus_summary, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(Alignment_length))


test.sum <- readRDS("/PATH_TO_DIRECTORY/PerLocus_summary_CLUSTER.RDS")

test.sum <- arrange(test.sum, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(Alignment_length))
head(test.sum)

