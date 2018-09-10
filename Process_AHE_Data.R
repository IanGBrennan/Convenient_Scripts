source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Extract_Data_from_AHE_labels.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Change_AHE_Taxa_Labels.R")

path = "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T474_Varanus/AHE_Data/Astral"
project.name = "Varanus"

# read in your tree of interest
tree <- read.tree(paste0(path, "/","AstralTree.tre"))

# extract the label information and write to a csv
label.out <- extract.labels(tree)
write.csv(label.out, file=paste0(path, "/", paste0(project.name, "_SampleInfo.csv")), row.names=F)

# intervening step: open the csv file you made above, add a new column if you want to
# change taxon names, call it 'new_genus_species'
# make any changes (fixed species names, haplotypes, etc)

# read you label information back into R
label.info <- read.csv("Neobatrachus_SampleInfo.csv")

# change taxa labels in alignments or trees according to your new information csv
setwd("/Users/Ian/Desktop/Clean_Alignments") # set your working directory
list.files(pattern=".tre", recursive=FALSE) # check to see what files you'll change
change.labels(label.info, "tre")


tree1 <- read.nexus("/Users/Ian/Desktop/Clean_Alignments/A1.tre")
tree2 <- read.nexus("/Users/Ian/Desktop/Clean_Alignments/R1.tre")
coplot <- cophylo(tree1, tree2, rotate=T)
plot(coplot, fsize=0.6, pts=F)
cophyloplot(tree1, tree2, fsize=0.2)





tr1<-pbtree(n=26,tip.label=LETTERS)
tr2<-pbtree(n=26,tip.label=sample(LETTERS))
obj<-cophylo(tr1,tr2)
plot(obj, fsize=0.2)





#
