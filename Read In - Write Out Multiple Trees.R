library(ape)
library(phytools)

# read in all the files you want to handle
#(designate the folder, then isolate the files by their standard ending)
files <- list.files(pattern="_species_final_align.trees", recursive=FALSE)

#################################################################
# loop through all the tree files in the folder
# use the treeannotator wrapper to make a consensus for each tree
#################################################################
for(i in 1:length(files)){
  input = files[i]
  burnin = 50 #in BEAST2, this is a percentage, might need to adjust this script for BEAST1
  exe = "/Applications/BEAST_2.4.2/bin/treeannotator" #this is the location of treeannotator in your current BEAST directory
  
  filename <- tail(unlist(strsplit(input, "/")), 1)
  datapath <- gsub(filename, "", input)
  filename <- head(unlist(strsplit(input, "\\.")), 1)
  output <- paste(filename, "con", sep = ".")	
  
  burn.in <- paste("-b",burnin)
  call <- paste(exe, burn.in, input, output)
  system(call)
  tr <- read.nexus(output)
  tr
}



#################################################################
# loop through all the tree files in the folder
# use the logcombiner wrapper to make a burned-in set of trees
#################################################################
for(i in 1:length(files)){
  input = files[i]
  burnin = 50 #in BEAST2, this is a percentage, might need to adjust this script for BEAST1
  exe = "/Applications/BEAST_2.4.2/bin/logcombiner" #this is the location of logcombiner in your current BEAST directory
  
  filename <- tail(unlist(strsplit(input, "/")), 1)
  datapath <- gsub(filename, "", input)
  filename <- head(unlist(strsplit(input, "\\.")), 1)
  output <- paste(filename, "burned", sep = ".")	
  
  burn.in <- paste("-b",burnin)
  call <- paste(exe, burn.in, "-log", input, "-o", output)
  system(call)
  tr <- read.nexus(output)
  tr
}


# read in all the files you want to handle
#(designate the folder, the isolate the files by their standard ending)
pb.files <- list.files(pattern=".burned", recursive=FALSE)

#################################################################
# loop through all the post-burnin (pb) tree files in the folder
# pull out 100 random trees, and write them to a tree file you'll use as bootstraps
#################################################################
for(i in 1:length(pb.files)){
  trees <- read.nexus(pb.files[i])
  treex <- sample(trees, size=100)
  write.tree(treex, paste(pb.files[i], "_random100.tre", sep=""))
}



# read in all the files you want to handle
#(designate the folder, the isolate the files by their standard ending)
con.files <- list.files(pattern=".con", recursive=FALSE)

#################################################################
# loop through all the consensus gene tree files in the folder
# write each one to the all.consensus.genetrees file for ASTRAL
#################################################################
for(i in 1:length(con.files)){
  input = read.nexus(con.files[i])
  write.tree(input, file="all.consensus.genetrees", append=TRUE)
}