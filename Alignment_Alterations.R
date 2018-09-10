library(dplyr)
# annoying, but you can't parallelize the processes that change things by opening and adjusting it, before closing it
# because in the split second to open and change, another core can't find the file

# example: name.changes <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T222_Elapidae/Elapid_macroevolution/Cichlids/Cichlid_name_changes.csv", header = T)
# example: setwd("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T428_Microhylidae/Alignments_Reduced") # set your working directory

name.changes <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T474_Varanus/AHE_Data/Varanus_SampleInfo.csv",
                         header = T)
setwd("/Users/Ian/Desktop/Species_Tree/Filtering_Trees") # set your working directory


# read in all the files you want to handle
#(designate the folder, then isolate the files by their standard ending)
files <- list.files(pattern=".tre", recursive=FALSE)

#################################################################
# loop through all the alignment files in the folder
# use sed to remove the prefixes on the names (only if the prefix is consistent!)
#################################################################
path = "/Users/Ian/Desktop/Renamed_415_alignments"
out.file <- ""
file.names <- dir(path, pattern = ".phy")
for(i in 1:length(file.names)){
  in.file <- paste(path, "/", file.names[i], sep="")
  input = file.names[i]

  filename <- tail(unlist(strsplit(input, ".phy")), 1)
  datapath <- gsub(filename, "", input)
  filecall <- head(unlist(strsplit(input, "ENSONIE0000")), 2)[2]
  #output <- paste(path, "/", "Locus_", filecall, sep="")

  pre.call <- paste("sed -i '' 's/^[^'",filename, "_","']*'",filename, "_","'//'", sep="")

  call <- paste(pre.call, in.file)
  system(call)
}

#################################################################
# loop through all the alignment files in the folder
# use unix to rename all the file names (if there's a consistent convention!)
#################################################################
for(y in 1:length(file.names)) {
  #in.file <- paste(path, "/", file.names[y], sep="")
  input = file.names[y]
  
  filename <- tail(unlist(strsplit(input, ".phy")), 1)
  datapath <- gsub(filename, "", input)
  filecall <- head(unlist(strsplit(input, "ENSONIE0000")), 2)[2]
  #output <- paste(path, "/", "Locus_", filecall, sep="")
  output <- paste("Locus_", filecall, sep="")
  
  rename <- paste("mv", input, output)
  system(rename)
  
  #output <- paste("Locus_", filecall, sep="")
  #renamecall <- paste("rename -f 's/", input, "/", output,"'")
  #system(renamecall)
}

#################################################################
# loop through all the alignment files in the folder
# use perl to change the taxon names according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/RAxML_conSeqs") # set your working directory
for(k in 1:nrow(name.changes)){
  current <- name.changes[k,]
  old.name <- current[,"tip_label"] # change as necessary
  new.name <- current[,"new_genus_species"] #change as necessary
  
  rename <- paste0("perl -pi -w -e 's/", old.name, "/", new.name, "/g;' *.tre")
  system(rename)
}


#################################################################
# loop through all the alignment files in the folder (phylip)
# use perl to drop some taxa according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/Species_Tree") # set your working directory
#to.drop <- name.changes[,3][1:78]
to.drop <- dplyr::filter(name.changes, species_tree == "No")$new_genus_species
for(j in 1:length(to.drop)){
  drop.it <- to.drop[j]

  drop.sample <- paste("sed -i '' -e '/", drop.it, "/d' *.phylip", sep="")
  system(drop.sample)
}
### when you're done, remember to change the number of samples in each alignment
#### at the top of your phylip files!
setwd("/Users/Ian/Desktop/Reduced_415_alignments") # set your working directory
call <- paste("perl -pi -w -e 's/ 147 / 69 /g;' *.phy")
system(call)


#################################################################
# if your alignments are in fasta format, try this instead:
# use perl to drop some taxa according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/Reduced_415_alignments") # set your working directory
name.changes <- name.changes[,c("new_genus_species", "dating_analysis")]
int.drop <- filter(name.changes, dating_analysis=="No")
to.drop <- int.drop[,1]
for(j in 1:length(to.drop)){
  drop.it <- to.drop[j]
  
  drop.sample <- paste0("/usr/local/bin/sed -i '' -e '/", drop.it,"/,+1 d' *.fasta")
  system(drop.sample)
  # this will probably give an error message "can't read : No such file or directory",
  # just ignore it, the taxon removal should still work
}


#################################################################
# finally, create a RAxML gene tree for each locus
# by calling RAxML for each alignment (change settings as you see fit)
#################################################################
setwd("/Users/Ian/Desktop/All_Alignments") # set your working directory
alignments <- list.files(pattern=".phy", recursive=FALSE)

for (z in 1:length(alignments)) {
  cat("tree", z, "of", length(alignments), "\n") #keep track of what tree/loop# we're on
  current <- alignments[z]
  locus.name <- tail(unlist(strsplit(current, ".phy")), 1)
  
  tree <- paste("/Users/Ian/Desktop/raxmlHPC-PTHREADS-SSE3  -f a -T 8 -m GTRGAMMA -s ", current, " -x 12345 -N 100 -k -n ", locus.name)
  system(tree)
}



#################################################################
# loop through all the tree files in the folder
# use perl to change the taxon names according to a CSV file (name.changes)
#################################################################
setwd("/Users/Ian/Desktop/RAxML_conSeqs") # set your working directory
for(k in 1:nrow(name.changes)){
  current <- name.changes[k,]
  old.name <- current[,"tip_label"] # change as necessary
  new.name <- current[,"subgenera_labels"] #change as necessary
  
  rename <- paste0("perl -pi -w -e 's/", old.name, "/", new.name, "/g;' *.tre")
  system(rename)
}

rename <- system(paste0("perl -pi -w -e 's/", "Soterosaurus_Varanus_bitatawa_KU_320000", "/", "Philippinosaurus_Varanus_bitatawa_KU_320000", "/g;' *.tre"))


#################################################################
# loop through all the trees in the folder
# drop tips to make simpler trees
#################################################################
keep.names <- c("Outgroup_Xenosaurus_grandis_MVZ_1377786",
                "Polydaedalus_Varanus_niloticus_PEMR_19410",
                "Solosaurus_Varanus_spinulosus_123428",
                "Euprepiosaurus_Varanus_prasinus_123719",
                "Psammosaurus_Varanus_griseus_123677",
                "Soterosaurus_Varanus_salvator_123682",
                "Euprepiosaurus_Varanus_indicus_MVZ_238226",
                "Soterosaurus_Varanus_bitatawa_KU_320000",
                "Varanae_Varanus_panoptes_123680",
                "Odatria_Varanus_tristis_55388a")

for(y in 1:length(files)) {
  #in.file <- paste(path, "/", file.names[y], sep="")
  input = files[y]
  current.tree <- read.tree(input)
  
  new.tree <- drop.tip(current.tree, tip=setdiff(current.tree$tip.label, keep.names))
  write.tree(new.tree, file=input)
}
