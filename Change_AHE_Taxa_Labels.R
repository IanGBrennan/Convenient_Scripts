#################################################################
# loop through all the alignment files in the folder
# use perl to change the taxon names according to a CSV file (name.changes)
#################################################################
#setwd("/Users/Ian/Desktop/Goodnames_415_alignments") # set your working directory

change.labels <- function(extract.out, file.extension) {
  for(k in 1:nrow(extract.out)){
    current <- as.data.frame(extract.out[k,])
    old.name <- current$tip_label
    new.name <- current$new_genus_species
    extension <- paste0(" *.", file.extension)
    
    rename <- paste("perl -pi -w -e 's/", old.name, "/", new.name, "/g;'", extension, sep="")
    system(rename)
  }
}
