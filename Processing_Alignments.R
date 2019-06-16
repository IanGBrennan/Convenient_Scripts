library(tidyverse)

#locus_summary <- read.table("/Users/Ian/Desktop/Species_Tree/FINAL_Thesis_Sampling/TOTAL_EVIDENCE_small_NEXUS/PerLocus_summary.txt", sep="\t", header=T)
locus_summary <- read.table("/Users/Ian/Documents/ANU_Finished/T203_Eulamprus/RAxML_2alleles/New_Renamed_Phased_Alignments/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- readRDS("/Users/Ian/Desktop/Species_Tree/FINAL_Thesis_Sampling/TOTAL_EVIDENCE_small_NEXUS/PerLocus_summary_CLUSTER.RDS")

lsum <- as.tibble(locus_summary)

#arrange(lsum, desc(Alignment_length))
lsum_sort <- arrange(lsum, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(AT_content))
#lsum_sort <- arrange(lsum, cluster.no, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(AT_content))
lsum_sort <- arrange(lsum, desc(No_of_taxa), Missing_percent)
View(lsum_sort)

setwd("/Users/Ian/Desktop/Geckomics_Alignments/Cleaned_Alignments")

pick.best <- function(ranking = NULL, loci.by.rank = 1:20, dir.name = "top_twenty", trimmed.alignments = FALSE, file.ext = ".nex"){
  if(is.null(ranking)) {stop("need a ranked data frame or tibble")}
  dir.path <- getwd(); call <- paste("mkdir", paste0(dir.path, "/", dir.name)); system(call)
  
  if(trimmed.alignments==FALSE){
    for(i in 1:length(loci.by.rank)){
      locus.name <- str_split(ranking$Alignment_name[i], "_")[[1]][1:2]
      locus.name <- paste0(locus.name[1], "_", locus.name[2], file.ext) # might have to add file extension!
      pre.call <- paste("cp", paste0(getwd(),"/",locus.name), paste0(dir.path, "/", dir.name))
      system(pre.call)
    }
  } else if(trimmed.alignments==TRUE){
    for(i in 1:length(loci.by.rank)){
      locus.name <- ranking$Alignment_name[i] # might have to add file extension!
      pre.call <- paste("cp", paste0(getwd(),"/",locus.name), paste0(dir.path, "/", dir.name))
      system(pre.call)
    }
  }

}
# EXPLANATION OF THE ARGUMENTS:
#### ranking: just the sorted tibble of our loci ranked from best (1) to worst (n)
#### loci.by.rank: how many loci do you want, specify which ones here 1:20 gives you the top 20, 21:40 would give you the second twenty
#### dir.name: what do you want to call the directory the alignments will be moved to?
#### trimmed.alignments: I can't remember why I put this in...
pick.best(ranking=lsum_sort, loci.by.rank = 1:10, dir.name = "Top10_Loci", trimmed.alignments = TRUE, file.ext=".fas.contree")
