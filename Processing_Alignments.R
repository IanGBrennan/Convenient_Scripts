library(tidyverse); library(RColorBrewer)
source("~/Google.Drive/R.Analyses/Convenient Scripts/ggplotRegression.R")

#locus_summary <- read.table("/Users/Ian/Desktop/Species_Tree/FINAL_Thesis_Sampling/TOTAL_EVIDENCE_small_NEXUS/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- read.table("/Users/Ian/Documents/ANU_Finished/T203_Eulamprus/RAxML_2alleles/New_Renamed_Phased_Alignments/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- read.table("~/Desktop/UCE_Alignments/Trimmed_Alignments/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- read.table("~/Desktop/Species_Tree/FINAL_FullSampling/Molecular_Alignments/Phylip_Trimmed/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- read.table("~/Desktop/GenomeStripper/Elapids/Existing_Alignments/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- read.table("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/Renamed_Alignments/Trimmed_Alignments/PerLocus_summary.txt", sep="\t", header=T)
#locus_summary <- readRDS("/Users/Ian/Desktop/Species_Tree/FINAL_Thesis_Sampling/TOTAL_EVIDENCE_small_NEXUS/PerLocus_summary_CLUSTER.RDS")
#locus_summary <- read.table("~/Google.Drive/ANU/AHE/T392_Neobatrachus/Combine_Projects/Combined_Final_Sampling_FASTA/REaligned/PerLocus_summary.txt", sep="\t", header=T)
locus_summary <- read.table("~/Google.Drive/ANU/AHE/T392_Neobatrachus/Diploid_Trees/Combined_Alignments/PerLocus_summary.txt", sep="\t", header=T)

# turn the table into a tibble
lsum <- as.tibble(locus_summary)

# check a couple common metrics
max(lsum$No_of_taxa); ggplot(lsum, aes(x=No_of_taxa)) + geom_density() + theme_classic()
max(lsum$Alignment_length); min(lsum$Alignment_length); ggplot(lsum, aes(x=Alignment_length)) + geom_density() + theme_classic()

#arrange(lsum, desc(Alignment_length))
lsum <- filter(lsum, Alignment_length > 500)
lsum <- filter(lsum, No_of_taxa > 20)
#lsum_sort <- arrange(lsum, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(AT_content))
#lsum_sort <- arrange(lsum, cluster.no, desc(No_of_taxa), desc(Proportion_parsimony_informative), desc(AT_content))
#lsum_sort <- arrange(lsum, desc(No_of_taxa), desc(No_variable_sites))
lsum_sort <- arrange(lsum, desc(No_variable_sites), desc(No_of_taxa))
#lsum_sort <- arrange(lsum, desc(No_of_taxa), Missing_percent)

# identify and color a subset of the loci (say the top 10)
#target.no <- 50; point.colors <- c(rep("#1F78B4",target.no), rep("#33A02C", (nrow(lsum_sort)-target.no)))
#colorz <- brewer.pal(4,"Spectral"); point.colors <- c(rep(colorz[1],20), rep(colorz[2],20), rep(colorz[3],20), rep(colorz[4],314))
clz <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- clz(nrow(lsum_sort))

View(lsum_sort)

# Visualize the relationships among variables of interest
no.var <- ggplotRegression(lm(lsum_sort$No_variable_sites ~ lsum_sort$Alignment_length), alpha=0.75, size=1, color=point.colors)
#p.var <- ggplotRegression(lm(lsum_sort$Proportion_variable_sites ~ lsum_sort$Alignment_length), alpha=0.75, size=1, color=point.colors)
at.var <- ggplotRegression(lm(lsum_sort$No_variable_sites ~ lsum_sort$AT_content), alpha=0.75, size=1, color=point.colors)

tax.len <- ggplotRegression(lm(lsum_sort$Alignment_length ~ lsum_sort$No_of_taxa), alpha=0.75, size=1, color=point.colors)
tax.var <- ggplotRegression(lm(lsum_sort$No_variable_sites ~ lsum_sort$No_of_taxa), alpha=0.75, size=1, color=point.colors)

no.par <- ggplotRegression(lm(lsum_sort$Parsimony_informative_sites ~ lsum_sort$Alignment_length), alpha=0.75, size=1, color=point.colors)
#p.par <- ggplotRegression(lm(lsum_sort$Proportion_parsimony_informative ~ lsum_sort$Alignment_length), alpha=0.75, size=1, color=point.colors)
at.par <- ggplotRegression(lm(lsum_sort$Parsimony_informative_sites ~ lsum_sort$AT_content), alpha=0.75, size=1, color=point.colors)

gridExtra::grid.arrange(no.var, at.var,
                        no.par, at.par, 
                        tax.len, tax.var, nrow=3)

# setwd("~/Desktop/Species_Tree/FINAL_FullSampling/Molecular_Alignments/Nexus_Full")
# setwd("~/Google.Drive/ANU/AHE/T392_Neobatrachus/Combine_Projects/Combined_Final_Sampling_FASTA/REaligned")
setwd("~/Google.Drive/ANU/AHE/T392_Neobatrachus/Diploid_Trees/Combined_Alignments")

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
#### file.ext: extension of the files, you could use this for trees if you've already inferred them!
pick.best(ranking=lsum_sort, loci.by.rank = 1:25, dir.name = "Top25_Loci", trimmed.alignments = TRUE, file.ext=".fasta") # picking alignments
#pick.best(ranking=lsum_sort, loci.by.rank = 1:100, dir.name = "Top100_Loci_Trees", trimmed.alignments = FALSE, file.ext="_trimmed.phy.contree") # picking trees
