# steps for downloading and combining a genome from NCBI
basic outline here: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/

1. go to: https://www.ncbi.nlm.nih.gov/assembly/
  search for 






#source("https://bioconductor.org/biocLite.R")
#biocLite("BSgenome")
library(BiocInstaller)
biocLite("BSgenome.Celegans.UCSC.ce2")
biocLite("BSgenome.Pmolurus")
library(BSgenome)
library(bamsignals)
biocLite("snpStats")
install.packages("hoardeR")
library(hoardeR)

species
summary(species)
python <- getAnnotation(species="Python bivittatus",
                        annotationFolder="/Users/Ian/Desktop")


installed.genomes()
available.genomes()
genome <- getBSgenome("BSgenome.Celegans.UCSC.ce2")

hits <- read.csv("/Users/Ian/Downloads/NT1PJ4ZR015-Alignment-HitTable.csv", header=F)
hit.table.columns <- c("query_id", "subject_ids", "query_acc.ver",
                       "subject_acc.ver", "%_identity", "alignment_length",
                       "mismatches", "gap_opens", "q.start", "q.end",
                       "s.start", "s.end", "evalue", "bit_score")
colnames(hits) <- hit.table.columns


# bootleg version: create loop to subset all results for each input locus, to filter out suboptimal matches
best.matches <- NULL
unique.targets <- unique(hits$query_id)
for (i in 1:length(unique.targets)) {
  current.target <- subset(hits, hits$query_id == unique.targets[i])
  keep <- current.target[1,]
  best.matches <- rbind(best.matches, keep)
}




# create loop to subset all results for each input locus, to filter out suboptimal matches
best.matches <- NULL
unique.targets <- unique(hits$query_id)
for (i in 1:length(unique.targets)) {
  current.target <- subset(hits, hits$query_id == unique.targets[i])
  for (j in 1:length(current.target)) {
    current.test <- current.target[j,]
    if(current.test$)
  }
}


unique(hits$subject_ids)
