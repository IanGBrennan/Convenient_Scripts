library(phangorn)
library(ips)

setwd("/Users/Ian/Documents/ANU_Finished/T203_Eulamprus/RAxML_2alleles/New_Renamed_Phased_Alignments")
files <- dir(getwd(), pattern=".phylip")

keep.names <- read.csv("Macro_TipLabels.txt", sep="\t", header=F)
keep.names <- as.vector(keep.names[,1])

# If you're using NEXUS alignments
for (i in 1:length(files)){
  int.align <- read.nexus.data(files[i])
  keep.align <- int.align[which(names(int.align) %in% keep.names)]
  write.nexus.data(keep.align, file=paste0("Extant_", files[i]), interleaved=F)
}

# If you're using PHYLIP alignments (make sure it's SPACES not TABS between taxon name and sequences)
for (i in 1:length(files)){
  int.align <- read.dna(files[i], format="sequential")
  keep.align <- int.align[which(rownames(int.align) %in% keep.names),]
  write.phy(keep.align, file=paste0("Reduced_", files[i], interleaved=F))
}

# If you're using FASTA alignments
for (i in 1:length(files)){
  int.align <- read.dna(files[i], format="fasta")
  keep.align <- int.align[which(rownames(int.align) %in% keep.names),]
  write.FASTA(keep.align, file=paste0("Reduced_", files[i]))
}


#########################################################################

# If you're using FASTA alignments and want to make a file of each locus for a given taxon
# I already have a better version of this function called "extract.targets" in "Condensing_Alignments.R" but it's for PHYLIP
all.locus.file <- function(directory.path, target.taxon, file.type=".fas", new.directory){
  setwd(directory.path)
  files <- dir(getwd(), pattern=file.type)
  system(paste("mkdir", new.directory))

  for (i in 1:length(files)){
    int.align <- read.dna(files[i], format="fasta")
    keep.align <- int.align[which(rownames(int.align) %in% target.taxon),]
    locus.name <- strsplit(files[i], "_")[[1]][1]
    full.name <- strsplit(files[i], file.type)[[1]][1]
    if(length(rownames(keep.align))==0){print(paste("taxon not found in alignment", locus.name)); next}
    rownames(keep.align) <- full.name
    write.FASTA(keep.align, file=paste0(new.directory, "/", locus.name, "_stripped.fas"))
  }
  cat.call <- paste("cat *_stripped.fas >> Combined_Loci.fas")
  system(paste("cd", new.directory, "&&", cat.call))
}

all.locus.file(directory.path = "/Users/Ian/Desktop/Gecko_AHE/Original_Single_line_FASTA/Cleaned_Alignments/Trimmed_Alignments/Copies",
               target.taxon = "I11974_LSUHC9463_Squamata_Gekkonidae_Gekko__subpalmatus_seq1", file.type=".fasta", new.directory="All_Loci")

testfas <- read.dna("/Users/Ian/Desktop/Geckomics_Alignments/TLR5_ex2_RELEC-190.fas", format="fasta")
keep <- testfas[which(rownames(testfas) == "Gekko_japonicus"),]
write.FASTA(keep, file="/Users/Ian/Desktop/Geckomics_Alignments/TLR5_TEST.fas")



#############################################################################
#### Choose only alignments which contain a set of taxa
setwd("/Users/Ian/Desktop/Gecko_Combined_Datasets/Combined_Aligned")
file.names <- dir(getwd(), pattern=".contree")
good.names <- NULL; bad.names <- NULL
for(i in 1:length(file.names)){
  curr.name <- file.names[[i]]
  curr.tree <- read.tree(file.names[[i]])
  if("Delma_butleri" %in% curr.tree$tip.label & "Lialis_burtonis" %in% curr.tree$tip.label & 
     "Aprasia_parapulchella" %in% curr.tree$tip.label & "Pygopus_nigriceps" %in% curr.tree$tip.label){print(paste(curr.name, "is good")); good.names <- append(good.names, curr.name)}
  else{print(paste(curr.name, "is bad")); bad.names <- append(bad.names, curr.name)}
}
good.names; bad.names

