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
 




testo <- read.nexus.data(files[1])
testa <- read.dna("/Users/Ian/Documents/ANU_Finished/T203_Eulamprus/RAxML_2alleles/Phased_Alignments/T203_L5.phylip", 
                  format="sequential")
rownames(testa)
int.test <- testa[which(rownames(testa) %in% keepers),]
write.phy(int.test, "XXX.phylip", interleave=F)
testa[,2]

    
testo <- read.tree("/Users/Ian/Desktop/MacroCascini.tre")
testo$node.label
nodelabels(t)
plot(tree)
testo$edge.length
testo$tip.label
View(tree)
