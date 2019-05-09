library(phangorn)

setwd("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/Alignments/Extended_Macropodoidea_Complete/Full_Nexus")
files <- dir(getwd(), pattern=".nex")

keep.names <- read.csv("Macro_TipLabels.txt", sep="\t", header=F)
    keep.names <- as.vector(keep.names[,1])

for (i in 1:length(files)){
  int.align <- read.nexus.data(files[i])
  keep.align <- int.align[which(names(int.align) %in% keep.names)]
  write.nexus.data(keep.align, file=paste0("Extant_", files[i]), interleaved=F)
}

    
testo <- read.tree("/Users/Ian/Desktop/MacroCascini.tre")
testo$node.label
nodelabels(t)
plot(tree)
testo$edge.length
testo$tip.label
View(tree)
