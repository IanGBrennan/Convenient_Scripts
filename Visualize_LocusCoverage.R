# alignments <- list.files("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Treefiles",
#                         pattern = "fasta_con.tre")
# setwd("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Treefiles")

alignments <- read.tree("/Users/Ian/Google.Drive/ANU/AHE/T392_Neobatrachus/Data_Consensus/T392_gene.trees")

all.names <- NULL
for (i in 1:length(alignments)){
  curr.tree <- alignments[[i]]
  all.names <- append(all.names, curr.tree$tip.label)
}
counts <- as.data.frame(table(all.names)); colnames(counts) <- c("Taxon", "Locus_Coverage")
counts$Taxon <- as.character(counts$Taxon)
library(ggplot2)
library(RColorBrewer)

short.name <- NULL
for (i in 1:nrow(counts)){
  splitted <- strsplit(counts[i,"Taxon"], "_")
  reordername <- paste0(splitted[[1]][3],"_",splitted[[1]][4],"_",
                        splitted[[1]][2])
  short.name <- append(short.name, reordername)
}
counts$New_Taxon <- short.name

#counts <- read.csv("Locus_counts.csv", header=T)

alphabetical <- ggplot(data=counts, aes(x=reorder(New_Taxon,desc(New_Taxon)), y=Locus_Coverage)) +
  geom_bar(stat="identity", aes(fill=cut(counts$Locus_Coverage, c(150,200,250,300,350,400,450))), show.legend = F)+ 
  theme_classic() +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #        panel.background = element_blank(), axis.line = element_line(colour = "black"),
  #        axis.text.x=element_text(angle=90, hjust=1)) +
  scale_fill_manual(values=brewer.pal(6, "YlGnBu")) +
  coord_flip() +
  geom_hline(yintercept = 100, linetype="dotted") +
  geom_hline(yintercept = 200, linetype="dotted") +
  geom_hline(yintercept = 300, linetype="dotted") 
  #scale_fill_manual(values=brewer.pal(11, "Paired"))

counts <- counts[order(counts$Subgenus),]
counts$Subgenus <- factor(counts$Subgenus, levels=unique(counts$Subgenus))
counts$Taxon <- factor(counts$Taxon, levels=counts$Taxon[order(counts$Subgenus)])

subgeneral <- ggplot(data=counts, aes(x=reorder(Taxon, desc(Taxon)), y=Locus_Coverage, fill=Subgenus)) +
  geom_bar(stat="identity") + 
  theme_classic() + scale_fill_manual(values=brewer.pal(11, "RdYlBu")) + coord_flip() + 
  geom_hline(yintercept = 100, linetype="dotted") +
  geom_hline(yintercept = 200, linetype="dotted") +
  geom_hline(yintercept = 300, linetype="dotted") 

gridExtra::grid.arrange(subgeneral, alphabetical, nrow=1)

ggplot(data=counts)+
  geom_bar(mapping = aes(x=Taxon, fill=Locus_Coverage))

ggplot(data=counts, aes(x=Taxon, y=Locus_Coverage)) +
  geom_bar(stat="identity")

library(RColorBrewer)
scales::show_col(brewer.pal(11, "Paired"))
