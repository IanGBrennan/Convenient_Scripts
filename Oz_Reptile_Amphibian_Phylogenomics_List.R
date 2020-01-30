library(phytools)
library(dplyr)
library(ggplot2)
library(RColorBrewer); library(scales)
library(gridExtra)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Extract_Data_from_AHE_labels.R")


python.tree <- read.nexus("/Users/Ian/Documents/ANU_Finished/T223_Pythonidae/FINALTREE_Python.tre")
goanna.tree <- read.nexus("/Users/Ian/Desktop/Species_Tree/FINAL_ExtantSampling/Operators/10Loci/UCLN_GTR/Extant_UCLN_GTR_10_con.tre")
    goanna.tree$tip.label
egernia.sampling <- read.csv("/Users/Ian/Google.Drive/ANU/AHE/T545_Egernia/Egernia_SampleInfo.csv", header=T)
    unique(egernia.sampling$full_genus_species)

sampling <- read.csv("~/Desktop/OzPhylogenomics_Sampling.csv", header=T)
    sampling <- select(sampling, -c("Member_Genera"))

sampling <- mutate(sampling, Percent_Species=Sampled_Species/Described_Species,
                   Percent_Genera=Sampled_Genera/Described_Genera)        
sampling

head(sampling)
samp_sort <- arrange(sampling, Group, Family)
samp_sort$Family <- factor(samp_sort$Family, levels = rev(samp_sort$Family))

#colorz <- colorRampPalette(brewer.pal(9, "RdYlBu"))
#samp_sort$Color <- colorz(23)
colorz <- c(rev(brewer.pal(nrow(filter(samp_sort, Group=="Anura")),"Reds")),
            #brewer.pal(nrow(filter(samp_sort, Group=="Crocodilia")),"Oranges"),
            "#FEE090",
            #rev(colorRampPalette(brewer.pal(9,"Greens"))(nrow(filter(samp_sort, Group=="Squamata_Lizards")))),
            rev(brewer.pal(nrow(filter(samp_sort, Group=="Squamata_Lizards")),"Greens")),
            rev(brewer.pal(nrow(filter(samp_sort, Group=="Squamata_Serpentes")),"Blues")),
            rev(brewer.pal(nrow(filter(samp_sort, Group=="Testudines")),"Purples")))
samp_sort$Color <- colorz
samp_sort$Color <- factor(samp_sort$Color, levels = rev(samp_sort$Color))


genera <- ggplot(samp_sort, aes(y=Percent_Genera, x=Family)) + 
  geom_bar(position="dodge", stat="identity", fill=rev(samp_sort$Color)) + 
  geom_text(aes(label=Described_Genera), position=position_dodge(width=1), vjust=-0.25) +
  theme_classic() + coord_flip()


species <- ggplot(samp_sort, aes(y=Percent_Species, x=Family)) + 
  geom_bar(position="dodge", stat="identity", fill=rev(samp_sort$Color)) + 
  geom_text(aes(label=Described_Species), position=position_dodge(width=1), vjust=-0.25) +
  theme_classic() + coord_flip()

grid.arrange(genera, species, nrow=1)


sqcl.tree <- read.nexus("~/Downloads/Fenker SQcl Brazi-aus.tre")
sqcl.tips <- extract.labels.SqCL(sqcl.tree); sqcl.tips <- data.frame(sqcl.tips)
sqcl.tips$Data <- "SqCL"
sqcl.tips$Organizer <- "Rabosky_Singhal_Moritz"
sqcl.tips$Project <- "Brazil_Oz_Convergence"
write.csv(sqcl.tips, "~/Desktop/SqCL_Coverage.csv")

elap.tree <- read.tree("/Users/Ian/Desktop/GenomeStripper/Elapids/Existing_Alignments/combined_alignments/Elapidae_NEW_ASTRAL.tre")
elap.tips <- extract.labels.basic(elap.tree); elap.tips <- data.frame(elap.tips)
elap.tips$Data <- "AHE"
elap.tips$Organizer <- "Keogh_Donnellan"
elap.tips$Family <- "Elapidae"
elap.tips$Project <- "Elapid_Phylogeny"
write.csv(elap.tips, "~/Desktop/AHE_Elapidae.csv")

goanna.tree <- read.nexus("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Treefiles/ASTRAL_Komodo.tre")
goanna.tips <- extract.labels.basic(goanna.tree); goanna.tips <- data.frame(goanna.tips)
goanna.tips$Data <- "AHE"
goanna.tips$Organizer <- "Keogh_Donnellan"
goanna.tips$Family <- "Varanidae"
goanna.tips$Project <- "Varanus_Phylogeny"
write.csv(goanna.tips, "~/Desktop/AHE_Varanidae.csv")

gex.tree <- read.tree("/Users/Ian/Google.Drive/ANU/Pygopodidae_ST/Pygopodoidea_Skipwith2019.tre")
gex.tips <- extract.labels.basic(gex.tree); gex.tips <- data.frame(gex.tips)
gex.tips$Data <- "UCE"
gex.tips$Organizer <- "Skipwith_Oliver"
gex.tips$Family <- "Diplodactylidae"
gex.tips$Project <- "Pygopodoidea_Phylogeny"
write.csv(gex.tips, "~/Desktop/UCE_Pygopodoidea.csv")
