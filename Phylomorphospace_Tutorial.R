library(phytools)
library(dplyr)
library(devtools)
# devtools::install_github("awhstin/awtools") # only need to do this once
library(awtools)

# start by reading in the tree of interest 
# trees are usually in one of two formats:
# if it's 'newick', use 'read.tree'
# if it's 'nexus', use 'read.nexus'
all.tree <- read.tree("/Users/Ian/Desktop/Olori.Pygopodoidea.tre")

# now read in the data you'd like to use
all.data <- read.csv("/Users/Ian/Desktop/Olori.Pygopodidae.Data.csv", header=T)

# You might notice the tree is for all the pygopodoid geckos 
# (Carphodatylidae, Diplodactylidae, Pygopodidae), but the data 
# is just for the Pygopodidae. We'll trim the data and tree down to match one another.
keepers <- intersect(all.tree$tip.label, all.data$Name_in_Tree)

# start by using 'drop.tip' and 'setdiff' to trim the tree
tree <- drop.tip(all.tree, setdiff(all.tree$tip.label, keepers))
# and then use 'filter' from 'dplyr' to reduce the data frame, 
# '%in%' searches for items of one matrix in another
pdata <- filter(all.data, Name_in_Tree %in% keepers)

# log transform the data to remove effects of massive differences in size, and normality
pdata[,5:7] <- log(pdata[,5:7]) # pick the traits of interest, here you'd use your first two PCs
rownames(pdata) <- pdata$Name_in_Tree # provide names for the data

morpho <- pdata[,c("SVL", "Head_Length")]
genus <- pdata$Genus; names(genus) <- pdata$Name_in_Tree

mycol <- character(length(genus))
mycol[genus=="Aprasia"] <- mpalette[2]
mycol[genus=="Delma"] <- mpalette[1]
mycol[genus=="Lialis"] <- mpalette[3]
mycol[genus=="Ophidiocephalus"] <- mpalette[4]
mycol[genus=="Paradelma"] <- mpalette[5]
mycol[genus=="Pygopus"] <- mpalette[5]
mycol[genus=="Pletholax"] <- mpalette[6]

pdata$color <- mycol
pdata <- pdata[match(tree$tip.label, pdata$Name_in_Tree),]
group.colors <- pdata$color
names(group.colors) <- 1:Ntip(tree)
nodecols <- rep("black", tree$Nnode)
names(nodecols) <- (Ntip(tree)+1) : (Ntip(tree)+tree$Nnode)
colorz <- c(group.colors, nodecols)

phylomorphospace(tree, morpho, 
                 label = "horizontal",
                 node.size = c(1, 3),
                 control=list(col.node=colorz),
                 xlab = "SVL", ylab = "Tail Length")




### Let's repeat the process coloring our tips by their ecology

#use 'filter' from 'dplyr' to reduce the data frame, '%in%' searches for items of one matrix in another
pdata <- filter(all.data, Name_in_Tree %in% keepers)

# log transform the data to remove effects of massive differences in size, and normality
pdata[,7:8] <- log(pdata[,7:8])
rownames(pdata) <- pdata$Name_in_Tree

morpho <- pdata[,c("SVL", "Tail_Length")]
ecology <- pdata$Ecology; names(ecology) <- pdata$Name_in_Tree

mycol <- character(length(ecology))
mycol[ecology=="arboreal"] <- mpalette[2]
mycol[ecology=="terrestrial"] <- mpalette[1]
mycol[ecology=="fossorial"] <- mpalette[3]

pdata$color <- mycol
pdata <- pdata[match(tree$tip.label, pdata$Name_in_Tree),]
group.colors <- pdata$color
names(group.colors) <- 1:Ntip(tree)
nodecols <- rep("black", tree$Nnode)
names(nodecols) <- (Ntip(tree)+1) : (Ntip(tree)+tree$Nnode)
colorz <- c(group.colors, nodecols)

phylomorphospace(tree, morpho, 
                 label = "horizontal",
                 node.size = c(1, 3),
                 control=list(col.node=colorz),
                 xlab = "SVL", ylab = "Tail Length")
