library(OUwie)
library(phytools)
#### For 'OUwie' we'll need to set up our data according to the question we're asking
#########################################################################
make.OUwie.input <- function(data, regime, level.chars, taxa, phy, trait) {
  regime.data <- regime
  regime.data <- as.character(regime.data)
  regime.legend <- list()
  regime.names <- unique(regime.data)
  for (i in 1:length(unique(regime.data))) {
    regime.legend[i] <- paste(regime.names[i], "=", level.chars[i])
    regime.data[regime.data == regime.names[i]] <- level.chars[i]
  }
  names(regime.data) <- taxa
  regime.simm <- make.simmap(phy, regime.data)
  plotSimmap(regime.simm, fsize=0)
  regime.trait <- as.data.frame(regime.data)
  regime.trait <- cbind(regime.trait, as.data.frame(trait))
  combined <- NULL
  combined <- as.data.frame(taxa)
  combined <- cbind.data.frame(combined, regime.trait)
  
  all.result <- list()
  all.result[[1]] <- regime.data
  all.result[[2]] <- regime.trait
  all.result[[3]] <- combined
  all.result[[4]] <- regime.simm
  all.result[[5]] <- regime.legend
  
  return(all.result)
}
#########################################################################

#########################################################################
# the function above helps to create input files for OUwie
# data = your whole data frame;
# regime = data$regime; 'regime' is a column with data as a word or phrase
# unique(data$regime); this will give you the order the regime vector is in
# level.chars = numeric matrix corresponding to the # of regimes; (0,1,2) or (0,1,0), etc
# taxa = data$taxon.names; ('tree_id' or 'name.in.tree' column)
# tree = phylo; the tree you're basing it off of
# trait = data[,"..."];

# the output values of the function are:
# output[[1]] is the regime matched to taxon names
# output[[2]] is a data frame of the regimes and trait values matched to the taxon names
# output[[3]] is the OUwie formatted data [taxon, regime, trait]
# output[[4]] is the Simmap object
# output[[5]] is a legend for the regime states

# usage: make.OUwie.input(morph.data, morph.data$Phylogeny, c(0,1,2,1,1,3,4,0,2,2,0,2,3), morph.data$Name.in.Tree, tree, morph.data$HL.Trunk)