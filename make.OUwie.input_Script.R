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