library(phytools)


# Read in your trees
trees <- read.nexus("/Users/Ian/Desktop/BP.fixed.trees")

# Loop through the trees, and select trees that fit a certain age criteria
good.trees <- list()
for (i in 1:length(trees)) {
  above.root <- findMRCA(trees[[i]], tips=c("Anilius_scytale", "Tropidophis_pardalis"), type="height") #gives height ABOVE the ROOT
  split.time <- max(nodeHeights(trees[[i]])) - above.root # so subtract that from the total depth!
  if (split.time < 85 && split.time > 40) { #designate the age range you're looking for
    good.trees[[length(good.trees)+1]] <- trees[[i]]
  }
}
class(good.trees) <- "multiPhylo" #set the class of the good trees
#write.tree(good.trees, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Boas.Pythons.100.trees")

# in case you need to thin out the set of trees further
treex <- sample(good.trees, size=100) # designate how many you want
write.nexus(treex, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Aniliidae.Tropidophiidae.100.trees")

###########################
keep <- c("Anilius_scytale",
          "Trachyboa_boulengeri",
          "Trachyboa_gularis",
          "Tropidophis_curtus_curtus",
          "Tropidophis_feicki",
          "Tropidophis_greenwayi",
          "Tropidophis_haetianus_haetianus",
          "Tropidophis_melanurus_melanurus",
          "Tropidophis_pardalis",
          "Tropidophis_taczanowskyi",
          "Tropidophis_wrighti")
all <- treex[[1]]$tip.label
drop <- setdiff(all, keep)
treem <- lapply(treex, drop.tip, tip=drop)
less <- sample(treem, size=100)
write.nexus(less, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Aniliidae.Tropidophiidae.100.trees")
