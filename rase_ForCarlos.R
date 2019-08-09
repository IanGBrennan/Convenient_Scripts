library(rase)
library(coda)
library(ggmap)
library(raster)
library(purrr); library(magick); install.packages("ImageMagick")
library(dplyr)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/process.rase.R")


# Working through the Goannas and Marsupials
agamid.tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.tre")
agamid.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Agamids.csv", header=T)
    agamid.dist <- agamid.dist[,c(1,4,5)]

# trim the tree and data to match
keep <- intersect(agamid.tree$tip.label, unique(agamid.dist$Name_in_Tree))
a.tree <- drop.tip(agamid.tree, setdiff(agamid.tree$tip, keep))
a.dist <- filter(agamid.dist, Name_in_Tree %in% keep)
  
# plot the distributions of extant (tip) taxa
tips <- plot.distmaps(a.dist, new.directory = "/Desktop/TEST_Agamid_Distributions")

# sort the data to make sure the order matches the tree appropriately
tree_poly <- name.poly(tips$OWin, a.tree, poly.names = unique(a.dist$Name_in_Tree))

# run the MCMC
res <- rase(a.tree, tree_poly, niter=1000, logevery = 10)

# extract and plot the MCMC output
resmc <- mcmc(res, start=(length(res[,1])*.2)) # remove 20% as burnin
par(mar=c(1,1,1,1))
plot(resmc)

# process the rase-output
rase.out <- process.rase(mcmc.object=resmc, new.directory="/Desktop/Agamid_Test_RASE",
                         distribution = a.dist, remove.extralimital = T, "/Users/Ian/Desktop/Map_ShapeFiles/Australia.shx")

# add in the information for the tips
rase.out$Tip_DistData <- a.dist
rase.out$Tip_SpatialPoints <- tips$SpatialPoints
rase.out$Tip_ConvexHulls <- tips$ConvexHulls
rase.out$Tip_OWin <- tips$OWin
  
saveRDS(rase.out, file="/Users/Ian/Desktop/Agamid_Test_RASE.RDS")
