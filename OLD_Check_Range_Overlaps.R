######################
## This was my old way of addressing species ranges and overlap:
  ## I used LoCoH to create convex hulls 
######################

## Next we need to determine (pairwise) if taxa overlap in their ranges
## this step only needs to be done once! (won't change with changes to the tree)
all.taxa <- unique(distribution$Name_in_Tree)
all.sp <- list()
all.hull <- list()
#all.poly <- list()
sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")

for (p in 202:length(unique(distribution$Name_in_Tree))) {
  current.taxon <- all.taxa[[p]]
  current.data <- subset(distribution, distribution$Name_in_Tree == current.taxon)
  if (nrow(current.data) > 1000) {
    current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
  } 
  points <- current.data[,c(2,3)]
  all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
  distribution.hull <- LoCoH.k(distribution.sp, k=length(distribution.sp)/2, duplicates="random")
  all.hull[[p]] <- gUnionCascaded(distribution.hull)
  #all.poly[[p]] <- distribution.poly <- mcp(distribution.sp, percent=100)
  
  # only plot the maps below if you need to check/clean the data
  pdf(paste("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Skink_Maps/", current.taxon, ".pdf", sep=""))
  print(ggmap(sq_map) + geom_point(data = current.data, mapping = aes(x = Longitude, y = Latitude), color="red"))
  dev.off()
}
names(all.sp) <- unique(distribution$Name_in_Tree)
names(all.hull) <- unique(distribution$Name_in_Tree)
saveRDS(all.sp,   file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.SpatialPoints.RDS")
saveRDS(all.hull, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.Polygons.RDS")

range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
arrange <- data.frame(species1 = unique.combo[1,], 
                      species2 = unique.combo[2,], 
                      range_overlap, stringsAsFactors=F)
for (t in 1:length(arrange[,1])) {
  overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
  arrange[t,3] <- overlap
}

all.TTR <- lapply(all.TTR, cbind, arrange$range_overlap)
all.TTR <- lapply(all.TTR, setNames, c("species1", "species2", "dist_trait", "dist_tree", "range_overlap"))
saveRDS(all.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.Matrices.RDS")