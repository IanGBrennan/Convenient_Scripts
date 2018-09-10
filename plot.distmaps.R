## This script will take an input data frame of distribution data
## column names "Name_in_Tree", "Latitude", "Longitude", and will:
# (1) convert lat/lon values to Spatial Points Objects (named by species)
# (2) create Convex Polygons of distribution for each species (named by species)
# (3) translate polygons into Observational Windows (OWin objects)
# (4) if provided with an output directory name, will plot the ranges


require(sp)
require(rworldmap)
require(ggmap)
require(rgeos)
require(adehabitatHR)
require(maptools)

plot.distmaps <- function(distribution.table, new.directory=NULL) {
  ## Next we need to determine (pairwise) if taxa overlap in their ranges
  ## this step only needs to be done once! (won't change with changes to the tree)
  all.taxa <- unique(distribution.table$Name_in_Tree)
  all.sp <- list()
  all.hull <- list()
  all.owin <- list()
  
  # create an ouput folder to catch the plots
  if (!is.null(new.directory)) {
    curr <- getwd()
    dir.path <- paste0(curr,"/", new.directory)
    call <- paste("mkdir", dir.path)
    system(call)
    
    sbbox <- make_bbox(lon = distribution.table$Longitude, lat = distribution.table$Latitude, f = .1)
    sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
  }
  
  for (p in 1:length(unique(distribution.table$Name_in_Tree))) {
    current.taxon <- all.taxa[[p]]
    cat("iteration", p, "of", length(unique(distribution.table$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
    
    current.data <- subset(distribution.table, distribution.table$Name_in_Tree == current.taxon)
    if (nrow(current.data) > 1000) {
      current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
    } 
    points <- current.data[,c(2,3)]
    all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
    all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)
    all.owin[[p]] <- as.owin(all.hull[[p]])
    
    # only plot the maps below if you need to check/clean the data
    if(!is.null(new.directory)) {
      fortified.data <- fortify(distribution.hull)
      pdf(paste(dir.path, "/", current.taxon, ".pdf", sep=""))
      print(ggmap(sq_map) 
            + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
            + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
      dev.off()
    }
  }
  names(all.sp) <- unique(distribution.table$Name_in_Tree)
  names(all.hull) <- unique(distribution.table$Name_in_Tree)
  names(all.owin) <- unique(distribution.table$Name_in_Tree)
  
  return(list(SpatialPoints=all.sp, ConvexHulls=all.hull, OWin=all.owin))
}

#test <- plot.distmaps(distribution)
