library(rase)
library(coda)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")


# start by creating an ultrametric tree of known height and branching pattern
newick1 <- "((((A:1,B:1):3,(C:3,D:3):1):2,E:6):1,((X:1.5,Y:1.5):3,Z:4.5):2.5);"
tree_1 <- read.tree(text=newick1)
plot(tree_1)

# read in distribution data that matches our tip labels
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/GMM_TEST_DATA.csv", header=T)
    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]

# plot the distributions of extant (tip) taxa
tips <- plot.distmaps(distribution)

# sort the data to make sure the order matches the tree appropriately
tree_poly <- name.poly(tips$OWin, tree_1, poly.names = unique(distribution$Name_in_Tree))

# run the MCMC
res <- rase(tree_1, tree_poly, niter=1000, logevery = 10)
# extract and plot the MCMC output
resmc <- mcmc(res, start=(length(res[,1])*.2)) # remove 20% as burnin
par(mar=c(1,1,1,1))
plot(resmc)

# here's a function for extracting the node distributions
# and turning it into Spatial Data that we can add to our existing records
process.rase <- function(mcmc.object, new.directory) {
  dist.frame <- NULL
  res.points <- as.data.frame(mcmc.object[,c(1:(ncol(mcmc.object)-2))])
  for (p in 1:(ncol(res.points)/2)) {
    int <- matrix(nrow=nrow(res.points), ncol=3)
    int <- as.data.frame(int)
    node <- test <- strsplit(colnames(res.points)[p], "_")
    int[,1] <- node[[1]][1]
    int[,2] <- res.points[,p]
    int[,3] <- res.points[,(p+(ncol(res.points)/2))]
    dist.frame <- rbind(dist.frame, int)
  }
  colnames(dist.frame) <- c("Name_in_Tree", "Latitude", "Longitude")
  
  # create an ouput folder to catch the plots
  curr <- getwd()
  dir.path <- paste0(curr,"/", new.directory)
  call <- paste("mkdir", dir.path)
  system(call)
  
  all.taxa <- unique(dist.frame$Name_in_Tree)
  all.sp <- list()
  all.hull <- list()
  all.owin <- list()
  sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
  sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
  
  for (p in 1:length(unique(dist.frame$Name_in_Tree))) {
    cat("iteration", p, "of", length(unique(dist.frame$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
    current.taxon <- all.taxa[[p]]
    current.data <- subset(dist.frame, dist.frame$Name_in_Tree == current.taxon)
    if (nrow(current.data) > 1000) {
      current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
    } 
    points <- current.data[,c(2,3)]
    all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
    all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)
    all.owin[[p]] <- as.owin(all.hull[[p]])
    
    # only plot the maps below if you need to check/clean the data
    fortified.data <- fortify(distribution.hull)
    pdf(paste(dir.path,"/", current.taxon, ".pdf", sep=""))
    print(ggmap(sq_map) 
      + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
      + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
    dev.off()
  }
  names(all.sp) <- unique(dist.frame$Name_in_Tree)
  names(all.hull) <- unique(dist.frame$Name_in_Tree)
  
  plot.location.call <- paste("plots printed to directory:", dir.path)
  print(plot.location.call)
  
  return(list(DistData=dist.frame, SpatialPoints=all.sp, ConvexHulls=all.hull, OWin=all.owin))
}

rase.out <- process.rase(resmc, "doodies")

remove.extralimital <- 




df3 <- data.for.3d(resmc, tree_1, tree_poly)
phylo.3d(df3, z.scale=10, pts=T)
add.polygons(df3)
add.dens(df3, resmc, z.scale=10, col=c(2:8))


test1 <- as.owin(all.hull[[1]])
plot(test1)
test1

shape <- readOGR(dsn = "/Users/Ian/Desktop/Australia.shp", layer = "SHAPEFILE")

library(raster)
oz <- shapefile("/Users/Ian/Desktop/Australia.dbf")
plot(oz)

oz = readOGR(dsn=".", layer="Australia")
oz@data$id = rownames(oz@data)
oz.points = fortify(oz, region="id")
oz.df = join(oz.points, oz@data, by="id")

plot(oz.df)

(ggplot(oz.df) + 
  aes(long,lat,group=group) + 
  geom_polygon() +
  geom_path(color="white") +
  coord_equal() +
  theme_bw())


all.hull[[1]]

taxon1 <- filter(distribution, Name_in_Tree == "A")
taxon1.pts <- taxon1[,c(3,2)] # change order to Lon/Lat
spoints <- SpatialPoints(taxon1.pts)
test <- !is.na(over(spoints, as(oz, "SpatialPolygons")))
mean(test)
over(spoints, as(oz, "SpatialPolygons"))
oz.poly <- as(oz, "SpatialPolygons")
over(spoints, oz.poly)


all.X <- list();
all.hull.X <- list()

# make a loop that creates distributions of all tips FOR THE FIRST TREE
for (p in 1:Ntip(tree_1)) {
  current.taxon <- tree_1$tip.label[[p]]
  current.data <- subset(distribution, distribution[,1] == current.taxon)
  if (nrow(current.data) > 1000) {
    current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
  } 
  points <- current.data[,c(2,3)]
  all.X[[p]] <- distribution.sp <- SpatialPoints(points)
  all.hull.X[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=1)
}
names(all.X) <- tree_1$tip.label
names(all.hull.X) <- tree_1$tip.label
#

test <- paste("mkdir", "/Users/Ian/Desktop/Butt_Maps")
system(test)
