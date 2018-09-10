library(rase)
library(coda)
library(ggmap)
library(raster)
library(purrr); library(magick); install.packages("ImageMagick")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")


# start by creating an ultrametric tree of known height and branching pattern
newick1 <- "((((A:1,B:1):3,(C:3,D:3):1):2,E:6):1,((X:1.5,Y:1.5):3,Z:4.5):2.5);"
tree_1 <- read.tree(text=newick1)
plot(tree_1)

# read in distribution data that matches our tip labels
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/GMM_TEST_DATA.csv", header=T)
    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]

# plot the distributions of extant (tip) taxa
tips <- plot.distmaps(aprasia.dist, new.directory = "TESTPLOTMAPS")

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
process.rase <- function(mcmc.object, distribution, new.directory=NULL,
                         remove.extralimital=NULL, range.shape) {
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
  
  if (remove.extralimital==TRUE){
    spatial.frame <- SpatialPointsDataFrame(coords = dist.frame[, c("Longitude", "Latitude")], 
                                            data = dist.frame)
    input.shp <- shapefile(range.shape)
    inside.data <- spatial.frame[!is.na(over(spatial.frame, as(input.shp, "SpatialPolygons"))), ]
    dist.frame <- as.data.frame(inside.data[,1:3])
  }
  
  # create an ouput folder to catch the plots
  if(!is.null(new.directory)) {
    curr <- getwd()
    dir.path <- paste0(curr,"/", new.directory)
    call <- paste("mkdir", dir.path)
    system(call)
  }
  
  all.taxa <- unique(dist.frame$Name_in_Tree)
  all.sp <- list()
  all.hull <- list()
  all.owin <- list()
  
  if(!is.null(new.directory)) {
    sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
    sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
  }
  
  for (p in 1:length(unique(dist.frame$Name_in_Tree))) {
    cat("iteration", p, "of", length(unique(dist.frame$Name_in_Tree)), "\n") #keep track of what tree/loop# we're on
    current.taxon <- all.taxa[[p]]
    current.data <- subset(dist.frame, dist.frame$Name_in_Tree == current.taxon)
    if (nrow(current.data) > 1000) {
      current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
    } 
    points <- current.data[,c(2,3)]
    all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
    all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=0.5)
    all.owin[[p]] <- as.owin(all.hull[[p]])
    
    # only plot the maps below if you need to check/clean the data
    if(!is.null(new.directory)) {
    fortified.data <- fortify(distribution.hull)
    pdf(paste(dir.path,"/", current.taxon, ".pdf", sep=""))
    print(ggmap(sq_map) 
      + geom_polygon(data = fortified.data, aes(x = lat, y = long, group=group, fill="tomato"), alpha=0.3)
      + geom_point(data = current.data, mapping = aes(x=Longitude, y=Latitude), color="tomato3"))
    dev.off()
    }
  }
  names(all.sp) <- unique(dist.frame$Name_in_Tree)
  names(all.hull) <- unique(dist.frame$Name_in_Tree)
  
  if(!is.null(new.directory)) {
    plot.location.call <- paste("plots printed to directory:", dir.path)
    print(plot.location.call)
  }
  
  return(list(DistData=dist.frame, SpatialPoints=all.sp, ConvexHulls=all.hull, OWin=all.owin))
}
# mcmc.object: the output from 'rase' (after burnin via 'mcmc')
# distribution: the lat/long data frame used as an input
# new directory: if you want plots of the nodes, give a name for the new directory
# remove.extralimital: would you like to drop points outside of a boundary?
# range.shape: if 'yes' above, provide the path to a .shx shape file (e.g."/Users/Ian/Desktop/Australia.shx")
rase.out <- process.rase(mcmc.object=resmc, distribution = distribution, new.directory="TEST",
                          remove.extralimital = T, "/Users/Ian/Desktop/Australia.shx")

# plot the rase results as a 3D tree range-polygon thingy
df3 <- data.for.3d(resmc, tree_1, tree_poly)
phylo.3d(df3, z.scale=10, pts=T)
add.polygons(df3)
add.dens(df3, resmc, z.scale=10, col=c(2:8))


# Empirical Example with Aprasia
aprasia <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.tre")
pygo.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Aprasia.csv", header=T)
pygo.dist <- pygo.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
### trim down distributional data to match the tree
aprasia.dist <- filter(pygo.dist, Name_in_Tree %in% aprasia$tip.label)

# plot the distributions of extant (tip) taxa
tips <- plot.distmaps(lerista.dist, new.directory = NULL)

# sort the data to make sure the order matches the tree appropriately
tree_poly <- name.poly(tips$OWin, lerista, poly.names = unique(lerista.dist$Name_in_Tree))

# run the MCMC
res <- rase(lerista, tree_poly, niter=1000, logevery = 10)
# extract and plot the MCMC output
resmc <- mcmc(res, start=(length(res[,1])*.2)) # remove 20% as burnin
par(mar=c(1,1,1,1))
plot(resmc)

# process the rase-output
rase.out <- process.rase(mcmc.object=resmc, new.directory="/Desktop/Lerista_TEST",
                         distribution = lerista.dist, remove.extralimital = T, "/Users/Ian/Desktop/Australia.shx")
saveRDS(rase.out, file="/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Lerista.RASE.RDS")

# plot the rase results as a 3D tree range-polygon thingy
df3 <- data.for.3d(resmc, lerista, tree_poly)
phylo.3d(df3, z.scale=10, pts=T)
add.polygons(df3)
add.dens(df3, resmc, z.scale=10, col=c(2:8))






test1 <- as.owin(all.hull[[1]])
plot(test1)
test1

oz.shape <- readOGR(dsn = "/Users/Ian/Desktop/Australia.shp", layer = "SHAPEFILE")

library(raster); library(rgdal)
oz <- shapefile("/Users/Ian/Desktop/Australia.shx")
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


# geo.object: is a result of "CreateGeoObject" or "CreateGeoObject_SP"
    # e.g.: lerista.geo <- CreateGeoObject_SP(lerista, lerista.dist)
# for some reason this function isn't working properly yet:
    # so don't use the function as it won't plot the node ranges after a certain point
    # but you can just work through the parts of the function by yourself and it's ok
tree.through.time <- function(rase.input, geo.object, phy, map, 
                              country="Australia", new.directory=NULL){
  plot.times <- geo.object$times
  #taxa.times <- NULL
  
  #rangemap <- get_map(country, zoom=4, maptype="terrain",source="google")
  # rangemap <- get_googlemap(center = "Australia", zoom = 4,
  #                           style = 'feature:all|element:labels|visibility:off')
  rangemap <- australia
  
  # create an ouput folder to catch the plots
  if (!is.null(new.directory)) {
    #curr <- getwd()
    dir.path <- paste0(getwd(),"/", new.directory)
    call <- paste("mkdir", dir.path)
    system(call)
  }
  for (k in 1:length(geo.object$geography.object)) {
    taxa.to.plot <- NULL
    extant.names <- rownames(geo.object$geography.object[[k]])
    for (j in 1:length(extant.names)) {
      if(extant.names[[j]] %in% phy$tip.label){
        #current.taxon <- filter(map, Name_in_Tree == extant.names[[j]])
        taxa.to.plot <- rbind(taxa.to.plot, filter(map, Name_in_Tree == extant.names[[j]]))
      } 
      else if(!extant.names[[j]] %in% phy$tip.label) {
        node.name <- paste0("n", geo.object$name.matrix[which(geo.object$name.matrix[,2]==extant.names[[j]]),3])
        #node.name <- paste0("n",node.name)
        taxa.to.plot <- rbind(taxa.to.plot, filter(rase.input$DistData, Name_in_Tree == node.name)[,1:3])
      } else {
        stop("current taxon does not have distributional data")
      }
    }
    #current.time <- paste0(round((max(nodeHeights(phy))-plot.times[[k+1]]),digits=1)," MYA")
    # remember R will go "1,10,11,...,2,20,21.." so rename anything less than 10 as 01,02...
    pdf(paste(dir.path, "/", if(k<10){paste0("0",k)}else{k}, ".pdf", sep=""))
    print(ggmap(rangemap) + 
      stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                     geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none') + annotate("text", x=120, y=-40, label= paste0(round(max(plot.times) - plot.times[[k+1]],digits=1)," MYA"), color="white", cex=10))
    dev.off()
  }
  # Step 2: List those Plots, Read them in, and then make animation
  list.files(path = paste0(getwd(),"/",new.directory,"/"), pattern = "*.pdf", full.names = T) %>% 
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=1) %>% # animates, can opt for number of loops
    image_write(paste0(getwd(),"/",new.directory,"/",new.directory,".gif")) # write to current dir
}

# don't use this at the moment, because it's not working
tree.through.time(lerista.rase, lerista.geo, lerista.dist, lerista, 
                          country="Australia", new.directory = "Lerista_Ancestors")

###################
## Below is just practices stuff for the mapping function above
###################

# define the area of the map, call it from 'ggmap'
#australia <- get_map("Australia", zoom=4, maptype="watercolor",source="stamen", color="bw")
australia <- get_googlemap(center = 'Australia', zoom = 4,
              style = 'feature:all|element:labels|visibility:off')
# plot the distributions as density ranges
ggmap(australia) + 
  stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                 geom = "polygon", data = lerista.rase$DistData) + theme(legend.position = 'none') + annotate("text", x=8, y=13000, label= "ALL TIME")

ggmap(australia) + 
  stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                 geom = "polygon", data = filter(lerista.dist))

for (h in 68:73){
#for (h in 1:length(unique(lerista.dist$Name_in_Tree))){
  print(paste(h, unique(lerista.dist$Name_in_Tree)[h]))
  print(ggmap(australia) + 
    stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                   geom = "polygon", data = filter(lerista.dist, Name_in_Tree == unique(lerista.dist$Name_in_Tree)[h])))
}

# Create a spatstat ppp object
oz.owin <- as.owin.SpatialPolygons(oz)
pts <- as.ppp(aprasia.dist[,3:2], W=oz.owin)
plot(pts)
# Plot a a density surface
plot(density(pts))
