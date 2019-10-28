## This script will take an input data frame of distribution data
## column names "Name_in_Tree", "Latitude", "Longitude", and will:
# (1) convert lat/lon values to Spatial Points Objects (named by species)
# (2) create Convex Polygons of distribution for each species (named by species)
# (3) translate polygons into Observational Windows (OWin objects)
# (4) if provided with an output directory name, will plot the ranges

# Google Maps API: https://developers.google.com/maps/documentation/maps-static/styling

# Gonna have to register an API:
# https://stackoverflow.com/questions/52565472/get-map-not-passing-the-api-key-http-status-was-403-forbidden
# Google Maps API key: AIzaSyBK4zmzavlL3Xh3223BqD22s8ui7cZ2rGw
#register_google(key = "AIzaSyBK4zmzavlL3Xh3223BqD22s8ui7cZ2rGw", write = TRUE)


require(sp)
require(rworldmap)
require(ggmap)
require(rgeos)
require(adehabitatHR)
require(maptools)
require(spatstat)

plot.distmaps <- function(distribution.table, new.directory=NULL, base.map=NULL, point.width=1) {
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
    
    if (is.null(base.map)) {
      sbbox <- make_bbox(lon = distribution.table$Longitude, lat = distribution.table$Latitude, f = .1)
      sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
    } else {sq_map <- base.map}
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
    all.hull[[p]] <- distribution.hull <- gBuffer(distribution.sp, width=point.width)
    all.owin[[p]] <- as.owin(all.hull[[p]])
    
    # only plot the maps below if you need to check/clean the data
    if(!is.null(new.directory)) {
      fortified.data <- fortify(distribution.hull)
      pdf(paste(dir.path, "/", current.taxon, ".pdf", sep=""), paper="a4")
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

ssp.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/WA_Gecko_Subspecies/Sample_Localities.csv", header=T)
# my API key is in nano ~/.bash_profile
rangemap <- get_googlemap(center = "Australia", zoom = 4, style = 'feature:all|element:labels|visibility:off')
WAmap    <- get_googlemap(center = c(120.5,-24.5), zoom = 5, style = 'feature:all|element:labels|visibility:off')
WAGmap   <- get_googlemap(center = c(120.5,-26.5), zoom = 5, style = 'https://maps.googleapis.com/maps/api/staticmap?&center=-33.9,151.14999999999998&zoom=12&format=png&maptype=roadmap&style=element:geometry%7Ccolor:0xf5f5f5&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360')
graymap  <- get_googlemap(center = "Australia", zoom = 4, style = 'https://maps.googleapis.com/maps/api/staticmap?&center=-33.9,151.14999999999998&zoom=12&format=png&maptype=roadmap&style=element:geometry%7Ccolor:0xf5f5f5&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360')

testo <- get_googlemap(center = c(120,-24.2), zoom = 5, style="https://maps.googleapis.com/maps/api/staticmap?center=-24.229810153347742,-242.6166746271749&zoom=5&format=png&maptype=roadmap&style=element:labels%7Cvisibility:off&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:transit%7Cvisibility:off&size=480x360")

# designed here: https://mapstyle.withgoogle.com/
ggmap(rangemap); ggmap(graymap)
ggmap(WAmap); ggmap(WAGmap)


col.pal <- brewer.pal(8, "Paired")
cols <- brewer.pal(4, "Paired"); names(cols) <- c("Pletholax gracilis gracilis", "Pletholax gracilis edelensis",
                                                  "Nephrurus wheeleri wheeleri", "Nephrurus wheeleri cinctus")
ggmap(WAmap) + geom_point(aes(x = Longitude, y = Latitude, fill=Scientific.Name), colour = "black", 
                          data = ssp.dist, size = 2, pch = 21) + theme(legend.position="bottom") + scale_fill_brewer(palette="Paired", "Taxon")

ssp.geno <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/WA_Gecko_Subspecies/Subspecies_Localities.csv", header=T)
ssp.geno <- dplyr::filter(ssp.geno, Geno. == "Yes")
pg.ssp <- filter(ssp.geno, GENUS=="Pletholax")
nw.ssp <- filter(ssp.geno, GENUS=="Nephrurus")

pletholax.tree <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/WA_Gecko_Subspecies/Pletholax/Trees/FINAL_Pletholax_Partitions_July23.tre")
    ptree <- drop.tip(pletholax.tree, c("Aprasia.inaurita", "Aprasia.parapulchella", "Delma.butleri", 
                                        "Delma.tincta", "Lialis.burtonis", "Ophidiocephalus.taeniatus",
                                        "Paradelma.orientalis", "Pygopus.lepidopodus", "Pygopus.nigriceps"))
dist.data <- pg.ssp[,c("Name_in_Tree", "DDLatitude", "DDLongitude")]
opoints <- dist.data[,2:3]
rownames(opoints) <- dist.data$Name_in_Tree; colnames(opoints) <- c("lat", "long") 

cols <- c(rep("red", 22), rep("green", 3)); names(cols) <- ptree$tip.label
obj <- phylo.to.map(ptree, opoints, rotate=TRUE, database="worldHires",
                    regions="Australia", plot=FALSE)
plot(obj, colors=cols, direction="rightwards", ftype="off",cex.points=c(0,1),pts=FALSE)

neph.tree <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/WA_Gecko_Subspecies/Nephrurus/Trees/FINAL_Nephrurus_Partitions_July23.tre")
ntree <- drop.tip(neph.tree, c("Carphodactylus.laevis", "Nephrurus.amyae", "Nephrurus.asper",
                                    "Nephrurus.sheai.A", "Nephrurus.deleani", "Nephrurus.laevissimus",
                                    "Nephrurus.vertebralis","Nephrurus.levis.levis", "Nephrurus.levis.occidentalis",
                                    "Nephrurus.stellatus", "Underwoodisaurus.milii","Uvidocolis.sphyrurus",
                                    "Orraya.occultus","Phyllurus.platurus","Saltuarius.swaini"))
dist.data <- nw.ssp[,c("Name_in_Tree", "DDLatitude", "DDLongitude")]
opoints <- dist.data[,2:3]
rownames(opoints) <- dist.data$Name_in_Tree; colnames(opoints) <- c("lat", "long") 

cols <- c(rep("red", 2), rep("green", 5)); names(cols) <- ntree$tip.label
obj <- phylo.to.map(ntree, opoints, rotate=TRUE, database="worldHires",
                    regions="Australia", plot=FALSE)
plot(obj, colors=cols, ftype="off",cex.points=c(0,1),pts=FALSE)



#

