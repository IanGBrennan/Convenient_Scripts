library(phytools); library(ALA4R); library(raster); library(rgeos); library(sp); library(vegan)
library(FD); library(raster); library(rgdal); library(ggmap); library(broom); library(dplyr)

source("~/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")

# remember to do the 'register_google(key="...", write=TRUE)' bit with the proper key from NANO


# Google Maps API key: AIzaSyBK4zmzavlL3Xh3223BqD22s8ui7cZ2rGw # this is old and won't work
#register_google(key = "AIzaSyBK4zmzavlL3Xh3223BqD22s8ui7cZ2rGw", write = TRUE) 

# Working through the Goannas and Marsupials
goanna.tree <- read.tree("~/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Goanna.tre")
goanna.dist <- read.csv("~/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Goanna_Distributions.csv", header=T)
marsupial.tree <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Marsupial.tre")
marsupial.dist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Marsupial_Distributions.csv", header=T)
joint.dist <- rbind(marsupial.dist, goanna.dist)

# plot the distributions of extant (tip) taxa
goanna.maps <- plot.distmaps(goanna.dist, new.directory = NULL, point.width = 0.5)
marsupial.maps <- plot.distmaps(marsupial.dist, new.directory = NULL, point.width = 0.5)
joint.maps <- c(goanna.maps$ConvexHulls, marsupial.maps$ConvexHulls)

# Read in the Trait Data we'll use
joint.trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/Goanna.Marsupial.DATA.csv", header=T)
goanna.trait <- filter(joint.trait, Name_in_Tree %in% goanna.tree$tip.label)


# ALA4R instructions: https://atlasoflivingaustralia.github.io/ALA4R/articles/ALA4R.html
shape.australia <- shapefile("/Users/Ian/Desktop/Map_Shapefiles/Australia.shx")
wkt.oz <- writeWKT(shape.australia)
#wkt.oz <- "POLYGON((152.5 -35,152.5 -32,140 -32,140 -35,152.5 -35))"
x <- occurrences(taxon="genus:Varanus", wkt=wkt.oz, qa="none", download_reason_id="testing")

# Create a tibble from the distribution data, turn it into Site X Species tibble
ygridded <- goanna.dist %>% # if you wanted to, you could change this to just the goannas and look at them instead (goanna.dist) or (joint.dist)
  ## discard genus- and higher-level records
#  dplyr::filter(rank %in%
#                  c("species", "subspecies", "variety", "form", "cultivar")) %>%
  
  ## bin into 0.5-degree bins
  dplyr::mutate(longitude=round(Longitude*2)/2, latitude=round(Latitude*2)/2) %>%
  
  #  ## average environmental vars within each bin
  group_by(longitude,latitude) %>%
  #  mutate(precipitationAnnual=mean(precipitationAnnual, na.rm=TRUE),
  #         temperatureAnnualMaxMean=mean(temperatureAnnualMaxMean, na.rm=TRUE)) %>%
  
  ## subset to vars of interest
  dplyr::select(longitude, latitude, Name_in_Tree) %>%
  
  ## take one row per cell per species (presence)
  distinct() %>%
  
  ## calculate species richness
  dplyr::mutate(richness=n()) %>%
  
  ## convert to wide format (sites by species)
  dplyr::mutate(present=1) %>%
  do(tidyr::spread(data=., key=Name_in_Tree, value=present, fill=0)) %>%
  ungroup()

# Plot diversity as a function of latitude
ggplot(ygridded, aes(latitude, richness)) + geom_point() + theme_bw()

# Change the tibble to a data frame we can work with normally
gridded.dist <- as.data.frame(ygridded)
gridded.dist[is.na(gridded.dist)] <- 0 # make NAs 0
gridded.dist <- filter(gridded.dist, !richness==1) # remove sites with just one taxon
gridded.dist <- filter(gridded.dist, latitude <= -11); gridded.dist <- filter(gridded.dist, longitude >= 113.5)
gdist <- gridded.dist[ , 4:ncol(gridded.dist)]


# make the order of the trait dataframe match the order of the Site X Species DF
joint.trait <- joint.trait[match(colnames(gdist), joint.trait$Name_in_Tree),]
goanna.trait <- goanna.trait[match(colnames(gdist), goanna.trait$Name_in_Tree),]
trait.frame <- as.data.frame(joint.trait$Body_Length); rownames(trait.frame) <- joint.trait$Name_in_Tree
goanna.frame <- as.data.frame(goanna.trait$Body_Length); rownames(goanna.frame) <- goanna.trait$Name_in_Tree

# Run the Functional Diversity function and extract just the Rao's Quadratic value
test <- dbFD(trait.frame, gdist)
RQ.scores <- test$RaoQ
res.table <- cbind.data.frame(latitude=gridded.dist$latitude, longitude=gridded.dist$longitude, RaoQ=test$RaoQ)

best <- dbFD(goanna.frame, gdist)
RQ.scores <- best$RaoQ
res.table <- cbind.data.frame(latitude=gridded.dist$latitude, longitude=gridded.dist$longitude, RaoQ=best$RaoQ)


## Read in your shapefile
oz <- shapefile("/Users/Ian/Desktop/Map_Shapefiles/Australia.shp")
plot(oz)

## Set up a raster "template" for a 0.5 degree grid
ext <- extent(113.2244, 153.6242, -43.64806, -10.70667)
gridsize <- 0.5
r <- raster(ext, res=gridsize)

## Rasterize the shapefile
rr <- rasterize(oz, r)

## Plot raster
plot(rr)

rr.cells <- xyFromCell(rr, 1:length(rr)); rr.cells <- as.data.frame(rr.cells)
rr.cells$x <- round(rr.cells$x*2)/2; rr.cells$y <- round(rr.cells$y*2)/2
colnames(rr.cells) <- c("longitude", "latitude")

#values(rr)

combo.Q <- left_join(rr.cells, res.table, by=c("longitude", "latitude"))
FDras <- rr
values(FDras) <- combo.Q$RaoQ
plot(FDras)

combo.R <- left_join(rr.cells, gridded.dist[,1:3], by=c("longitude", "latitude"))
RICHras <- rr
values(RICHras) <- combo.R$richness
plot(RICHras)

combo.SES <- left_join(rr.cells, ses.table, by=c("latitude", "longitude"))
SESras <- rr
values(SESras) <- combo.SES$SES
plot(SESras)


plot(oz, add=T)

hiFD <- res.table[which.max(res.table$RaoQ),]
filter(gridded.dist, latitude==hiFD$latitude & longitude==hiFD$longitude)

hiRICH <- gridded.dist[which.max(gridded.dist$richness),]


# my API key is in nano ~/.bash_profile
rangemap <- get_googlemap(center = "Australia", zoom = 4, style = 'feature:all|element:labels|visibility:off')
graymap <- get_googlemap(center = "Argentina", zoom = 4, style = 'https://maps.googleapis.com/maps/api/staticmap?&center=-33.9,151.14999999999998&zoom=12&format=png&maptype=roadmap&style=element:geometry%7Ccolor:0xf5f5f5&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360')
# designed here: https://mapstyle.withgoogle.com/
ggmap(rangemap); ggmap(graymap)

RICHpoly <- rasterToPolygons(RICHras); max.colors <- length(unique(RICHpoly$layer)); filled <- rep(RICHpoly$layer, each=5) # 'each' is important, otherwise the polygon values get screwed up
FDpoly <- rasterToPolygons(FDras); max.colors <- length(unique(FDpoly$layer)); filled <- rep(FDpoly$layer, each=5) # 'each' is important, otherwise the polygon values get screwed up
SESpoly <- rasterToPolygons(SESras); max.colors <- length(unique(SESpoly$layer)); filled <- rep(SESpoly$layer, each=5) # 'each' is important, otherwise the polygon values get screwed up

pal.length <- abs(min(ses.raster@data@values) - max(ses.raster@data@values)) * 10
myBreaks <- c(seq(min(ses.raster@data@values), 0, length.out=ceiling(pal.length/2) + 1), 
              seq(max(ses.raster@data@values)/pal.length, max(ses.raster@data@values), length.out=floor(pal.length/2)))

ggmap(graymap) + geom_polygon(data = SESpoly, 
                              aes(x = long, y = lat, group = group, fill = filled), size = 0, alpha = 1)  +
                              #scale_fill_gradientn("RasterValues", colors = wes_palette("Zissou1", max.colors, type="continuous")) + # use this for FDpoly and RICHpoly for the basic color ramp
                              #scale_fill_gradient2(low = "#F21A00", mid = "white", high = "#3B9AB2", midpoint = 0, breaks = myBreaks) +
                              scale_fill_gradientn(values=scales::rescale(c(min(myBreaks), -0.5, 0, 0.5, max(myBreaks))), colours=c("#2f7b8e","#78B7C5","#EBCC2A","#E1AF00","#F21A00")) + # Use this for SESpoly, I removed the second orange to fit my scale better ("#E1AF00")
                              #scale_fill_gradientn(values=scales::rescale(c(min(res.table$RaoQ), 0.4 , mean(res.table$RaoQ), 1, max(res.table$RaoQ))), colours=c("#2f7b8e","#78B7C5","#EBCC2A","#E1AF00","#F21A00")) + # Use this for FDpoly
                              #scale_fill_gradientn(values=scales::rescale(c(min(gridded.dist$richness), (mean(gridded.dist$richness) + min(gridded.dist$richness))/2 , mean(gridded.dist$richness), (mean(gridded.dist$richness) + max(gridded.dist$richness))/2, max(gridded.dist$richness))), colours=c("#2f7b8e","#78B7C5","#EBCC2A","#E1AF00","#F21A00")) + # Use this for RICHpoly
                              theme_classic()

#RICHmap, FDmap, SESmap
multiplot(RICHmap, FDmap, SESmap)


## Create a function to calculate the confidence interval of the SES
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error, 
              "error" = error, "mean" = vec_mean, "sd" = vec_sd, "N" = n)
  return(result)
}

# Can also be calculated as:
## upper = mean + (error * 1.96)
## lower = mean - (error * 1.96)


test <- confidence_interval(RQ$ses, 0.95)

data(finch.ind)

res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
                    sp = sp.finch, nperm = 9)
ses(res.finch$Tstats$T_IP.IC, res.finch$Tstats$T_IP.IC_nm)


ses(RQ$emp.val, RQ$sim.mean, val.quant=c(0.025, 0.975))

  
#

aus <- get_map("australia")

fort.oz <- fortify(rr)

# convert spatial object to a ggplot ready data frame
oz_df <- tidy(oz)
# make sure the shapefile attribute table has an id column
oz$id <- rownames(oz@data)
# join the attribute table from the spatial object to the new data frame
oz_df <- left_join(oz_df, oz@data, by = "id")


test_spdf <- as(rr, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")
ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8)  
  geom_polygon(data=oz, aes(x=long, y=lat), 
               fill=NA, color="grey50", size=0.25) +
  scale_fill_viridis() +
  coord_equal() +
  theme_classic() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"))


straya <- c(left=113.2244, bottom=-43.64806, right=153.6242,  top=-10.70667)
testo <- get_googlemap(straya)

sq_map <- get_googlemap(location = straya, maptype = "terrain", source = "google")
ggmap(sq_map)
ggplot(rr)

us <- c(left = -125, bottom = 25.75, right = -67, top = 49)
get_stamenmap(us, zoom = 5, maptype = "toner-lite") %>% ggmap() 
testo <- get_googlemap(center="Australia", zoom=4)
plot(testo)

ext <- extent(113.2244, 153.6242, -43.64806, -10.70667)

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



## where a species was not present, it will have NA: convert these to 0
sppcols <- setdiff(names(xgridded),
                   c("longitude", "latitude", "precipitationAnnual", "temperatureAnnualMaxMean",
                     "richness"))
xgridded <- xgridded %>% mutate_at(sppcols, function(z) ifelse(is.na(z), 0, z))

ss <- sites_by_species(taxon="rk_genus:Varanus",wkt="POLYGON((144 -43,148 -43,148 -40, 144 -40,144 -43))",gridsize=0.1,verbose=TRUE)




if (remove.extralimital==TRUE){
  spatial.frame <- SpatialPointsDataFrame(coords = dist.frame[, c("Longitude", "Latitude")], 
                                          data = dist.frame)
  input.shp <- shapefile(range.shape)
  inside.data <- spatial.frame[!is.na(over(spatial.frame, as(input.shp, "SpatialPolygons"))), ]
  dist.frame <- as.data.frame(inside.data[,1:3])
}

input.shp <- shapefile("/Users/Ian/Desktop/Australia.shx")
plot(input.shp)
plot(point1, add=T, col="red")


#### This was a function I created to get the FD from a number of localities, but it's not great
sample.FD <- function(traitDF, distMAP, sample.points, n.iter){
  FDres <- NULL
  
  randos <- sample(1:nrow(sample.points), n.iter)
  for(p in 1:n.iter){
    curr.point <- SpatialPoints(sample.points[randos[p],c(2,3)])
    overlaps <- unlist(lapply(distMAP, function(x) gIntersects(curr.point, x)))
    overlap.res <- as.data.frame(rownames=names(overlaps), overlaps)
    
    filter.overlaps <- overlap.res %>%
      rownames_to_column('taxon') %>%
      filter(overlap.res, overlaps==T) %>%
      column_to_rownames('taxon')
    
    filter.size <- filter(traitDF, Name_in_Tree %in% rownames(filter.overlaps))
    curr.locale <- as.data.frame(filter.size$Body_Length); rownames(curr.locale) <- filter.size$Name_in_Tree
    curr.FD <- dbFD(curr.locale)
    
    rando.size <- traitDF[sample(traitDF, nrow(filter.size))]
    joint.trait[sample(1:nrow(joint.trait), 3),]
    rando.size <- traitDF[sample(1:nrow(traitDF), nrow(filter.size))]
    
    
    curr.res <- cbind(sample.points[randos[p],c(2,3)], Rao.Q=curr.FD$RaoQ)
    curr.res <- cbind(curr.res, Taxa=paste(filter.size$Name_in_Tree, collapse=" "))
    FDres <- rbind(FDres, curr.res)
    
  }
  return(FDres)
}
testo <- sample.FD(joint.trait, joint.maps, goanna.dist, n.iter=10)

gIntersects(point1, goanna.maps$ConvexHulls$Varanus.acanthurus)
gIntersects(point1, goanna.maps$ConvexHulls$Varanus.komodoensis)

overlaps <- unlist(lapply(joint.maps, function(x) gIntersects(point1, x)))
overlap.res <- as.data.frame(rownames=names(overlaps), overlaps)

filter.overlaps <- overlap.res %>%
  rownames_to_column('taxon') %>%
  filter(overlap.res, overlaps==T) %>%
  column_to_rownames('taxon')

filter.size <- filter(joint.trait, Name_in_Tree %in% rownames(filter.overlaps))
curr.locale <- as.data.frame(filter.size$Body_Length); rownames(curr.locale) <- filter.size$Name_in_Tree
curr.FD <- dbFD(curr.locale)

plot(filter.size$Body_Length)


overlap.res[joint.trait$Name_in_Tree,]
overlap.res[order(joint.trait$Name_in_Tree),]

filter.size <- rbind(filter.size, joint.trait[32,])


####################################################################################
## Plot some community composition information (body size, but could be any trait)
####################################################################################
unique(gridded.dist$richness); pick.richness <- 14
a.cell <- gridded.dist[which(gridded.dist$richness == pick.richness),]; 
a.cell <- a.cell[sample(1:nrow(a.cell),1),]
a.cell <- a.cell[ , which((!a.cell[1,]==0))]; 
cell.taxa <- names(a.cell)[4:length(names(a.cell))]; cell.taxa

# Read in the Trait Data (get in RPANDA format)
joint.trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/Goanna.Marsupial.DATA.csv", header=T)
cell.traits <- filter(joint.trait, Name_in_Tree %in% cell.taxa); cell.traits <- cell.traits[order(cell.traits$Body_Length),]
cell.traits$Name_in_Tree <- factor(cell.traits$Name_in_Tree, levels=cell.traits$Name_in_Tree)
cell.traits$barcol <- colorRampPalette(brewer.pal(9, "YlGnBu"))(nrow(cell.traits)+1)[2:(nrow(cell.traits)+1)]

i14 <- ggplot(cell.traits, aes(x=Name_in_Tree, y=Body_Length)) +
  geom_bar(stat="identity", fill=colorRampPalette(brewer.pal(9, "YlGn"))(nrow(cell.traits))) +
  theme_classic() +
  geom_text(aes(label=Name_in_Tree, vjust=-.5)) +
  labs(title = paste(a.cell$latitude, a.cell$longitude))

gridExtra::grid.arrange(a3, a4, a5, a6, b7, c8, d9, e10, f11, g12, h13, i14, j15)

# Ideally we'd make a loop to just plot each cell that has richness > n



