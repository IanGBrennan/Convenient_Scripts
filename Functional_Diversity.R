library(phytools); library(ALA4R); library(raster); library(rgeos); library(sp); library(vegan)
library(FD)

# Working through the Goannas and Marsupials
goanna.tree <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Goanna.tre")
goanna.dist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Goanna_Distributions.csv", header=T)
marsupial.tree <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Marsupial.tre")
marsupial.dist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/CoEvo_Marsupial_Distributions.csv", header=T)
joint.dist <- rbind(marsupial.dist, goanna.dist)

# plot the distributions of extant (tip) taxa
goanna.maps <- plot.distmaps(goanna.dist, new.directory = NULL, point.width = 0.5)
marsupial.maps <- plot.distmaps(marsupial.dist, new.directory = NULL, point.width = 0.5)
joint.maps <- c(goanna.maps$ConvexHulls, marsupial.maps$ConvexHulls)


point1 <- SpatialPoints(goanna.dist[2,c(2,3)])

gIntersects(point1, goanna.maps$ConvexHulls$Varanus.acanthurus)
gIntersects(point1, goanna.maps$ConvexHulls$Varanus.komodoensis)

overlaps <- unlist(lapply(joint.maps, function(x) gIntersects(point1, x)))
overlap.res <- as.data.frame(rownames=names(overlaps), overlaps)

joint.trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/Goanna.Marsupial.DATA.csv", header=T)
trait.frame <- as.data.frame(joint.trait$Body_Length); rownames(joint.trait$Name_in_Tree)

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

shape.australia <- shapefile("/Users/Ian/Desktop/Australia.shx")
wkt.oz <- writeWKT(shape.australia)
wkt.oz <- "POLYGON((152.5 -35,152.5 -32,140 -32,140 -35,152.5 -35))"
x <- occurrences(taxon="genus:Varanus", wkt=wkt.oz, qa="none", download_reason_id="testing")

xgridded <- x$data %>%
  ## discard genus- and higher-level records
  dplyr::filter(rank %in%
                  c("species", "subspecies", "variety", "form", "cultivar")) %>%
  
  ## bin into 0.5-degree bins
  mutate(longitude=round(longitude*2)/2, latitude=round(latitude*2)/2) %>%
  
#  ## average environmental vars within each bin
  group_by(longitude,latitude) %>%
#  mutate(precipitationAnnual=mean(precipitationAnnual, na.rm=TRUE),
#         temperatureAnnualMaxMean=mean(temperatureAnnualMaxMean, na.rm=TRUE)) %>%
  
  ## subset to vars of interest
  dplyr::select(longitude, latitude, scientificName) %>%
  
  ## take one row per cell per species (presence)
  distinct() %>%
  
  ## calculate species richness
  mutate(richness=n()) %>%
  
  ## convert to wide format (sites by species)
  mutate(present=1) %>%
  do(tidyr::spread(data=., key=scientificName, value=present, fill=0)) %>%
  ungroup()

xgridded <- as.data.frame(xgridded)
xgridded[is.na(xgridded)] <- 0
xgridded <- as_tibble(xgridded)

library(ggplot2)
ggplot(xgridded, aes(longitude, richness)) + geom_point() + theme_bw()



ygridded <- joint.dist %>%
  ## discard genus- and higher-level records
#  dplyr::filter(rank %in%
#                  c("species", "subspecies", "variety", "form", "cultivar")) %>%
  
  ## bin into 0.5-degree bins
  mutate(longitude=round(Longitude*2)/2, latitude=round(Latitude*2)/2) %>%
  
  #  ## average environmental vars within each bin
  group_by(longitude,latitude) %>%
  #  mutate(precipitationAnnual=mean(precipitationAnnual, na.rm=TRUE),
  #         temperatureAnnualMaxMean=mean(temperatureAnnualMaxMean, na.rm=TRUE)) %>%
  
  ## subset to vars of interest
  dplyr::select(longitude, latitude, Name_in_Tree) %>%
  
  ## take one row per cell per species (presence)
  distinct() %>%
  
  ## calculate species richness
  mutate(richness=n()) %>%
  
  ## convert to wide format (sites by species)
  mutate(present=1) %>%
  do(tidyr::spread(data=., key=Name_in_Tree, value=present, fill=0)) %>%
  ungroup()

ygridded <- as.data.frame(ygridded)
ygridded[is.na(ygridded)] <- 0
ygridded <- as_tibble(ygridded)

ggplot(ygridded, aes(longitude, richness)) + geom_point() + theme_bw()







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
