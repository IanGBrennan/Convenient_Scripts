# geo.object: is a result of "CreateGeoObject" or "CreateGeoObject_SP"
# e.g.: lerista.geo <- CreateGeoObject_SP(lerista, lerista.dist)
# for some reason this function isn't working properly yet:
# so don't use the function as it won't plot the node ranges after a certain point
# but you can just work through the parts of the function by yourself and it's ok

rangemap <- get_googlemap(center = "Australia", zoom = 4, style = 'feature:all|element:labels|visibility:off')

register_google(key = "your_API_key") 

tree.through.time <- function(rase.input, geo.object, distribution.data, phy, map, 
                              country="Australia", new.directory=NULL){
  plot.times <- geo.object$times
  #taxa.times <- NULL
  
  if (is.null(map)){
    rangemap <- get_map(country, zoom=4, maptype="terrain",source="osm")
    # rangemap <- get_googlemap(center = "Australia", zoom = 4,
    #                           style = 'feature:all|element:labels|visibility:off')
    #rangemap <- australia 
  }

  
  # create an ouput folder to catch the plots
  if (!is.null(new.directory)) {
    #curr <- getwd()
    dir.path <- paste0(getwd(),"/", new.directory)
    call <- paste("mkdir", dir.path)
    system(call)
  }
  for (k in 1:length(geo.object$geography.object)) {
    current.time <- round(max(plot.times) - plot.times[[k+1]],digits=1)
    print(paste("plotting distributions at", current.time, "time units ago"))
    taxa.to.plot <- NULL
    extant.names <- rownames(geo.object$geography.object[[k]])
    for (j in 1:length(extant.names)) {
      if(extant.names[[j]] %in% phy$tip.label){
        #current.taxon <- filter(map, Name_in_Tree == extant.names[[j]])
        taxa.to.plot <- rbind(taxa.to.plot, filter(distribution.data, Name_in_Tree == extant.names[[j]]))
      } 
      else if(!extant.names[[j]] %in% phy$tip.label) {
        node.name <- paste0("n", geo.object$name.matrix[which(geo.object$name.matrix[,2]==extant.names[[j]]),3])
        #node.name <- paste0("n",node.name)
        taxa.to.plot <- rbind(taxa.to.plot, filter(rase.input$DistData, Name_in_Tree == node.name)[,1:3])
      } else {
        stop("current taxon does not have distributional data")
      }
    }
    
    curr.color <- subset(colorz, names(colorz) %in% taxa.to.plot$Name_in_Tree)
    #current.time <- paste0(round((max(nodeHeights(phy))-plot.times[[k+1]]),digits=1)," MYA")
    # remember R will go "1,10,11,...,2,20,21.." so rename anything less than 10 as 01,02...
    pdf(paste(dir.path, "/", if(k<10){paste0("0",k)}else{k}, ".pdf", sep=""))
    print(ggmap(rangemap) + 
            stat_density2d(aes(x = Longitude, y = Latitude, fill = curr.color), alpha = .3, bins=4,
                           geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none') + annotate("text", x=120, y=-40, label= paste0(current.time," MYA"), color="white", cex=10))
    dev.off()
  }
  # Step 2: List those Plots, Read them in, and then make animation
  list.files(path = paste0(getwd(),"/",new.directory), pattern = "*.pdf", full.names = T) %>% 
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps=1) %>% # animates, can opt for number of loops
    image_write(paste0(getwd(),"/",new.directory,"/",new.directory,".gif")) # write to current dir
}

# don't use this at the moment, because it's not working
tree.through.time(lerista.rase, lerista.geo, lerista.dist, lerista, 
                  country="Australia", new.directory = "Lerista_Ancestors")

# resmc, goanna.geo.object
tree.through.time(rase.input = goanna.rase, geo.object = goanna.geo.object, phy = oz.g$phy, 
                  distribution.data = goanna.dist, map=rangemap, new.directory="Goanna_Ancestors")


list.files(path = "/Users/Ian/Documents/GitHub/MonitorPhylogenomics/Goanna_Ancestors", pattern = "*.pdf", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=1) %>% # animates, can opt for number of loops
  image_write(paste0(getwd(),"/",new.directory,"/",new.directory,".gif"))

testo <- ggmap(rangemap) + 
  stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                 geom = "polygon", data = curr.taxon) + theme(legend.position = 'none') + annotate("text", x=120, y=-40, label= paste0(10.3," MYA"), color="white", cex=10)

ggplotly(testo)

curr.taxon <- filter(goanna.dist, Name_in_Tree == "Varanus.eremius")
ggmap(rangemap) + stat_density2d(curr.taxon, aes(x="Longitude", y="Latitude"))


all.names <- data.frame(Names = c(unique(goanna.rase$DistData$Name_in_Tree), oz.g$phy$tip.label))
colorz <- colorRampPalette(brewer.pal(9, "YlGnBu"))(nrow(all.names))
names(colorz) <- all.names$Names

curr.color <- subset(colorz, names(colorz) %in% taxa.to.plot$Name_in_Tree)
ggmap(rangemap) + 
        stat_density2d(aes(x = Longitude, y = Latitude, fill = curr.color), alpha = .3, bins=4,
                       geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none')

testo <- taxa.to.plot$Name_in_Tree
curr.color[which(testo[1] == names(curr.color))]

apply(taxa.to.plot$Name_in_Tree, FUN=function(x) curr.color[which(testo[x] == names(curr.color))])
lapply(taxa.to.plot, function(x) curr.color[which(taxa.to.plot[x,"Name_in_Tree"] == names(curr.color))])


testa <- NULL
for (i in 1:nrow(taxa.to.plot)){
  valuez <- curr.color[which(taxa.to.plot[i, "Name_in_Tree"] == names(curr.color))]
  testa <- append(testa, valuez)
}

taxa.to.plot$color <- testa
ggmap(rangemap) + 
  stat_density2d(aes(x = Longitude, y = Latitude, fill = Name_in_Tree), alpha = .3, bins=4,
                 geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none') +
  scale_color_brewer(palette="Spectral")
opts <- options()
options(ggplot2.continuous.fill = "viridis")
options(ggplot2.continuous.colour="viridis")



