library(plotly)
library(transformr)
library(gganimate)

Sys.setenv("plotly_username" = "igbrennan")
Sys.setenv("plotly_api_key" = "M0PqBWJaZcwriWRnd45O")

# apparently the Plotly API is limited to 50 calls a day, so be judicious


########################################################
### An example of the density map function with earthquakes
########################################################
quakes = read.csv('https://raw.githubusercontent.com/plotly/datasets/master/earthquakes-23k.csv')
p <- quakes %>%
  plot_ly(
    type = 'densitymapbox',
    lat = ~Latitude,
    lon = ~Longitude,
    coloraxis = 'coloraxis',
    radius = 10) %>%
  layout(
    mapbox = list(
      style="stamen-toner",
      center= list(lon=180)), coloraxis = list(colorscale = "Viridis"))

### and if you wanted to try it for real data
vg <- filter(taxa.to.plot, Name_in_Tree == "Varanus.glebopalma")

p <- vg %>%
  plot_ly(
    type = "densitymapbox",
    lat = ~Latitude,
    lon = ~Longitude,
    coloraxis = "coloraxis",
    radius = 10) %>%
  layout(
    mapbox = list(
      style="stamen-terrain",
      center= list(lon=180)), coloraxis = list(colorscale = "Viridis"))



library(gapminder)
library(gganimate) # might be useful to come back to this at some point.
library(gifski)

goannas <- filter(alldata, Group == "Varanoidea" && Status == "Extant")
goannas <- filter(alldata, Group == "Varanoidea", is.na(Tail)==FALSE, is.na(SVL)==FALSE)
goannas <- mutate(goannas, TailRatio = Tail/SVL)


########################################################
### Interactive plot of monitor body sizes
########################################################
sizes <- plot_ly(goannas, 
        x = ~SVL, 
        y = ~Tail, 
        text = ~paste("Species: ", Name_in_Tree,
                      "<br> Habitat: ", Habitat), 
        type = 'scatter', 
        mode = 'markers', 
        color = ~TailRatio, 
        colors = 'YlGnBu',
        #frame = ~Habitat,
        marker = list(size = ~TailRatio*20, opacity = 0.75, sizemode = "diameter",
                      line = list(width = 2, color = "black"))) %>%
  layout(title = 'Monitor Lizard Body and Tail Length',
         xaxis = list(type = "log"),
         yaxis = list(type = "log"))

api_create(sizes, filename = "Varanus_SVL_vs_Tail")


########################################################
### Use ggplotly and ggtree to plot an interactive tree
########################################################

p1 <- ggtree(oz.g$phy)

id <- oz.g$phy$tip.label
set.seed(42)
grp <- sample(LETTERS[1:4], size = Ntip(oz.g$phy), replace = T)
dat <- tibble::tibble(id = id,
                      grp = grp)

metat <- p1$data %>%
  dplyr::inner_join(dat, c('label' = 'id'))

p2 <- p1 +
  geom_point(data = metat,
             aes(x = x,
                 y = y,
                 colour = grp,
                 label = id))

plotly::ggplotly(p2, tooltip = "label")


########################################################
### Plot the distribution of multiple taxa through time
########################################################
vtest <- filter(taxa.to.plot, Name_in_Tree == "Varanus.glebopalma" | Name_in_Tree == "n40")
vtest <- mutate(vtest, time = 10)
vtest2 <- vtest
vtest2$Latitude <- vtest2$Latitude-2; vtest2$Longitude <- vtest2$Longitude-1; vtest2$time <- vtest2$time-5

vtestall <- rbind(vtest, vtest2)

vtestall %>%
  plot_ly(
    type = "densitymapbox",
    lat = ~Latitude,
    lon = ~Longitude,
    color = ~Name_in_Tree,
    frame = ~time,
    #coloraxis = "coloraxis",
    radius = 10) %>%
  layout(
    mapbox = list(
      style="stamen-toner",
      center= list(lon=180)), coloraxis = list(colorscale = "Viridis"))


library(datasauRus)

ggplot(datasaurus_dozen, aes(x=x, y=y))+
  geom_point()+
  theme_minimal() +
  transition_states(dataset, 3, 1) + 
  ease_aes('cubic-in-out')


#############################################################
### Animate the relationship between variables among groups
#############################################################

gramp <- colorRampPalette(brewer.pal(3, "YlGnBu"))
gcolors <- gramp(nrow(goannas)); 
#goannas$Name_in_Tree <- factor(goannas$Name_in_Tree, levels = goannas$Name_in_Tree[order(goannas$SVL)])
goannas <- goannas[order(goannas$SVL),]
goannas <- mutate(goannas, plotcols = gcolors) # could've also just done 'goannas$plotcols <- gcolors'

TailSVL <- ggplot(goannas, aes(x=log(SVL), y=log(Tail), size=TailRatio)) + # using 'color=Name_in_Tree' will lose smooth animations
  geom_point(alpha=0.75, show.legend=FALSE, color="black", fill=goannas$plotcols, shape=21) + # designate color here to keep the smoothness
  #scale_colour_manual(values = gcolors) +
  theme_bw() +
  scale_size(range=c(5,20)) +
  labs(title="Habitat: {closest_state}") + # this will label the plot by the groupings
  transition_states(Habitat, transition_length=5, state_length=3) +
  ease_aes('cubic-in-out')
anim_save(TailSVL, file="~/Google.Drive/ANU/AHE/T474_Varanus/Goanna_Figures/TailSVL.gif")



ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  ease_aes('linear')


# Plot distributions through time

# I'm going to need to start by going through the geo.object
## for each geo.object[i], I need to get the taxa, and bind
## them to a new data frame with their distribution data 
## and most importantly, the time they are at!

# start by combining 
alldist <- rbind(goanna.dist, goanna.rase$DistData[,1:3])
phy <- oz.g$phy

plot.times <- tutorial$geo.object$times
geo.object <- tutorial$geo.object
rase.input <- tutorial$rase.data
colramp <- colorRampPalette(brewer.pal(6, "Reds"))
      colorz <- colramp(Ntip(phy) + (Ntip(phy)-1))
          names(colorz) <- c(phy$tip.label, names(tutorial$rase.data$ConvexHulls))



#for (k in 1:length(geo.object$geography.object)) {
thedist <- NULL
for (k in 1:10) {
  current.time <- round(max(plot.times) - plot.times[[k+1]],digits=1)
  time.minus1  <- round(max(plot.times) - plot.times[[k]],digits=1)
  #print(paste("plotting distributions at", current.time, "time units ago"))
  taxa.to.plot <- NULL
  extant.names <- rownames(geo.object$geography.object[[k]])
  
## THIS DOESN'T YET WORK BECAUSE I NEED TO RECONCILE THE RPANDA NAMES (.AAA, .AAB)
## WITH THE NODE NAMES I'VE GIVEN THEM (n67, n68)
## I SHOULD BE ABLE TO DO THIS STEALING SCRIPT FROM CreateGeoObject_SP
#  for (j in 1:length(extant.names)){
#    current.taxon <- filter(alldist, Name_in_Tree == extant.names[[j]])
#    current.taxon$time <- current.time
#    taxa.to.plot <- rbind(taxa.to.plot, current.taxon)
#  }
  
## I THINK I CAN JUST USE THE CODE I STOLE FROM tree.through.time.R   
  for (j in 1:length(extant.names)) {
    if(extant.names[[j]] %in% phy$tip.label){
      current.taxon <- filter(alldist, Name_in_Tree == extant.names[[j]])
      current.taxon$time <- -current.time
      current.taxon$color <- subset(colorz, names(colorz) == extant.names[[j]])
      #taxa.to.plot <- rbind(taxa.to.plot, filter(distribution.data, Name_in_Tree == extant.names[[j]]))
      taxa.to.plot <- rbind(taxa.to.plot, current.taxon)
      
    # I added this bit as a potential resolution, but it didn't work how I thought it would  
    #  ancestor <- paste0("n", geo.object$name.matrix[which(geo.object$name.matrix[,2]==extant.names[[j]]),1])
    #  ancestor.taxon <- filter(rase.input$DistData, Name_in_Tree == ancestor)[,1:3]
    #  ancestor.taxon$time <- -time.minus1
    #  ancestor.taxon$color <- subset(colorz, names(colorz) == extant.names[[j]])
    #  taxa.to.plot <- rbind(taxa.to.plot, ancestor.taxon)
    #  print(paste("k=",k, "j=",j))
    } 
    else if(!extant.names[[j]] %in% phy$tip.label) {
      node.name <- paste0("n", geo.object$name.matrix[which(geo.object$name.matrix[,2]==extant.names[[j]]),3])
      #node.name <- paste0("n",node.name)
      current.taxon <- filter(rase.input$DistData, Name_in_Tree == node.name)[,1:3]
      current.taxon$time <- -current.time
      current.taxon$color <- subset(colorz, names(colorz) == node.name)
      taxa.to.plot <- rbind(taxa.to.plot, current.taxon)
      
    # I added this bit as a potential resolution, but it didn't work how I thought it would  
    #  ancestor <- paste0("n", geo.object$name.matrix[which(geo.object$name.matrix[,2]==extant.names[[j]]),1])
    #  if(ancestor == paste0("n", Ntip(phy))){
    #    ancestor.taxon <- filter(rase.input$DistData, Name_in_Tree == paste0("n",(Ntip(phy)+1)) | Name_in_Tree == paste0("n",(Ntip(phy)+2)))[,1:3]
    #    }
    #  else if(!ancestor == paste0("n", Ntip(phy))) {
    #    ancestor.taxon <- filter(rase.input$DistData, Name_in_Tree == ancestor)[,1:3]
    #  }
    #  ancestor.taxon$Name_in_Tree <- ancestor
    #  ancestor.taxon$time <- -time.minus1
    #  ancestor.taxon$color <- subset(colorz, names(colorz) == ancestor)
    #  taxa.to.plot <- rbind(taxa.to.plot, ancestor.taxon)
    #  #print(paste("k=",k, "j=",j))
    } else {
      stop("current taxon does not have distributional data")
    }
  }
  thedist <- rbind(thedist, taxa.to.plot)
}
#thedist$time <- -thedist$time
head(thedist); tail(thedist)

#curr.color <- subset(colorz, names(colorz) %in% taxa.to.plot$Name_in_Tree)
##current.time <- paste0(round((max(nodeHeights(phy))-plot.times[[k+1]]),digits=1)," MYA")
## remember R will go "1,10,11,...,2,20,21.." so rename anything less than 10 as 01,02...
#pdf(paste(dir.path, "/", if(k<10){paste0("0",k)}else{k}, ".pdf", sep=""))
#print(ggmap(rangemap) + 
#        stat_density2d(aes(x = Longitude, y = Latitude, fill = curr.color), alpha = .3, bins=4,
#                       geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none') + annotate("text", x=120, y=-40, label= paste0(current.time," MYA"), color="white", cex=10))
#dev.off()

#ggmap(rangemap) + 
#        stat_density2d(aes(x = Longitude, y = Latitude, fill = curr.color), alpha = .3, bins=4,
#                       geom = "polygon", data = taxa.to.plot) + theme(legend.position = 'none') + annotate("text", x=120, y=-40, label= paste0(current.time," MYA"), color="white", cex=10)

rangemap <- get_googlemap(center = "Australia", zoom = 4, style = 'feature:all|element:labels|visibility:off')

ggmap(rangemap) +
  stat_density2d(aes(x = Longitude, y = Latitude, fill = color), # using 'fill = color' breaks the smooth transitions
                 alpha = 0.3, bins = 4,
                 geom = "polygon", data = thedist) + theme(legend.position = 'none') +
  labs(title="Time: {closest_state}") + # this will label the plot by the groupings
  transition_states(time, transition_length=5, state_length=3) +
  ease_aes('cubic-in-out')
  
ggmap(rangemap) +
  stat_density2d(aes(x = Longitude, y = Latitude, fill=Name_in_Tree), # using 'fill = color' breaks the smooth transitions
                 alpha = 0.3, bins = 4,
                 geom = "polygon", data = filter(alldist, Name_in_Tree == "Varanus.glebopalma")) + 
  scale_fill_manual(values = "blue") +
  theme(legend.position = 'none')

### Below this is just a test of the method
### It's the best looking version really, two ranges shrinking
tempdist <- thedist
td1 <- filter(tempdist, time == -18.3 & Name_in_Tree == "n35")
td1$Name_in_Tree <- "Varanus.glebopalma"; td1$color <- colorz[which(names(colorz)=="Varanus.glebopalma")]
td2 <- filter(tempdist, time == -18.3 & Name_in_Tree == "n40")
td3 <- filter(tempdist, time == -16 & Name_in_Tree == "Varanus.glebopalma")
td4 <- filter(tempdist, time == -16 & Name_in_Tree == "n37")
td4$Name_in_Tree <- "n40"; td4$color <- td2[1,"color"]
#td5 <- filter(tempdist, time == )

testdist <- rbind(td1, td2, td3, td4)

ggmap(rangemap) +
  stat_density2d(aes(x = Longitude, y = Latitude, fill = color), alpha = 0.3, bins = 4,
                 geom = "polygon", data = testdist) + theme(legend.position = 'none') +
  scale_fill_manual(values = unique(testdist$color)) +
  labs(title="Time: {closest_state}") + # this will label the plot by the groupings
  transition_states(time, transition_length=5, state_length=3) +
  ease_aes('cubic-in-out')

### I could try to follow each taxon from the root to its tip
### Mapping the distribution, but it could be really messy

# we can start with an example for the acanthurus group
acantree <- extract.clade(oz.g$phy, node=43)

ggmap(rangemap) +
  stat_density2d(aes(x = Longitude, y = Latitude, fill = Color), alpha = 0.3, bins = 4,
                 geom = "polygon", data = ranges) + theme(legend.position = 'none') +
  scale_fill_manual(values = unique(ranges$Color)) +
  labs(title="Time: {closest_state}") + # this will label the plot by the groupings
  transition_states(time, transition_length=5, state_length=3) +
  ease_aes('cubic-in-out')


