library(ape) #call up APE
library(geiger) #call up Geiger, both are necessary to mess with the phylogeny
library(phytools)#will also be useful for some tree manipulation
library(laser)
library(BAMMtools)#fire that bitch up!
library(TESS)
library(diversitree)


pygopodoidea<-read.tree("Pygopodoidea.tre") #assign your tree and load it
svl<- read.csv("body.biome.csv", row.names=1, header=TRUE) #read in the data you want to address as an object "svl"
data<-setNames(svl[,5], rownames(svl)) #set your data(here: svl) and column ([,5]), apply rownames, and call it object "data"


## what is the length of the current color ramp?
n<-length(ACR$cols)
## change to grey scale
ACR$cols[1:n]<-grey(0:(n-1)/(n-1))
plot(ACR)## change to blue -> red
ACR$cols[1:n]<-colorRampPalette(c("red", "pink", "white", "lightblue", "navy"), space="Lab")(n)
plot(setMap(ACR,invert=TRUE), fsize = 0.4, lwd=3, legend = TRUE, type = "phylogram", outline=TRUE, direction="rightwards")

### To invert colors:
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}

plot(setMap(ACR,invert=TRUE), fsize = 0.3, lwd=3, legend = TRUE, type = "phylogram", outline=TRUE, direction="rightwards")
plot(setMap(ACR,invert=TRUE), fsize = 0.4, lwd=5, legend = TRUE, type = "fan" , outline=TRUE, direction="rightwards")




ACR <- contMap(pygopodoidea, data, method = "fastAnc", res=200, fsize = 0.7, lwd=7, legend = TRUE, type = "phylogram", plot = "TRUE", direction="rightwards")
