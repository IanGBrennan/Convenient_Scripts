library(phytools)
library(ggplot2)

#####################################################################
#### Plot a Distribution as Binned Histogram
#####################################################################

data<-read.csv("ecological.matrix.csv", row.names=1, header=TRUE)
svl<-data["SVL"]
ggplot(svl, aes(SVL)) + 
  geom_histogram(binwidth=5) + 
  scale_x_continuous(breaks=c(0,50,100,150,200,250,300))
#plot the SVL data as binned histogram with x-axis ticks denoted by scale_x_continuous