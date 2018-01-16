library(ggplot2)
data("mtcars")
mtcars
data("diamonds")
tail(diamonds)

data<-read.csv("Pygopodoidea.inward.outward.Transitions.csv")
freq.out<-data$freq.out
freq.in<-data$freq.in
biome.out<-data$biome.out
biome.in<-data$biome.in

ggplot(data, aes(biome.out, fill=freq.out)) + geom_bar()
ggplot(data, aes(biome.in, fill=freq.in)) + geom_bar()
