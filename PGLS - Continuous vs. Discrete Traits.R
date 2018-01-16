library(caper)
library(phytools)

data <- read.csv("ecological.matrix.csv", row.names=1, header=TRUE) #data
pygo.data<-read.csv("ecological.matrix.csv", header=TRUE)
tree<-read.tree("pygopodoidea.tre") #assign your tree and load it
logsvl<-setNames(data[,13], rownames(data)) #choose your data file (here: svl) and column ([,5]), apply rownames
ecology<-setNames(data[,6], rownames(data)) #choose your data file (here: svl) and column ([,5]), apply rownames
species<-pygo.data[,1]

pygo<-comparative.data(phy=tree, data=pygo.data, names.col=speciesName, vcv=FALSE, na.omit=FALSE)
model.pgls<-pgls(logSVL~x.diel, data=pygo, lambda="ML")
summary(model.pgls)

plot(logSVL~combo.eco, data=pygo$data)
abline(model.pgls)


pygo.data$x.diel
fish<-as.factor(pygo.data$x.diel)
