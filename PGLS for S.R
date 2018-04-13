library(caper)
library(phytools)

## read in your data as a CSV (comma separated values) file
### this is assuming (header=T) that you have column names indicating the variable
rawdata<-read.csv("/PATH_TO_YOUR_DATA/frogdata.csv", header=TRUE)

## read in your tree file
### this is assuming your tree is in newick/nexus format
tree<-read.tree("/PATH_TO_YOUR_TREE/frogs.tre")

## use caper to set up your data properly, this will take into account your phylogeny
### when assessing your trait values
data<-comparative.data(phy=tree, data=rawdata, names.col=speciesName, vcv=TRUE, na.omit=FALSE)

## run the actual analysis, here, logSVL is the dependent variable, ecology is the independent
### meaning length (SVL) is dictated by the type of ecology
model.pgls<-pgls(logSVL~x.ecology, data=data, lambda="ML")

## have a look at the model outputs, more information can be found in the 
### 'PGLS' and 'Cooper - 2014 - Comparative Methods' documents I sent
summary(model.pgls)
plot(model.pgls)