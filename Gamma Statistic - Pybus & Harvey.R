library(ape)
library(laser)
library(geiger)

#Easier way in LASER (can ignore everything after this)
################################################
tree<-read.tree("Pygopodoidea.tre")
gamStat(branching.times(tree), return.list=TRUE)
################################################
tree<-read.tree("oz.marsupials.tre")
gamStat(branching.times(tree), return.list=TRUE)





#### Calculate the Gamma statistic (Pybus & Harvey)
#then compare it to a null hypothesis (constant rate of diversification)
tree<-read.tree("Pygopodoidea.tre")
gammaStat(tree) #return the gamma stat
pygogamma<-2*(1-pnorm(abs(gammaStat(tree)))) #two tailed T test
pygogamma1<-1-pnorm(abs(gammaStat(tree))) #one tailed T test

pygo<-read.tree("pygopodidae.tre")
gammaStat(pygo)
pygamma<-2*(1-pnorm(abs(gammaStat(pygo))))
pygamma1<-1-pnorm(abs(gammaStat(pygo)))

carpho<-read.tree("carphodactylidae.tre")
gammaStat(carpho)
cargamma<-2*(1-pnorm(abs(gammaStat(car))))
cargamma1<-1-pnorm(abs(gammaStat(carpho)))

diplo<-read.tree("diplodactylidae.tre")
gammaStat(diplo)
pygamma<-2*(1-pnorm(abs(gammaStat(pygo))))
pygamma1<-1-pnorm(abs(gammaStat(diplo)))

cdiplo<-read.tree("core.diplodactylidae.tre")
gammaStat(cdiplo)
pygamma<-2*(1-pnorm(abs(gammaStat(pygo))))
pygamma1<-1-pnorm(abs(gammaStat(cdiplo)))

carpygo<-read.tree("carphodactylidae.pygopodidae.tre")
gammaStat(carpygo)
pygamma<-2*(1-pnorm(abs(gammaStat(carpygo))))
pygamma1<-1-pnorm(abs(gammaStat(carpygo)))
