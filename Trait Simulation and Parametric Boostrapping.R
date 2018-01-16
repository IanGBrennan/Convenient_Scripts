library(TreeSim)
library(diversitree)

tree<-read.tree("pygopodoidea.tre")

char.rates<-list(rbind(c(-.5,.5),c(.5,-.5)))
pygosimchar<-sim.char(tree, char.rates, nsim=1, model="discrete", root=1)




phy<-read.tree("Pygopodoidea.tre") #we've already loaded our tree, but why not just make sure
s<- read.csv("biome.simulated.csv", header=TRUE) #reads in the csv with header as the top row
sim<-s$arid #call your data column
names(sim)<-s$Species #attach the rownames (Species) to your arid data
p<-starting.point.geosse(phy) #make a starting pass for the GeoSSE analysis, this is gonna be dirty
p #have a quick look to make sure things worked, and your rough starting estimates
simfull.lik<-make.geosse(phy, sim) #estimate the full likelihood for the 
lik.con<-constrain(simfull.lik, dA~dB, sAB~0, xA~xB) #constrain dispersal directional dispersal as equal
lik.sp<-constrain(simfull.lik, sA~.125, sB~.035)

simm.full<-find.mle(simfull.lik, p) #lnLik = -685.7991 #this estimates the marginal likelihood of our model
ml.con<-find.mle(lik.con, p[argnames(lik.con)]) #lnLik = -702.6337 #get the marginal likelihood of alternate schemes (here lik1), it's worse(ml.d)
ml.sp<-find.mle(lik.sp, p[argnames(lik.sp)])

coef(simm.full)

simm.full$lnLik-ml.con$lnLik
simm.full$lnLik-ml.sp$lnLik
