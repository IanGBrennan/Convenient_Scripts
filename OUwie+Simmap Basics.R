library(ape)
library(phytools)
library(OUwie)

tree <- read.nexus("BT.Pygopodoidea.tre")
data <- read.csv("OU.Pygopodoidea.Relics.csv", header=T) # read it in for OUwie format

## make.simmap can be annoying, follow directions
ecology <- data[,2] # pull out the character coding, from the input file
names(ecology) <- data[,1] # attach names to the data from the input file

sim.tree <- make.simmap(tree, ecology, model="ER", nsim=1)
colorz <- c("blue", "green", "red")
names(colorz) <- c(0,1,2)
plotSimmap(sim.tree, colorz, pts=F, ftype="off", lwd=1)

test.bms <- OUwie(sim.tree, data, model=c("BMS"), root.station=T, simmap.tree = T)
test.ouma<- OUwie(sim.tree, data, model=c("OUMA"), root.station=T, simmap.tree = T)
test.oum <- OUwie(sim.tree, data, model=c("OUM"), root.station=T, simmap.tree = T)
test.oumv<- OUwie(sim.tree, data, model=c("OUMV"), root.station=T, simmap.tree = T)
test.oumva<- OUwie(sim.tree, data, model=c("OUMVA"), root.station=T, simmap.tree = T)
aicc <- NULL
aicc <- as.data.frame(c(test.bms$AICc, test.ouma$AICc, test.oum$AICc, test.oumv$AICc, test.oumva$AICc))
rownames(aicc) <- c("bms", "ouma", "oum", "oumv", "oumva"); colnames(aicc) <- "AICc"
aicc
# for the Pygopodoidea, the best fitting model is 'OUMV' allowing different rates and optima per family
# for comparison of Aprasia to the rest of the Pygopodoidea, the best model is equally 'BMS' or 'OUMVA' (need bigger Aprasia tree)
# for comparison of Crenas  to the rest of the Pygopodoidea, the best model is 
# for comparison of relic geckos to the rest of the Pygopodoidea, the best model is 

### Now test Aprasia differ from other pygopods
plot(sim.tree, cex=0.3); nodelabels(cex=0.2, frame="circle")
pygo.tree <- extract.clade(tree, 192)
sim.pygo <- make.simmap(pygo.tree, ecology, model="ER", nsim=1)
pygo.data <- subset(data, data["regime"]==1)
pygo.data <-rbind.data.frame(pygo.data, (subset(data, data["regime"]==3)))

bms  <- OUwie(sim.pygo, pygo.data, model=c("BMS"), root.station=T, simmap.tree=T)
oum  <- OUwie(sim.pygo, pygo.data, model=c("OUM"), root.station=T, simmap.tree=T)
ouma <- OUwie(sim.pygo, pygo.data, model=c("OUMA"), root.station=T, simmap.tree=T)
oumv <- OUwie(sim.pygo, pygo.data, model=c("OUMV"), root.station=T, simmap.tree=T)
oumva <- OUwie(sim.pygo, pygo.data, model=c("OUMVA"), root.station=T, simmap.tree=T)

bms$AICc; oum$AICc; ouma$AICc; oumv$AICc; oumva$AICc

