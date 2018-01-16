tree <- read.nexus("BT.Pygopodoidea.tre")
carpygo <- findMRCA(tree, c("P_Delma_australis", "P_Nephrurus_vertebralis"))
diplo <- findMRCA(tree, c("P_Pseudothecadactylus_australis", "P_Strophurus_congoo"))
dip.tree <- extract.clade(tree, diplo)
carpygo.tree <- extract.clade(tree, carpygo)

add.carpygo <- 150 - max(nodeHeights(carpygo.tree))
add.dip     <- 150 - max(nodeHeights(dip.tree))


new.carpygo <- addroot(carpygo.tree, add.carpygo)
new.diplo   <- addroot(dip.tree, add.dip)
long.tree <- bind.tree(new.carpygo, new.diplo, where="root")
plot(long.tree, show.tip.label=F)
write.nexus(long.tree, file="Long.stem.Pygopodoidea.tre")

data <- read.csv("BT.Pygopodoidea.logSVL.csv", header=F, row.names=1)
BM <- fitContinuous(tree, data, model="BM")
longBM <- fitContinuous(long.tree, data, model="BM")

eco.data <- read.csv("OU.Pygopodoidea.Relics.csv", header=T) # read it in for OUwie format
ecology <- eco.data[,2]
names(ecology) <- eco.data[,1]

sim.tree <- make.simmap(tree, ecology, model="ER", nsim=1)
colorz <- c("blue", "green", "red")
names(colorz) <- c(0,1,2)
plotSimmap(sim.tree, colorz, pts=F, ftype="off", lwd=1)

sim.long.tree <- make.simmap(long.tree, ecology, model="ER", nsim=1)
colorz <- c("blue", "green", "red")
names(colorz) <- c(0,1,2)
plotSimmap(sim.long.tree, colorz, pts=F, ftype="off", lwd=1)

test.bms <- OUwie(sim.tree, eco.data, model="BMS", root.station=T, simmap.tree=T)
test.long.bms <- OUwie(sim.long.tree, eco.data, model="BMS", root.station=T, simmap.tree=T)

test.ouma <- OUwie(sim.tree, eco.data, model=c("OUMA"), root.station=T, simmap.tree = T)
test.long.ouma <- OUwie(sim.long.tree, eco.data, model=c("OUMA"), root.station=T, simmap.tree = T)

test.oum <- OUwie(sim.tree, eco.data, model=c("OUM"), root.station=T, simmap.tree = T)
test.long.oum <- OUwie(sim.long.tree, eco.data, model=c("OUM"), root.station=T, simmap.tree = T)









