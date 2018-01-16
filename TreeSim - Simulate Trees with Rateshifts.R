library(TreeSim)
library(TESS)

#################################################################
# Simulate Tree with Rateshift Conditioned on #taxa
#################################################################

#sim.rateshift.taxa(n, numbsim, lamda, mu, frac, times, complete)
rfish.sim<-sim.rateshift.taxa(126, 2, c(0.21,0.21), c(0.27,0.1), c(1,0.5), c(0,0.5), complete=FALSE)
plot(rfish.sim[[1]])
write.tree(try, file="try.tre")
View(try)

x<-runif(100, 1, 56)
sort(x)

?sim.rateshift.taxa.help


rainbowfish<-read.tree("oz.rainbowfish.tre")
spiders<-read.tree("oz.archaeidae.tre")
crayfish<-read.tree("oz.crayfish.tre")
hydroporini<-read.tree("oz.diving.beetles.tre")
fungi<-read.tree("oz.fungi.tre")
galaxias<-read.tree("oz.galaxias.tre")
honeyeaters<-read.tree("oz.honeyeaters.pardalotes.tre")
bees<-read.tree("oz.hylaeus.tre")
marsupials<-read.tree("oz.marsupials.tre")
skinks<-read.tree("oz.sphenomorphines.tre")
geckos<-read.tree("Pygopodoidea.tre")
mltt.plot(rainbowfish, spiders, crayfish, hydroporini, fungi, galaxias, honeyeaters, bees, marsupials, skinks, geckos, log='y')
mltt.plot(rainbowfish, spiders, crayfish, hydroporini, fungi, galaxias, honeyeaters, bees, marsupials, skinks, geckos)




