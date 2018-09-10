library(ratematrix); library(phytools)
data(anoles)

cols <- setNames(c("greenyellow", "plum"), c("mainland", "island"))
plotSimmap(tree=anoles$phy, colors=cols, fsize=0.6, lwd=3)
estimateTimeMCMC(data=anoles$data[,2:3], phy=anoles$phy, gen=100000)
handle <- ratematrixMCMC(data=anoles$data[,2:3], phy=anoles$phy, prior="empirical_mean", start="mle"
                         , gen=100000, outname="anoles_example")

plot(anoles$phy)
astral.raxml <- read.tree ("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Final_Trees/Marsupials_ASTRAL_RAxML.trees")
trait.data <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Macro_Inference/Trait_Data/Marsupials.ASTRAL.RAxML.PhyloResiduals.rds")

estimateTimeMCMC(data=trait.data[[1]], phy=astral.raxml[[1]], gen=100000)
handle <- ratematrixMCMC(data=trait.data[[1]], phy=astral.raxml[[1]], prior="empirical_mean", start="mle"
                         , gen=100000, outname="mars_example")
post <- readMCMC(handle=handle, burn=0.1, thin=100)
checkConvergence(post)
logAnalyzer(handle=handle, burn=0.1, thin=100)
plotRatematrix(chain=post, colors=c("greenyellow", "plum"), l.cex=1, set.xlim=c(0,0.02))
