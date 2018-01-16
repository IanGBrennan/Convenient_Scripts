library(phytools)
library(btw)
library(ggplot2)

# tree manipulating and clade extraction #
tree<-read.tree("Pygopodoidea.tre")
keep<-c("Saltuarius.salebrosus", "Underwoodisaurus.seorsus",
        "Delma.australis", "Aprasia.picturata", "Pseudothecadactylus.australis",
        "Crenadactylus.ocellatus.KM.A", "Lucasium.alboguttatum", "Diplodactylus.klugei",
        "Strophurus.wilsoni", "Hesperoedura.reticulata")
plot(tree); nodelabels(cex=0.3)
pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, keep));
write.tree(pruned.tree, file="PGLS.pygo.tre")
plot(pruned.tree)
pleth<-read.tree("Plethodontidae.tre")
plot(pleth)

tree<-read.nexus("MSY.pygo.tre")
plot(tree, cex=0.3); nodelabels(cex=0.3)
msydip<-extract.clade(tree, node=210)
plot(msydip)
gammaStat(msydip)
1-pnorm(abs(gammaStat(msydip))) #one tailed T test


leaftails<-extract.clade(tree, node=159)
nephs<-extract.clade(tree, node=170)
delma<-extract.clade(tree, node=184)
pygos<-extract.clade(tree, node=205)
pseudo<-extract.clade(tree, node=224)
crena<-extract.clade(tree, node=227)
luc<-extract.clade(tree, node=240)
dip<-extract.clade(tree, node=255)
stroph<-extract.clade(tree, node=280)
oed<-extract.clade(tree, node=297)

  
#### Calculate the Gamma statistic (Pybus & Harvey)
tree<-read.tree("strophurus.tre")
gammaStat(tree) #return the gamma stat
1-pnorm(abs(gammaStat(tree))) #one tailed T test

#### Start by checking for correlation between traits #####
nocorrC <- Continuous(tree, data, tc = TRUE)
corrC <- Continuous(tree, data)
lrtest(nocorrC, corrC)

tree<-read.tree("PGLS.pygo.tre")
data<-read.csv("Gamma.vs.Size.csv")
pgls<-Continuous(tree, data, mode="ML", mlt=1000, regression=TRUE)
plot(data)

ggplot(data, aes(x=gamma, y=size.rate)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)




#### get MDI estimates for a discrete period
svl <- read.csv("body.biome.csv", row.names=1, header=TRUE) #data
tree<-read.tree("diplodactylus.tre") #assign your tree and load it
data<-setNames(svl[,5], rownames(svl)) #choose your data file (here: svl) and column ([,5]), apply rownames
source("dtt.full with confidence limits and p value.R"); ## now source the function from your WD
# and run the analysis. i recommend using ~10,000 sims to get a stable P value
##if you've used the regular dtt.full, you'll notice that this function doesn't immediately plot up the empirical dtt curve. don't worry - it's just that the way I get the nice shaded 95% region requires a blank plot at first with the empirical and median curves overlain on the shaded region.
dttOut<-dttFullCIs(tree, data, nsims=1000)
dttOut$MDI; dttOut$Pvalue
#dttOut<-dttFullCIs(pygopodoidea, data, nsims=1000, mdi.range=c(0,0.5)); #designate a time slice to measure the MDI, here first half of the tree
dttOut<-dttFullCIs(tree, data, nsims=1000, mdi.range=c(0,0.5)); #designate a time slice to measure the MDI, here second half of the tree
dttOut$MDI; dttOut$Pvalue
dttOut<-dttFullCIs(tree, data, nsims=1000, mdi.range=c(0.5,1)); #designate a time slice to measure the MDI, here second half of the tree
dttOut$MDI; dttOut$Pvalue

#### Start by checking for correlation between traits #####
tree<-read.tree("PGLS.pygo.tre")
data<-read.csv("Gamma.vs.MDI.csv")
nocorrC <- Continuous(tree, data, tc = TRUE)
corrC <- Continuous(tree, data)
lrtest(nocorrC, corrC)
ggplot(data, aes(x=second.MDI, y=gamma)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)
