# function 'bind.at.depth' to choose where a new tip goes into the phylogeny
# might be useful for 

# extinction is an additionally interesting derivative of missing sampling. 
# in studies that aren't looking at temporal patterns in trait change, you could
# simply add tips to an existing phylogeny ('bind.tip' or 'bind.at.depth').
# however, an extinct taxa is sampled in a given time period, then may not 
# be sampled in the subsequent period, changing the dynamics. 

library(phytools)
library(phangorn)
library(geiger)
library(mvtnorm)
library(MuMIn)
library(BioGeoBEARS)
library(diversitree)
library(Rmisc)
library(paleotree)
library(ggplot2); library(wesanderson)
library(mvMORPH); library(parallel)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Sim.Fossil.Trees.SOURCE.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD


#######################################################
# Read in the empirical trees you'll use
######################################################
agam <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.tre") # 100 tips
mars <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.tre") # 133 tips
bird <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagides.tre") # 149 tips
pygo <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Pygopodoidea.tre") # 189 tips
skink<- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Skinks.1.tre") # 240 tips
# or read in fossilized trees you've already made
stochastic <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Stochastic.Fossil.trees")
miocene    <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Miocene.Fossil.trees")
pliopleis  <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.PlioPleistocene.Fossil.trees")

#######################################################
# Attach extinct tips to your empirical phylogeny
######################################################
base.tree <- skink
## simulate trees under the "stochastic extinction" concept
#sim.extinct <- sim.fossil.trees(phy=base.tree, time.frame=c(0.1,max(nodeHeights(base.tree)),0.1), 
#                                num.taxa=(length(base.tree$tip.label))/2, num.trees=3) # output is called 'output.trees'
sim.extinct <- alt.sim.fossil.trees(phy=base.tree, time.frame=c(10,max(nodeHeights(base.tree)),0.1),
                                    die=c(5,10,0.1), 
                                    num.taxa=(length(base.tree$tip.label))*.50, num.trees=10) # output is called 'output.trees'

## Designate the trees you'd like to use for the data simulations and model fitting
###################################################################################
input.trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.PlioPleistocene.Fossil.number.trees")


######################################################
## Set empty object to hold your results
sim.traits.geiger <- list(); sim.traits.ouwie <- list() # make trait lists for all trees in geiger and ouwie data format
save.sim.traits <- NULL
extant.data <- list()
# make sure to the outputs depending on if you're simulating different times, or alphas
# this means changing the 'm' and the 'm[[2]]' objects below!
num.sims <- 1 # designate the number of simulated data sets you like to create per tree


######################################################################################
# Simulate data under many different processes using mvMORPH
# simulate onto the extinct trees, then fit models! Drop tips/data the fit again!
######################################################################################
# or read in fossilized trees you've already made
stochastic <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Stochastic.Fossil.trees")
miocene    <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Miocene.Fossil.trees")
pliopleis  <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.PlioPleistocene.Fossil.trees")

trees <- miocene
# trees <- multiphylo.era.simmap(trees, shift.time) # if simulating with mvMORPH!

bm.pars <- 0.1 # set the diffusion parameter of the BM process
trend.pars <- sample(seq(from=0.1, to=0.5, by=0.01), 100, replace=T)
shift.time <- 10

all.models.lik <- NULL
model.fit <- NULL

group.name <- "Extinct"
save.path <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/"

# if simulating using Phytools for BMtrend
simulated.traits <- mclapply(1:length(trees), function(x){
  (fastBM(trees[[x]], mu=trend.pars[[x]], sig2=bm.pars, a=0, nsim=1))})
# if simulating using mvMORPH for BMOUi or BMOU
simulated.traits <- mclapply(1:length(trees), function(x){
  mvSIM(trees[[x]], nsim=1, model="BMOUi", param=agam.i$BMOUi[[x]])}, mc.cores=8)


# Fit and Process the OU model
OU1_list <- mclapply(1:length(trees), function(x){
  mvOU(trees[[x]], simulated.traits[[x]], model="OU1", method="sparse", # error=data.me^2
       diagnostic=F, echo=F)}, mc.cores = 8)
OU_res <- NULL; for(y in 1:length(OU1_list)) {
  OU_temp <- as.data.frame(t(c(OU1_list[[y]]$LogLik,OU1_list[[y]]$AIC,OU1_list[[y]]$AICc, y, "OU1", NA, group.name)))
  OU_res <- rbind(OU_res, OU_temp)
}
all.models.lik <- rbind(all.models.lik, OU_res)
model.fit[["OU"]] <- OU1_list

# Fit and Process the  BM model
BM1_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], simulated.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
BM_res <- NULL; for(y in 1:length(BM1_list)) {
  BM_temp <- as.data.frame(t(c(BM1_list[[y]]$LogLik,BM1_list[[y]]$AIC,BM1_list[[y]]$AICc, y, "BM1", NA, group.name)))
  BM_res <- rbind(BM_res, BM_temp)
}
all.models.lik <- rbind(all.models.lik, BM_res)
model.fit[["BM"]] <- BM1_list

# Fit and Process the EB model
EB_list <- mclapply(1:length(trees), function(x){
  mvEB(trees[[x]], simulated.traits[[x]], method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
EB_res <- NULL; for(y in 1:length(EB_list)) {
  EB_temp <- as.data.frame(t(c(EB_list[[y]]$LogLik,EB_list[[y]]$AIC,EB_list[[y]]$AICc, y, "EB", NA, group.name)))
  EB_res <- rbind(EB_res, EB_temp)
}
all.models.lik <- rbind(all.models.lik, EB_res)
model.fit[["EB"]] <- EB_list

# Fit and Process the BMOU model
BMOU_list <- mclapply(1:length(trees), function(x){
  #mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-cheese.time[[x]]))), simulated.traits[[x]], model="BMOU", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
  mvSHIFT(trees[[x]], simulated.traits[[x]], model="BMOU", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)

BMOU_res <- NULL; for(y in 1:length(BMOU_list)) {
  #BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", cheese.time[[y]], BMOU_list[[y]]$shift.time, group.name)))
  BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", shift.time, BMOU_list[[y]]$shift.time, group.name)))
  BMOU_res <- rbind(BMOU_res, BMOU_temp)
}
all.models.lik <- rbind(all.models.lik, BMOU_res)
model.fit[["BMOU"]] <- BMOU_list

# Fit and Process the BMOUi model
BMOUi_list <- mclapply(1:length(trees), function(x){
  #mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-cheese.time[[x]]))), simulated.traits[[x]], model="BMOUi", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
  mvSHIFT(trees[[x]], simulated.traits[[x]], model="BMOUi", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)

BMOUi_res <- NULL; for(y in 1:length(BMOUi_list)) {
  #BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", cheese.time[[y]], BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", shift.time, BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_res <- rbind(BMOUi_res, BMOUi_temp)
}
all.models.lik <- rbind(all.models.lik, BMOUi_res)
model.fit[["BMOUi"]] <- BMOUi_list

# Fit and Process the BMtrend model
Trend_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], simulated.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F, param=list(trend=TRUE))}, mc.cores = 8)
Trend_res <- NULL; for(y in 1:length(Trend_list)) {
  Trend_temp <- as.data.frame(t(c(Trend_list[[y]]$LogLik,Trend_list[[y]]$AIC,Trend_list[[y]]$AICc, y, "BMtrend", NA, group.name)))
  Trend_res <- rbind(Trend_res, Trend_temp)
}
all.models.lik <- rbind(all.models.lik, Trend_res)
model.fit[["Trend"]] <- Trend_list

######################################################################################
### Now do the same, but removing all the extinct tips from the trees and data!
######################################################################################

group.name <- "Extant"
    
#drops <- is.extinct(phy, tol=0.0001) # find out which tips are extinct
extant.trees<-NULL; for (k in 1:length(trees)){
  extant.trees[[k]] <- drop.extinct(trees[[k]], tol=0.00001)
} 
class(extant.trees) <- "multiPhylo"; trees <- extant.trees
# trees <- multiphylo.era.simmap(trees, shift.time) # if simulating with mvMORPH!

extant.traits<-NULL; for (j in 1:length(simulated.traits)){
  #extant.traits[[j]] <- subset(simulated.traits[[j]], names(simulated.traits[[j]]) %in% trees[[j]]$tip.label)
  extant.traits[[j]] <- subset(simulated.traits[[j]], rownames(simulated.traits[[j]]) %in% trees[[j]]$tip.label)
}
    
# Fit and Process the OU model
OU1_list <- mclapply(1:length(trees), function(x){
  mvOU(trees[[x]], extant.traits[[x]], model="OU1", method="sparse", # error=data.me^2
       diagnostic=F, echo=F)}, mc.cores = 8)
OU_res <- NULL; for(y in 1:length(OU1_list)) {
  OU_temp <- as.data.frame(t(c(OU1_list[[y]]$LogLik,OU1_list[[y]]$AIC,OU1_list[[y]]$AICc, y, "OU1", NA, group.name)))
  OU_res <- rbind(OU_res, OU_temp)
}
all.models.lik <- rbind(all.models.lik, OU_res)
model.fit[["OU"]] <- OU1_list

# Fit and Process the  BM model
BM1_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], extant.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
BM_res <- NULL; for(y in 1:length(BM1_list)) {
  BM_temp <- as.data.frame(t(c(BM1_list[[y]]$LogLik,BM1_list[[y]]$AIC,BM1_list[[y]]$AICc, y, "BM1", NA, group.name)))
  BM_res <- rbind(BM_res, BM_temp)
}
all.models.lik <- rbind(all.models.lik, BM_res)
model.fit[["BM"]] <- BM1_list

# Fit and Process the EB model
EB_list <- mclapply(1:length(trees), function(x){
  mvEB(trees[[x]], extant.traits[[x]], method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
EB_res <- NULL; for(y in 1:length(EB_list)) {
  EB_temp <- as.data.frame(t(c(EB_list[[y]]$LogLik,EB_list[[y]]$AIC,EB_list[[y]]$AICc, y, "EB", NA, group.name)))
  EB_res <- rbind(EB_res, EB_temp)
}
all.models.lik <- rbind(all.models.lik, EB_res)
model.fit[["EB"]] <- EB_list

# Fit and Process the BMOU model
BMOU_list <- mclapply(1:length(trees), function(x){
  #mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-cheese.time[[x]]))), extant.traits[[x]], model="BMOU", method="sparse", diagnostic=F, echo=F)}, mc.cores = 8)
  mvSHIFT(trees[[x]], extant.traits[[x]], model="BMOU", method="sparse", diagnostic=F, echo=F)}, mc.cores = 8)
BMOU_res <- NULL; for(y in 1:length(BMOU_list)) {
  #BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", cheese.time[[y]], BMOU_list[[y]]$shift.time, group.name)))
  BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", shift.time, BMOU_list[[y]]$shift.time, group.name)))
  BMOU_res <- rbind(BMOU_res, BMOU_temp)
}
all.models.lik <- rbind(all.models.lik, BMOU_res)
model.fit[["BMOU"]] <- BMOU_list

# Fit and Process the BMOUi model
BMOUi_list <- mclapply(1:length(trees), function(x){
  #mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-cheese.time[[x]]))), extant.traits[[x]], model="BMOUi", method="sparse", diagnostic=F, echo=F)}, mc.cores = 8)
  mvSHIFT(trees[[x]], extant.traits[[x]], model="BMOUi", method="sparse", diagnostic=F, echo=F)}, mc.cores = 8)
BMOUi_res <- NULL; for(y in 1:length(BMOUi_list)) {
  #BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", cheese.time[[y]], BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", shift.time, BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_res <- rbind(BMOUi_res, BMOUi_temp)
}
all.models.lik <- rbind(all.models.lik, BMOUi_res)
model.fit[["BMOUi"]] <- BMOUi_list

# Fit and Process the BMtrend model
Trend_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], extant.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F, param=list(trend=TRUE))}, mc.cores = 8)
Trend_res <- NULL; for(y in 1:length(Trend_list)) {
  Trend_temp <- as.data.frame(t(c(Trend_list[[y]]$LogLik,Trend_list[[y]]$AIC,Trend_list[[y]]$AICc, y, "BMtrend", NA, group.name)))
  Trend_res <- rbind(Trend_res, Trend_temp)
}
all.models.lik <- rbind(all.models.lik, Trend_res)
model.fit[["Trend"]] <- Trend_list

######################################################################################

save.all.liks <- paste0(save.path, "BMOUi_",group.name,"Extinct_Miocene_ModelFit.csv")
    write.csv(all.models.lik, file=save.all.liks, row.names=F)
# unique(all.models.lik$V5); length(unique(all.models.lik$V5))
# filter(all.models.lik, V4==49)

# names(model.fit); length(names(model.fit))
save.all.models <- paste0(save.path, "BMOUi_", group.name,"Extinct_Miocene_ModelObjects.RDS")
    saveRDS(model.fit, file = save.all.models)
# model.fit <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Pygopodoidea_ModelObjects.RDS")   
    
    
total.results <- all.models.lik
# total.results <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/ExtantExtinct_Stochastic_ModelFit.csv")
colnames(total.results) <- c("lnL", "AIC", "AICc", "tree.no", "model", "shift.time", "group")
total.results$tree.no <- as.numeric(total.results$tree.no) # set the tree numbers as actual numbers
extinct <- dplyr::filter(total.results, group=="Extinct")
extant <- dplyr::filter(total.results, group=="Extant")

extant.first <- subset(extant, extant$tree.no < 10) # subset 1-9
  extinct.first <- subset(extinct, extinct$tree.no < 10) # subset 1-9

extant.first <- arrange(extant.first, tree.no) # sort this subset
  extinct.first <- arrange(extinct.first, tree.no)
  
extant.second <- subset(extant, extant$tree.no > 9) # subset 10-100
  extinct.second <- subset(extinct, extinct$tree.no > 9)

extant.second <- arrange(extant.second, tree.no) # sort this subset, this is done to avoid counting errors
  extinct.second <- arrange(extinct.second, tree.no)

sorted.extant <- rbind.data.frame(extant.first, extant.second) # now bind them together
  sorted.extinct <- rbind.data.frame(extinct.first, extinct.second) # now bind them together

weight.delta.extant <- NULL; weight.delta.extinct <- NULL
for (j in 1:100){
  targetdata1 <- dplyr::filter(sorted.extant, tree.no==j)
  res <- geiger::aicw(as.numeric(as.character(targetdata1$AICc)))
  weight.delta.extant <- rbind.data.frame(weight.delta.extant, res)
  
  targetdata2 <- dplyr::filter(sorted.extinct, tree.no==j)
  res <- geiger::aicw(as.numeric(as.character(targetdata2$AICc)))
  weight.delta.extinct <- rbind.data.frame(weight.delta.extinct, res)
}
weight.delta.extant <- within(weight.delta.extant, rm(fit))
  weight.delta.extinct <- within(weight.delta.extinct, rm(fit))

sorted.extant <- cbind.data.frame(sorted.extant, weight.delta.extant)
  sorted.extinct <- cbind.data.frame(sorted.extinct, weight.delta.extinct)

#sub <- subset(sorted.results, sorted.results$tree.no < 101)
extant.outz <- summarySE(sorted.extant, measurevar="w", groupvars="model")
extinct.outz <- summarySE(sorted.extinct, measurevar="w", groupvars="model")



## If you want to read in data to plot
########################################
input.res <- read.csv(file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/PlioPleistocene.SRC.results.csv")
extinct <- subset(input.res, stree.res$tree.type == "fossil tree")
extant  <- subset(input.res, stree.res$tree.type == "extant tree")
outz <- summarySE(extant, measurevar="weight.w", groupvars="model")
colnames(outz) <- c("model", "N", "w", "sd", "se", "ci")


## Otherwise Just Plot the model results
########################################
outz <- extant.outz
outz$model <- factor(outz$model, levels=c("BM1", "EB", "OU1", "BMOU", "BMOUi", "BMtrend"))

myplot <- (ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + theme(axis.text.x=element_text(angle=45, hjust=1))
  + scale_fill_manual(values=wes_palette("Royal2", 6, "continuous")))
extant.BMOUi.Miocene <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
multiplot(extinct.BMOUi.stochastic, extant.BMOUi.stochastic, 
          extinct.BMOUi.stochastic, extant.BMOUi.stochastic,
          extinct.BMOUi.Miocene,    extant.BMOUi.Miocene, cols=1)
# multiplot(extinct.SRC, extinct.low.BM, extinct.hi.BM, extinct.OU,
#           extant.SRC,  extant.low.BM, extant.hi.BM, extant.OU,
#           cols=2)
multiplot(extinct.hi.BM, extant.hi.BM, 
          extinct.hi.BM, extant.hi.BM,
          extinct.hi.BM, extant.hi.BM,
          extinct.hi.BM, extant.hi.BM,cols=2)

## Plot the composite bar graphs
#group.outz$model <- factor(group.outz$model, levels=c("BM", "delta", "kappa", "lambda", "gamma", "EB", "DensDep", "OU", "BMS", "TS", "OUS", "OUMA", "OUMV", "OUMVA", "SRC", "TRC")) # this re-orders the models in the legend
outz$model <- factor(outz$model, levels=c("BM", "EB", "OU", "SRC", "TRC")) # this re-orders the models in the legend

(ggplot(outz)
  + geom_bar(aes(y=w, x=model, fill=model), stat="identity")
  + theme(axis.text.x=element_text(angle=25, hjust=1), panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual( values=wes_palette("Zissou", 5, "continuous")))

# if you just want to do a single plot
myplot <- (ggplot(outz)
           + geom_bar(aes(1, y=w, fill=outz$model), stat="identity")
           + scale_fill_manual( values=wes_palette("Zissou", 12, "continuous")))

myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))






#####################################################################################
### Use PhyTools to plot trait differences between two trees (extinct vs. extant)
#####################################################################################
tree <- input.trees[[21]]
treedata <- save.sim.traits[[21]]
datanames <- setNames(treedata[,1], rownames(treedata))
plotTree.barplot(tree, datanames, list(col="red"))
extanttree <- drop.extinct(tree, tol=0.00001)
plotTree.barplot(extanttree, datanames)
par(mfrow=c(2,4))




record <- simFossilRecord(p=0.1, q=0.1, r=0, nruns=1, nTotalTaxa=50, plot=TRUE)


#####################################################################################
### Use PaleoTree to plot the differences in LTT/trees for the extinct/extant trees
#####################################################################################
phyloDiv(sim.extinct[[1]], int.length=0.1, plotLogRich=T)
phyloDiv(drop.extinct(sim.extinct[[1]], tol=0.00001), int.length=0.1, plotLogRich=T)
par(mfrow=c(2,1))

trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Miocene.Fossil.trees")
trees <- trees[c(1,21,41,61,81)]

#####################################################################################
### Use PaleoTree to plot the differences in LTT/trees for the extinct/extant trees
#####################################################################################
b <- 5
phyloDiv(trees[[b]], int.length=0.1, plotLogRich=T)
phyloDiv(drop.extinct(trees[[b]], tol=0.00001), int.length=0.1, plotLogRich=T)
par(mfrow=c(2,1))






#### If you need to find tips that died within a certain age.
test <- input.trees
targets <- list()
for (k in 1:length(test$tip.label)) {
  test <- input.trees[[14]] # designate the tree of interest
  mnh <- max(nodeHeights(test)) # Max Node Height
  tnh <- nodeheight(test, k) # Target Node Height
  if ((mnh - tnh < 11.5) & (mnh - tnh > 10.5)) { # set the time range you're interested in
    targets <- append(targets, k)
  }
}

tist <- NULL
for (p in 1:length(input.trees)) {
  test <- input.trees[[p]]
  for (k in 1:length(test$tip.label)) {
    mnh <- max(nodeHeights(test)) # Max Node Height of the tree
    tnh <- nodeheight(test, k) # Target Node Height (height above root of the node/tip of interest)
    if (mnh-tnh < 10.5 & mnh-tnh > 9.5) { # set the time range you're interested in
      tist <- rbind(tist, c(mnh-tnh, k, p))
    }
  }
}
colnames(tist) <- c("extinction.time", "tip.number", "tree.number")
tost <- as.data.frame(tist)
tist.out <- summarySE(tost, measurevar="extinction.time", groupvars="tree.number")
# 14,16,17,21,30,31,33,37,38,40,50,54,62,63,64,



# add.lengths to appropriate edges of the tree
res.counts <- table(shift.positions.list) # make a table of the shift frequencies
shifted.tips <- as.data.frame(res.counts) # turn it into a data frame
all.node.numbers <- as.data.frame(trees[[1]]$tip.label) # get all the nodes of the tree and the numbers
all.node.numbers[,"tip.no"] <- rownames(all.node.numbers); colnames(all.node.numbers) <- c("tip.name", "tip.no") # make a column that shows the tip number
target.numbers <-  all.node.numbers[all.node.numbers$tip.name %in% shifted.tips$shift.positions.list,] # subset all the tips, so show just the node numbers of the shifted tips
target.numbers[,2] <- as.numeric(target.numbers[,2]) # make it numeric
target.numbers <- target.numbers[order(match(target.numbers$tip.name, shifted.tips$shift.positions.list)),] # match the new frame to the output (shift) frame
target.numbers[,"shift.freq"] <- cbind(shifted.tips$Freq) # add on the frequencies of shifts


for (i in 1:length(target.numbers[,2])){
  rownames(tree$edge) <- c(1:length(tree$edge[,1])) # give the tree edge frame rownames
  target.almost <- subset(tree$edge, tree$edge[,2]==target.numbers[,2][i]) # pull out the ancestor and descendant nodes of the target edge
  interim.target <- subset(target.numbers, target.numbers[,2]==target.numbers[,2][i]) # subset the data frame to just the current tip of interest (descendant node)
  target.edge <- as.numeric(rownames(target.almost)) # get the number of the target edge
  tree$edge.length[[target.edge]] <- tree$edge.length[[target.edge]]+(0.01*interim.target[,3]) # add the desired length to the branch, per shift (here, 0.01)
}
shift.fit$tree <- tree # set the tree in the 'shift.fit' object to our rescaled tree
#plot(tree) # have a look to see if it worked
plot(shift.fit, cex=1)

tree$tip.label



