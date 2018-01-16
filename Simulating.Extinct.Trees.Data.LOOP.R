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

keep <- sim.extinct
t <- length(keep)
for (i in 1:length(sim.extinct)){
  keep[[t+i]] <- sim.extinct[[i]]
}

miocene.agam <- keep
miocene.mars <- keep
miocene.bird <- keep
miocene.pygo <- keep
miocene.skink <- keep

stochastic.agam <- keep
stochastic.mars <- keep
stochastic.bird <- keep
stochastic.pygo <- keep
stochastic.skink <- keep
## simulate trees under the "old extinction" concept
stochastic.trees <- list(); miocene.trees <- list()
stochastic.trees <- c(stochastic.agam, stochastic.mars, stochastic.bird, stochastic.pygo, stochastic.skink)
miocene.trees <- c(miocene.agam, miocene.mars, miocene.bird, miocene.pygo, miocene.skink)
class(miocene.trees) <- "multiPhylo"
write.tree(miocene.trees, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Miocene.Fossil.trees")

## Designate the trees you'd like to use for the data simulations and model fitting
###################################################################################
input.trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.PlioPleistocene.Fossil.number.trees")

## Set your parameters relevant to your empirical parameter estimates
######################################################
# alpha
#alpha = 2 # either a static value
diff.alpha <- seq(from=0.5, to=5, by=0.1) # or set it as a vector of sampled values
diff.alpha <- sample(diff.alpha, size=100, replace=T) # or set it as a vector of sampled values
# pre/post shift sigma
preshift.sigma  = 1
postshift.sigma = 5
# and shift time
sim.shifts <- seq(from=1, to=20, by=1)
sim.shifts <- sample(sim.shifts, size=100, replace=T)
######################################################
## Set empty object to hold your results
sim.traits.geiger <- list(); sim.traits.ouwie <- list() # make trait lists for all trees in geiger and ouwie data format
save.sim.traits <- NULL
extant.data <- list()
# make sure to the outputs depending on if you're simulating different times, or alphas
# this means changing the 'm' and the 'm[[2]]' objects below!
num.sims <- 1 # designate the number of simulated data sets you like to create per tree
### if trees are sensitive to shift dates: (this is for the Plio-Pleistocene trees)
######################################################
cheese.time <- c(10,10,10,11,10,10,10,10,10,10,
                 10,10,10,12,10,11,12,10,10,10,
                 12,12,10,11,10,11,10,10,10,11,
                 12,10,11,10,10,10,11,11,10,11,
                 10,10,10,10,10,10,11,10,10,11,
                 10,10,10,11,10,10,10,10,10,10,
                 10,11,11,11,10,10,10,10,10,10,
                 10,10,10,10,10,10,10,10,10,10,
                 10,10,10,10,10,10,10,10,10,10,
                 10,10,10,10,10,10,10,10,10,10)
######################################################
stree.res <- NULL

##########################################################################################################
## LOOP 1:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either the SRC ('Single-rate-constraint', BM-OU, 1 sig2) or the 
### TRC ('Two-rate-constraint', BM-OU, 2 sig2) model, then comparatively fit a set of standard models 
### (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a different model use the second loop, to simulate under a Brownian
### motion or Ornstein-Uhlenbeck process, using Diversitree (although you could also use Geiger).

### This loop will simulate data onto the tree with fossil tips, fit a series of models, then
### drop the extinct tips and associated data, refit the same models, and provide a summary for each
##########################################################################################################
for (z in 91:length(input.trees)) {
  traits.geiger <- list(); traits.ouwie <- list() # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- input.trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(input.trees), "\n") #keep track of what tree/loop# we're on
  
  for (i in 1:num.sims) {
    # Option A (comment out top when simulating different shift times, comment out bottom when simulating different alphas)
    m <- split.matrices <- split.vcv(phy, cheese.time[z]) # divide the vcv matrix at a static time
    #m <- split.matrices <- split.vcv(phy, sim.shifts[i]) # or at differing times
    
    # Option B (adjust to change simulating model, comment out both to get the SRC model)
    #m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
    #m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
    
    # Option C (comment out top when simulating different times, comment out bottom when simulating different alphas)
    m[[2]] <- ouMatrix(m[[2]], alpha=diff.alpha[z]) # transform the second half of the vcv matrix according to your alpha
    #m[[2]] <- ouMatrix(m[[2]], alpha=alpha) # transform the second half of the vcv matrix according to your alpha
    
    m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
    
    # OR, do it like the 'ecological release' model
    # m <- lapply(m, function(x) x*0.1)
    # m.rev.rel <- m[[1]] + m[[2]]
    
    # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
    traitz[[i]] <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), 
                                    sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))
    t.ouwie <- NULL
    t.ouwie <- as.data.frame(names(traitz[[i]]))
    t.ouwie[,2] <- as.data.frame(t(traitz[[i]]))
    #traitz.ouwie[[i]] <- t.ouwie
    
    t.geiger <- t.ouwie; names <- t.ouwie[,1]
    rownames(t.geiger) <- names; t.geiger[,1] <- NULL
    #traitz.geiger[[i]] <- t.geiger
    
    traits.geiger[[i]] <- t.geiger
    traits.ouwie[[i]]  <- t.ouwie
  }
  sim.traits.geiger[[z]] <- traits.geiger
  #sim.traits.ouwie[[z]] <- traits.ouwie
  
  
  save.sim.traits[[z]] <- as.data.frame(sim.traits.geiger[[z]]); 
  save.sim.traits[[z]][,"tree.num"] <- z
  save.sim.traits[[z]][,"gen.model"] <- "SRC"
  save(save.sim.traits, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Traits.PlioPleistocene.SRC.RData")
  
  for (i in 1:num.sims) {
    tree <- input.trees[[z]] # change this to match the tree size you want
    data <- sim.traits.geiger[[z]][[i]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[z]][[i]] # change this to match the proper sized tree
    
    bmfit    <- fitContinuous_paleo(tree, data, model="BM")
    ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
    oufit    <- fitContinuous(tree, data, SE=NA, model="OU", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=cheese.time[z])
    SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=cheese.time[z], bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    #OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
    #OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
    #OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
    
    
    #####################################################
    ###### Summarize and Compare Model Fitting ##########
    #####################################################
    results.names <- list(ebfit$opt, oufit$opt)
    results <- NULL
    for (k in 1:length(results.names)) {
      x <- as.data.frame(results.names[k])
      results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
    }
    
    results.nonstan <- NULL
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(bmfit$Trait1$lnl, bmfit$Trait1$aic, bmfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
    
    ## combine both
    results <- rbind.data.frame(results, results.nonstan)
    model <- c("EB", "OU", "BM", "TRC", "SRC") #, "OUM", "OUMA", "OUMV")
    tree.type <- paste("fossil tree")
    results[,"tree.type"] <- tree.type; results[,"tree.num"] <- z
    colnames(results) <- c("lnL", "AIC", "AICc", "tree.type", "tree.num")
    results <- cbind(results, model)
    
    ## Use AIC weights to determine best fitting model and model contributions
    weight <- aicw(results$AICc)
    results <- cbind(results, weight$delta)
    results <- cbind(results, weight$w)
    stree.res <- rbind.data.frame(stree.res, results)
  }
  
  drops <- is.extinct(phy, tol=0.0001) # find out which tips are extinct
  extant.sim.tree <- drop.extinct(phy, tol=0.00001)
  extant.data[[z]] <- subset(sim.traits.geiger[[z]][[i]], !rownames(sim.traits.geiger[[z]][[i]]) %in% drops)
  
  for (i in 1:num.sims) {
    tree <- extant.sim.tree # change this to match the tree size you want
    data <- extant.data[[z]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[z]][[i]] # change this to match the proper sized tree
    
    bmfit    <- fitContinuous_paleo(tree, data, model="BM")
    ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
    oufit    <- fitContinuous(tree, data, SE=NA, model="OU", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=cheese.time[z])
    SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=cheese.time[z], bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    #OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
    #OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
    #OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
    
    
    #####################################################
    ###### Summarize and Compare Model Fitting ##########
    #####################################################
    results.names <- list(ebfit$opt, oufit$opt)
    results <- NULL
    for (k in 1:length(results.names)) {
      x <- as.data.frame(results.names[k])
      results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
    }
    
    results.nonstan <- NULL
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(bmfit$Trait1$lnl, bmfit$Trait1$aic, bmfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
    
    ## combine both
    results <- rbind.data.frame(results, results.nonstan)
    model <- c("EB", "OU", "BM", "TRC", "SRC") #, "OUM", "OUMA", "OUMV")
    tree.type <- paste("extant tree")
    results[,"tree.type"] <- tree.type; results[,"tree.num"] <- z
    colnames(results) <- c("lnL", "AIC", "AICc", "tree.type", "tree.num")
    results <- cbind(results, model)
    
    ## Use AIC weights to determine best fitting model and model contributions
    weight <- aicw(results$AICc)
    results <- cbind(results, weight$delta)
    results <- cbind(results, weight$w)
    stree.res <- rbind.data.frame(stree.res, results)
  }
}
extinct <- subset(stree.res, stree.res$tree.type == "fossil tree")
extant  <- subset(stree.res, stree.res$tree.type == "extant tree")
outz <- summarySE(extant, measurevar="weight$w", groupvars="model")
colnames(outz) <- c("model", "N", "w", "sd", "se", "ci")

write.csv(stree.res, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/PlioPleistocene.SRC.results.csv")


##########################################################################################################
## LOOP 2:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either Brownian Motion (BM) or Ornstein-Uhlenbeck (OU) models then comparatively fit a set 
### of standard models (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a mode-variable model, use the first loop.
##########################################################################################################

bm.pars <- 0.1 # set the diffusion parameter of the BM process
ou.pars <- c(0.5, sample(diff.alpha, 1), 1) # set the diffusion parameter, the alpha, and the optimum

for (z in 13:length(input.trees)) {
  traits.geiger <- list(); traits.ouwie <- list() # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- input.trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(input.trees), "\n") #keep track of what tree/loop# we're on
  
  for (i in 1:num.sims) {
    simulated.traits <- NULL
    #simulated.traits <- as.data.frame(sim.character(phy, model="bm", bm.pars))
    simulated.traits <- as.data.frame(sim.character(phy, model="ou", ou.pars))
    traits.geiger[[i]] <- simulated.traits
  }
  
  sim.traits.geiger[[z]] <- traits.geiger
  #sim.traits.ouwie[[z]] <- traits.ouwie
  
  save.sim.traits[[z]] <- as.data.frame(sim.traits.geiger[[z]]); 
  save.sim.traits[[z]][,"tree.num"] <- z
  save.sim.traits[[z]][,"gen.model"] <- "OU" # OU, lowBM, hiBM, SRC
  save(save.sim.traits, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Traits.PlioPleistocene.OU.RData")
  
  for (i in 1:num.sims) {
    tree <- input.trees[[z]] # change this to match the tree size you want
    data <- sim.traits.geiger[[z]][[i]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[z]][[i]] # change this to match the proper sized tree
    
    bmfit    <- fitContinuous_paleo(tree, data, model="BM")
    ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
    oufit    <- fitContinuous(tree, data, SE=NA, model="OU", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=cheese.time[z])
    SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=cheese.time[z], bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    #OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
    #OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
    #OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
    
    
    #####################################################
    ###### Summarize and Compare Model Fitting ##########
    #####################################################
    results.names <- list(ebfit$opt, oufit$opt)
    results <- NULL
    for (k in 1:length(results.names)) {
      x <- as.data.frame(results.names[k])
      results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
    }
    
    results.nonstan <- NULL
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(bmfit$Trait1$lnl, bmfit$Trait1$aic, bmfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
    
    ## combine both
    results <- rbind.data.frame(results, results.nonstan)
    model <- c("EB", "OU", "BM", "TRC", "SRC") #, "OUM", "OUMA", "OUMV")
    tree.type <- paste("fossil tree")
    results[,"tree.type"] <- tree.type; results[,"tree.num"] <- z
    colnames(results) <- c("lnL", "AIC", "AICc", "tree.type", "tree.num")
    results <- cbind(results, model)
    
    ## Use AIC weights to determine best fitting model and model contributions
    weight <- aicw(results$AICc)
    results <- cbind(results, weight$delta)
    results <- cbind(results, weight$w)
    stree.res <- rbind.data.frame(stree.res, results)
  }
  
  drops <- is.extinct(phy, tol=0.0001) # find out which tips are extinct
  extant.sim.tree <- drop.extinct(phy, tol=0.00001)
  extant.data[[z]] <- subset(sim.traits.geiger[[z]][[i]], !rownames(sim.traits.geiger[[z]][[i]]) %in% drops)
  
  for (i in 1:num.sims) {
    tree <- extant.sim.tree # change this to match the tree size you want
    data <- extant.data[[z]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[z]][[i]] # change this to match the proper sized tree
    
    bmfit    <- fitContinuous_paleo(tree, data, model="BM")
    ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
    oufit    <- fitContinuous(tree, data, SE=NA, model="OU", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=cheese.time[z])
    SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=cheese.time[z], bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
    #OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
    #OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
    #OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
    
    
    #####################################################
    ###### Summarize and Compare Model Fitting ##########
    #####################################################
    results.names <- list(ebfit$opt, oufit$opt)
    results <- NULL
    for (k in 1:length(results.names)) {
      x <- as.data.frame(results.names[k])
      results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
    }
    
    results.nonstan <- NULL
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(bmfit$Trait1$lnl, bmfit$Trait1$aic, bmfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
    #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
    
    ## combine both
    results <- rbind.data.frame(results, results.nonstan)
    model <- c("EB", "OU","BM", "TRC", "SRC") #, "OUM", "OUMA", "OUMV")
    tree.type <- paste("extant tree")
    results[,"tree.type"] <- tree.type; results[,"tree.num"] <- z
    colnames(results) <- c("lnL", "AIC", "AICc", "tree.type", "tree.num")
    results <- cbind(results, model)
    
    ## Use AIC weights to determine best fitting model and model contributions
    weight <- aicw(results$AICc)
    results <- cbind(results, weight$delta)
    results <- cbind(results, weight$w)
    stree.res <- rbind.data.frame(stree.res, results)
  }
}
extinct <- subset(stree.res, stree.res$tree.type == "fossil tree")
extant  <- subset(stree.res, stree.res$tree.type == "extant tree")
outz <- summarySE(extant, measurevar="weight$w", groupvars="model")
colnames(outz) <- c("model", "N", "w", "sd", "se", "ci")

write.csv(stree.res, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/PlioPleistocene.OU.results.csv")



## If you want to read in data to plot
########################################
input.res <- read.csv(file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/PlioPleistocene.SRC.results.csv")
extinct <- subset(input.res, stree.res$tree.type == "fossil tree")
extant  <- subset(input.res, stree.res$tree.type == "extant tree")
outz <- summarySE(extant, measurevar="weight.w", groupvars="model")
colnames(outz) <- c("model", "N", "w", "sd", "se", "ci")


## Otherwise Just Plot the model results
########################################
myplot <- (ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + theme(axis.text.x=element_text(angle=45, hjust=1))
  + scale_fill_manual(values=wes_palette("Royal2", 5, "continuous")))
extant.SRC <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
multiplot(extinct.SRC, extinct.low.BM, extinct.hi.BM, extinct.OU,
          extant.SRC,  extant.low.BM, extant.hi.BM, extant.OU,
          cols=2)
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



