#install.packages("neldermead")
#install.packages("fBasics")
#install.packages("gsl")
#install.packages("stabledist")
#install.packages("statmod")
library("neldermead")
library("fBasics")
library("gsl")
library("stabledist")
library("statmod")
source("/Users/Ian/Google.Drive/R.Analyses/pulsR/levy_pruning_cf.r")
source("/Users/Ian/Google.Drive/R.Analyses/pulsR/levy_pruning_tools.r")
source("/Users/Ian/Google.Drive/R.Analyses/pulsR/levy_pruning_prob.r")
source("/Users/Ian/Google.Drive/R.Analyses/pulsR/levy_pruning_optim.r")
source("/Users/Ian/Google.Drive/R.Analyses/pulsR/levy_pruning_sim.r")
 # to source each script, you have to provide the path to their dependent source scripts (open each file and check)


#######################################################
# Read in all the data you'll use
######################################################
#tree<-read.nexus("BT.Pygopodoidea.tre") #read in desired tree
trees <- read.nexus("PB.Meliphagides.100.trees")
#drip <- c("S_Papuascincus_sp","S_Prasinohaema_virens","S_Sphenomorphus_jobiensis", "S_Sphenomorphus_muelleri","S_Sphenomorphus_solomonis")
#trees<-lapply(trees, drop.tip, tip=drip) #drop tips if necessary


#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data       <- read.csv("BT.Meliphagides.logMASS.csv", row.names = 1, header=F) #read in data file in GEIGER format
data.OUwie <- read.csv("BT.Meliphagides.logMASS.csv", header=F) #read in data file in OUwie format
data.fitenv<- data.OUwie[,2]; names(data.fitenv) <- data.OUwie[,1] #read in data file in RPANDA format

pulsr.bm <- fit_reml_levy(trees[[1]], data.fitenv, model="BM")
gebm <- fitContinuous(trees[[1]], data, model="BM")
trc <- fitContinuous_paleo(trees[[1]], data, model="radiate.constrain", shift.time=5)
trc2 <-fitContinuous_paleo(trees[[1]], data, model="TRC", shift.time=5)
pulsr.jn <- fit_reml_levy(trees[[1]], data.fitenv, model="JN")
pulsr.nig <- fit_reml_levy(trees[[1]], data.fitenv, model="NIG")
pulsr.bmjn <- fit_reml_levy(trees[[1]], data.fitenv, model="BMJN")
pulsr.bmnig <- fit_reml_levy(trees[[1]], data.fitenv, model="BMNIG")

pulsr.bm$AIC; gebm$opt$aic; trc$Trait1$aic; pulsr.jn$AIC; pulsr.nig$AIC

calc_AICc_vals(pulsr.nig$lnL, pulsr.nig$n_params, length(pulsr.nig$dat))
calc_AICc_vals(trc$Trait1$lnl, trc$Trait1$k, length(data[,1]))













