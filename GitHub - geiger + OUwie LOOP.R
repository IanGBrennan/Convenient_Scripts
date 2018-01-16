library(phytools)
library(geiger)
library(mvtnorm)
library(ggplot2)
library(OUwie)
library(plyr)
library(dplyr)
library(wesanderson)
library(Rmisc)
library(phylobase)
library(RPANDA)
# *warning 'calc_AIC_vals' is sourced from somewhere in BioGeoBEARS, I should find out where


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


#logBL<-setNames(data[,14], rownames(data)) #choose your data column (here: logBL) and apply rownames
name.check(trees[[3]], data); name.check(trees[[3]], data.fitenv) #check to make sure the tips match the data labels

total.results <- NULL
tree.model.comparison <- NULL
all.shift.model.estimates <- NULL
for (i in 1:length(trees)) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  ptm <- proc.time() # get an idea of how long each loop takes
  
  ################################################################
  # run Continuous model fitting on the data using standard models
  ################################################################
  bmfit<-fitContinuous(trees[[i]], data, SE=NA, model="BM", control=(niter=1000))
  ebfit<-fitContinuous(trees[[i]], data, SE=NA, model="EB", control=(niter=1000))
  deltafit<-fitContinuous(trees[[i]], data, SE=NA, model="delta", control=(niter=1000))
  kappafit<-fitContinuous(trees[[i]], data, SE=NA, model="kappa", control=(niter=1000))
  oufit<-fitContinuous(trees[[i]], data, SE=NA, model="OU", control=(niter=1000))
  lambdafit<-fitContinuous(trees[[i]], data, SE=NA, model="lambda", control=(niter=1000))
  
  ################################################################
  # run a model which accounts for Diversity Dependence (decline)
  ################################################################
  DDfit <- fitDiversityModel(trees[[i]], data, showTree=F) # fit the model
  DDloglik <- DDfit$logL # set the loglikelihood as an object
  DDAIC <- calc_AIC_vals(DDfit$logL, 2) # and estimate the AIC score
  DDAICc <- calc_AICc_vals(DDfit$logL, 2, nTips(trees[[i]])) # estimate the AICc ('calc_AICc_vals' is from BioGeoBEARS)
  DDres <- NULL
  DDres <- rbind(DDres, c(DDloglik, DDAIC, DDAICc)); colnames(DDres) <- c("lnL", "AIC", "AICc"); rownames(DDres)<-"DensDep"

  ######################################################
  ######## Slater Paleo Code from 2013 Mammal Paper ####
  ######################################################
  #the fxn 'fitContinuous_paleo' allows timed shifts and release of trait evolution
  source("https://github.com/IanGBrennan/MioceneAustralia/blob/master/Scripts/Phenotypic_Evolution/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your GitHub
  
  res.ts <- NULL # create a null data frame for the GEIGER time-shift model (BM)
  res.trc <- NULL # create a null data frame for the GEIGER Two-Rate Constraint model (BM+OU)
  res.src <- NULL # create a null data frame for the  GEIGER Single-Rate Constraint model (BM+OU)
  res.bms <- NULL # create a null data frame for the OUwie time-shift model (BM)
  res.ous <- NULL # create a null data frame for the OUwie multiple optima (OUM)
  res.ouma <- NULL # create a null data frame for the OUwie multiple attractions model (OUMA)
  res.oumv <- NULL # create a null data frame for the OUwie multiple variance model (OUMV)
  res.oumva <- NULL # createa a null data frame for the OUwie multiple attractions/variance model (OUMVA)
  
  best.ts <- NULL # create a null data frame for the BEST GEIGER time-shift model (BM)
  best.trc <- NULL # create a null data frame for the BEST GEIGER Two-Rate Constraint model (BM+OU)
  best.src <- NULL # create a null data frame for the BEST GEIGER Single-Rate Constraint model (BM+OU)
  best.bms <- NULL # create a null data frame for the BEST OUwie time-shift model (BM)
  best.ous <- NULL # create a null data frame for the BEST OUwie multiple optima (OU)
  best.ouma <- NULL # create a null data frame for the BEST OUwie multiple attractions model (OUMA)
  best.oumv <- NULL # create a null data frame for the BEST OUwie multiple variance model (OUMV)
  best.oumva <- NULL # createa a null data frame for the BEST OUwie multiple attractions/variance model (OUMVA)

  ## Fit a BM time shift model (this includes saving parameter estimates for each iteration)
  TS <- NULL
  for (h in 3:11) {
    time.shift<-fitContinuous_paleo(trees[[i]], data, model="timeshift", shift.time=h)
    res.ts <- rbind(res.ts, as.data.frame(t(c(time.shift$Trait1$lnl, time.shift$Trait1$aic, time.shift$Trait1$aicc))))
    ts.estimates <- NULL
    ts.estimates <- rbind(ts.estimates, as.data.frame(t(c(time.shift$Trait1$beta, time.shift$Trait$rate2))))
    ts.estimates <- rbind(ts.estimates, c(NA, NA))
    ts.estimates <- rbind(ts.estimates, c(time.shift$Trait1$root.state, NA))
    rownames(ts.estimates) <- c("sigma.sq", "alpha", "theta")
    ts.estimates [,3] <-paste("TS",h)
    TS <- rbind(TS, ts.estimates)
  }
  colnames(TS) <- c("before.shift", "after.shift", "model")
  
  ## Fit the Two Rate Constrain (BM/OU) model (this includes saving parameter estimates for each iteration)
  TRC <- NULL
  for (h in 3:11) {
    TRCfit<-fitContinuous_paleo(trees[[i]], data, model="radiate.constrain", shift.time=h)
    res.trc <- rbind(res.trc, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    TRC.estimates <- NULL
    TRC.estimates <- rbind(TRC.estimates, as.data.frame(t(c(TRCfit$Trait1$beta, TRCfit$Trait1$post.shift.scalar))))
    TRC.estimates <- rbind(TRC.estimates, c(NA, TRCfit$Trait1$alpha))
    TRC.estimates <- rbind(TRC.estimates, c(time.shift$Trait1$root.state, NA))
    rownames(TRC.estimates) <- c("sigma.sq", "alpha", "theta")
    TRC.estimates [,3] <-paste("TRC",h)
    TRC <- rbind(TRC, TRC.estimates)
  }
  colnames(TRC) <- c("before.shift", "after.shift", "model")
  
  ## Fit the Single Rate Constrain (BM/OU) model (this includes saving parameter estimates for each iteration)
  SRC <- NULL
  for (h in 3:11) {
    SRCfit<-fitContinuous_paleo(trees[[i]], data, model="SRC", shift.time=h)
    res.src <- rbind(res.src, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
    SRC.estimates <- NULL
    SRC.estimates <- rbind(SRC.estimates, as.data.frame(t(c(SRCfit$Trait1$beta, NA))))
    SRC.estimates <- rbind(SRC.estimates, c(NA, SRCfit$Trait1$alpha))
    SRC.estimates <- rbind(SRC.estimates, c(time.shift$Trait1$root.state, NA))
    rownames(SRC.estimates) <- c("sigma.sq", "alpha", "theta")
    SRC.estimates [,3] <-paste("SRC",h)
    SRC <- rbind(SRC, SRC.estimates)
  }
  colnames(SRC) <- c("before.shift", "after.shift", "model")
  
  ## Fit a BMS model with an unknown timeslice: 
  BMSfit <- OUwie.slice(trees[[i]],data.OUwie[,c(1,2)],model=c("BMS"), root.station=T, timeslices=c(NA))
  res.bms <- rbind(res.bms, as.data.frame(t(c(BMSfit$loglik, BMSfit$AIC, BMSfit$AICc))))
  obj <- as.data.frame(BMSfit["timeslices"]) # pull out the time of the shift
  BMS.shift <- max(nodeHeights(trees[[i]])) - obj[2,] # create an object with the shift time

  ## Fit a multiple optima OU model (OUM)
  OUMfit <- OUwie.slice(trees[[i]],data.OUwie[,c(1,2)],model=c("OUM"), root.station=T, timeslices=c(NA))
  res.ous <- rbind(res.ous, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
  obj <- as.data.frame(OUMfit["timeslices"]) # pull out the time of the shift
  OUM.shift <- max(nodeHeights(trees[[i]])) - obj[2,] # create an object with the shift time
  
  ## Fit a multiple attractions OU model (OUMA)
  OUMAfit <- OUwie.slice(trees[[i]],data.OUwie[,c(1,2)],model=c("OUMA"), root.station=T, timeslices=c(NA))
  res.ouma <- rbind(res.ouma, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
  obj <- as.data.frame(OUMAfit["timeslices"]) # pull out the time of the shift
  OUMA.shift <- max(nodeHeights(trees[[i]])) - obj[2,] # create an object with the shift time
  
  ## Fit a multiple variance OU model (OUMV)
  OUMVfit <- OUwie.slice(trees[[i]],data.OUwie[,c(1,2)],model=c("OUMV"), root.station=T, timeslices=c(NA))
  res.oumv <- rbind(res.oumv, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
  obj <- as.data.frame(OUMVfit["timeslices"]) # pull out the time of the shift
  OUMV.shift <- max(nodeHeights(trees[[i]])) - obj[2,] # create an object with the shift time
  
  ## Fit a multiple attractions/variance OU model (OUMVA)
  OUMVAfit <- OUwie.slice(trees[[i]],data.OUwie[,c(1,2)],model=c("OUMVA"), root.station=T, timeslices=c(NA))
  res.oumva <- rbind(res.oumva, as.data.frame(t(c(OUMVAfit$loglik, OUMVAfit$AIC, OUMVAfit$AICc))))
  obj <- as.data.frame(OUMVAfit["timeslices"]) # pull out the time of the shift
  OUMVA.shift <- max(nodeHeights(trees[[i]])) - obj[2,] # create an object with the shift time

  ## Set the column names of the results data frames
  colnames(res.ts)   <- c("lnL", "AIC", "AICc")
  colnames(res.trc)  <- c("lnL", "AIC", "AICc")
  colnames(res.src)  <- c("lnL", "AIC", "AICc")
  colnames(res.bms)  <- c("lnL", "AIC", "AICc")
  colnames(res.ous)  <- c("lnL", "AIC", "AICc")
  colnames(res.ouma) <- c("lnL", "AIC", "AICc")
  colnames(res.oumv) <- c("lnL", "AIC", "AICc")
  colnames(res.oumva)<- c("lnL", "AIC", "AICc")
  
  ## Set the row names of the results data frames
  rownames(res.ts) <- c(paste("TS",3:11))
  rownames(res.trc) <- c(paste("TRC",3:11))
  rownames(res.src) <- c(paste("SRC",3:11))
  rownames(res.bms) <- c(paste("BMS",BMS.shift))
  rownames(res.ous) <- c(paste("OUS",OUM.shift))
  rownames(res.ouma) <- c(paste("OUMA",OUMA.shift))
  rownames(res.oumv) <- c(paste("OUMV",OUMV.shift))
  rownames(res.oumva) <- c(paste("OUMVA",OUMVA.shift))
  
  ## Choose the best (lowest AIC) shift time from each model results data frame
  best.ts <- subset(res.ts, AICc == min(res.ts$AICc))
  best.trc <- subset(res.trc, AICc == min(res.trc$AICc))
  best.src <- subset(res.src, AICc == min(res.src$AICc))
  best.bms <- subset(res.bms, AICc == min(res.bms$AICc))
  best.ous <- subset(res.ous, AICc == min(res.ous$AICc))
  best.ouma <- subset(res.ouma, AICc == min(res.ouma$AICc))
  best.oumv <- subset(res.oumv, AICc == min(res.oumv$AICc))
  best.oumva <- subset(res.oumva, AICc == min(res.oumva$AICc))
  
  
  #####################################################
  ###### Summarize and Compare Model Fitting ##########
  #####################################################
  
  results <- NULL
  results.names <- list(bmfit$opt, ebfit$opt, oufit$opt,
                        kappafit$opt, deltafit$opt, lambdafit$opt)
  #non.standard <- list(time.shift$Trait1, TRCfit$Trait1, SRCfit$Trait1)
  for (k in 1:length(results.names)) {
    x <- as.data.frame(results.names[k])
    results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
  }
  rownames(results) <- c("BM", "EB", "OU", "kappa", "delta", "lambda")
  colnames(results) <- c("lnL", "AIC", "AICc")
  
  ## Add the density dependent model
  results <- rbind(results, DDres)
  
  ## Add the best non-standard GEIGER and OUwie model results to the results data frame
  results <- rbind(results, best.ts)
  results <- rbind(results, best.trc)
  results <- rbind(results, best.src)
  results <- rbind(results, best.bms)
  results <- rbind(results, best.ous)
  results <- rbind(results, best.ouma)
  results <- rbind(results, best.oumv)
  results <- rbind(results, best.oumva)
  
  ## Use AIC weights to determine best fitting model and model contributions
  weight <- aicw(results$AICc)
  results <- cbind(results, weight$delta)
  results <- cbind(results, weight$w)

  ## Make a data frame that holds ALL of the results
  total.results <- rbind(total.results, results)
  
  ## Make a data frame that holds only the best model(s) (and its AICc score) for each tree
  ### This includes any alternative models with delta values <2
  #score <- subset(results, AICc == min(results$AICc)) # this chooses only the best (lowest AICc) model
  score <- subset(results, weight$delta < 2) # this chooses all models with delta AICc < 2
  score["tree.no"] <- i
  tree.model.comparison <- rbind(tree.model.comparison, score)
  
  
  ##################################################################
  ###### Summarize and Compare Rate Shift Model Estimates ##########
  ##################################################################
  
  BM.box<-NULL
  BM.estimates <- as.data.frame(BMSfit["solution"])
  theta <- as.data.frame(BMSfit["theta"])
  BM.box<-cbind(BM.box,theta[1,1])
  BM.box<-cbind(BM.box,theta[2,1])
  BM.estimates<-rbind(BM.estimates, BM.box[1,])
  rownames(BM.estimates)<-c("alpha", "sigma.sq", "theta")
  BM.estimates[,3] <- "BMS"
  
  OUM.box<-NULL
  OUM.estimates <- as.data.frame(OUMfit["solution"])
  theta <- as.data.frame(OUMfit["theta"])
  OUM.box<-cbind(OUM.box,theta[1,1])
  OUM.box<-cbind(OUM.box,theta[2,1])
  OUM.estimates<-rbind(OUM.estimates, OUM.box[1,])
  rownames(OUM.estimates)<-c("alpha", "sigma.sq", "theta")
  OUM.estimates[,3] <- "OUM"
  
  OUMA.box<-NULL
  OUMA.estimates <- as.data.frame(OUMAfit["solution"])
  theta <- as.data.frame(OUMAfit["theta"])
  OUMA.box<-cbind(OUMA.box,theta[1,1])
  OUMA.box<-cbind(OUMA.box,theta[2,1])
  OUMA.estimates<-rbind(OUMA.estimates, OUMA.box[1,])
  rownames(OUMA.estimates)<-c("alpha", "sigma.sq", "theta")
  OUMA.estimates[,3] <- "OUMA"
  
  OUMV.box<-NULL
  OUMV.estimates <- as.data.frame(OUMVfit["solution"])
  theta <- as.data.frame(OUMVfit["theta"])
  OUMV.box<-cbind(OUMV.box,theta[1,1])
  OUMV.box<-cbind(OUMV.box,theta[2,1])
  OUMV.estimates<-rbind(OUMV.estimates, OUMV.box[1,])
  rownames(OUMV.estimates)<-c("alpha", "sigma.sq", "theta")
  OUMV.estimates[,3] <- "OUMV"
  
  OUMVA.box<-NULL
  OUMVA.estimates <- as.data.frame(OUMVAfit["solution"])
  theta <- as.data.frame(OUMVAfit["theta"])
  OUMVA.box<-cbind(OUMVA.box,theta[1,1])
  OUMVA.box<-cbind(OUMVA.box,theta[2,1])
  OUMVA.estimates<-rbind(OUMVA.estimates, OUMVA.box[1,])
  rownames(OUMVA.estimates)<-c("alpha", "sigma.sq", "theta")
  OUMVA.estimates[,3] <- "OUMVA"
  
  OU.model.estimates <- NULL
  OU.model.estimates <- rbind(OU.model.estimates, 
                              BM.estimates,
                              BM.estimates,
                              OUM.estimates,
                              OUMA.estimates,
                              OUMV.estimates,
                              OUMVA.estimates)
  colnames(OU.model.estimates) <- c("before.shift", "after.shift", "model")
  
  best.geiger <- NULL
  best.geiger <- rbind(best.geiger, subset(TS, TS["model"] == rownames(best.ts)))
  best.geiger <- rbind(best.geiger, subset(TRC, TRC["model"] == rownames(best.trc)))
  best.geiger <- rbind(best.geiger, subset(SRC, SRC["model"] == rownames(best.src)))
  
  all.est <- NULL
  all.est <- rbind(OU.model.estimates, best.geiger)
  all.est["tree.no"] <- i
  all.shift.model.estimates <- rbind(all.shift.model.estimates, all.est)

  print(proc.time() - ptm)
}

write.csv(tree.model.comparison, 
          file="/YOUR_DIRECTORY/Meliphagides.Model.Comparison.BEST.csv",
          quote=F)
write.csv(total.results,
          file="/YOUR_DIRECTORY/Meliphagides.Model.Comparison.FINALTOTAL.csv",
          quote=F)
write.csv(all.shift.model.estimates,
          file="/YOUR_DIRECTORY/Meliphagides.Shift.Model.FINALPARAM.ESTIMATES.csv",
          quote=F)


#####################################################################
## If you need to add an additional model...recalculate the AICcWt ##
## if you don't, then skip this step! ##
#####################################################################
# you need to run a loop of the model across all the trees, then add
# the AIC scores to the previous model.TOTAL and model.BEST files
# I'll run the Density Dependent ones for this example

#total.results <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Comparison TOTAL/Agamidae.Model.Comparison.TOTAL.csv", header=T)
total.results <- read.csv("/YOUR_DIRECTORY/FINAL.Meliphagides.Model.Comparison.TOTAL.csv", header=T)
total.results <- within(total.results, rm(X, delta, w)) #drop these 

## run a model which accounts for Diversity Dependence (decline)
for (z in 1:100) {
  DDfit <- fitDiversityModel(trees[[z]], data, showTree=F) # fit the model
  DDloglik <- DDfit$logL # set the loglikelihood as an object
  DDAIC <- calc_AIC_vals(DDfit$logL, 2) # and estimate the AIC score
  DDAICc <- calc_AICc_vals(DDfit$logL, 2, nTips(trees[[i]])) # estimate the AICc ('calc_AICc_vals' is from BioGeoBEARS)
  DDres <- NULL
  model <- paste("DensDep")
  treenum <- z
  timing <- NA
  DDres <- rbind(DDres, c(model, DDloglik, DDAIC, DDAICc, timing, treenum))
  colnames(DDres) <- c("model", "lnL", "AIC", "AICc", "timing", "tree.no")
  total.results <- rbind(total.results, DDres)
}

## run a model which tests the relation between trait evolution and an environmental variable (Clavel & Morlon, 2017)
data(InfTemp)
for (z in 1:100) {
  env.fit <- fit_t_env(trees[[z]], data.fitenv, env_data=InfTemp, df=50, scale=T, plot=F) # change the 'env_data' field if providing your data!
  env.loglik <- env.fit$LH; env.AIC <- env.fit$aic; env.AICc <- env.fit$aicc # set the loglikelihood as an object
  env.res <- NULL
  model <- paste("fit_env")
  treenum <- z
  timing <- NA
  env.res <- rbind(env.res, c(model, env.loglik, env.AIC, env.AICc, timing, treenum))
  colnames(env.res) <- c("model", "lnL", "AIC", "AICc", "timing", "tree.no")
  total.results <- rbind(total.results, env.res)
}

## If you want to drop some models before you determine AICWt [check with 'levels(total.results$model)']
total.results <- total.results[!(total.results$model=="kappa"),]
total.results <- total.results[!(total.results$model=="BM"),]
total.results <- total.results[!(total.results$model=="lambda"),]
total.results <- total.results[!(total.results$model=="BMS"),] # check it isn't "BMS " or "bms " or something similar
total.results <- total.results[!(total.results$model=="OU"),]

## Use AIC weights to determine best fitting model and model contributions
total.results$tree.no <- as.numeric(total.results$tree.no) # set the tree numbers as actual numbers
sorted.first <- subset(total.results, total.results$tree.no < 10) # subset 1-9
sorted.first <- arrange(sorted.first, tree.no) # sort this subset
sorted.second <- subset(total.results, total.results$tree.no > 9) # subset 10-100
sorted.second <- arrange(sorted.second, tree.no) # sort this subset, this is done to avoid counting errors
sorted.results <- rbind.data.frame(sorted.first, sorted.second) # now bind them together
weight.delta <- NULL
for (j in 1:100){
  targetdata <- subset(sorted.results, sorted.results$tree.no == j)
  res <- aicw(as.numeric(targetdata$AICc))
  weight.delta <- rbind.data.frame(weight.delta, res)
}
weight.delta <- within(weight.delta, rm(fit))
sorted.results <- cbind.data.frame(sorted.results, weight.delta)

write.csv(sorted.results,
          file="/YOUR_DIRECTORY/Trimmed.FINAL.Meliphagides.Model.Comparison.TOTAL.csv",
          quote=F)
#####################################################################



#############
# open the file in excel ("____.Model.Comparison.TOTAL.csv"):
# change first column name to "model", & remove numbers
# and last two columns to "delta" and "w"
#############
output.total <- read.csv("/YOUR_DIRECTORY/Model.Comparison TOTAL/Trimmed.FINAL.Meliphagides.Model.Comparison.TOTAL.csv", header=T)
output.total[,1] <- NULL # drop first column if it's 'X'
outz <- summarySE(output.total, measurevar="w", groupvars="model")

## If you want to plot a composite bar graph, subset each radiation, summarize model weights
outz.mars <- summarySE(subset(output.total, group=="marsupials"), measurevar="w", groupvars="model")
outz.agam <- summarySE(subset(output.total, group=="agamidae"), measurevar="w", groupvars="model")
outz.pygo <- summarySE(subset(output.total, group=="pygopodoidea"), measurevar="w", groupvars="model")
outz.skink <- summarySE(subset(output.total, group=="skinks"), measurevar="w", groupvars="model")
outz.bird <- summarySE(subset(output.total, group=="birds"), measurevar="w", groupvars="model")
## Then combine them back together
outz.mars["group"] <- "marsupial mammals"
outz.agam["group"] <- "agamidae lizards"
outz.pygo["group"] <- "pygopodoid geckos"
outz.skink["group"] <- "sphenomorphine skinks"
outz.bird["group"] <- "meliphagid birds"
group.outz <- rbind.data.frame(outz.mars, outz.agam, outz.pygo, outz.skink, outz.bird)

## If you want to summarize each model into a table (you probably don't)
###################################################
BMs <- subset(output.total, model == "BM")
EBs <- subset(output.total, model == "EB")
Deltas <- subset(output.total, model == "delta")
Kappas <- subset(output.total, model == "kappa")
OUs <- subset(output.total, model == "OU")
Lambdas <- subset(output.total, model == "lambda")
TSs <- subset(output.total, model == "ts")
TRCs <- subset(output.total, model == "trc")
SRCs <- subset(output.total, model == "src")
##################################################

## Otherwise Just Plot the model results
########################################
outz$model <- factor(outz$model, 
               levels=c("delta", "DensDep", "EB", "fit_env", "fit_bgb", 
                        "TS", "OUS", "OUMA", "OUMV", "OUMVA", 
                        "SRC", "TRC")) # re-order the models in the plot
(ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + theme(axis.text.x=element_text(angle=45, hjust=1)))

## Plot the composite bar graphs
#group.outz$model <- factor(group.outz$model, levels=c("BM", "delta", "kappa", "lambda", "gamma", "EB", "DensDep", "OU", "BMS", "TS", "OUS", "OUMA", "OUMV", "OUMVA", "SRC", "TRC")) # this re-orders the models in the legend
group.outz$model <- factor(group.outz$model, levels=c("delta","DensDep", "EB", "fit_env", "fit_bgb", "TS", "OUS", "OUMA", "OUMV", "OUMVA", "SRC", "TRC")) # this re-orders the models in the legend

(ggplot(group.outz)
  + geom_bar(aes(y=w, x=group, fill=group.outz$model), stat="identity")
  + theme(axis.text.x=element_text(angle=25, hjust=1), panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual( values=wes_palette("Zissou", 11, "continuous")))

# if you just want to do a single plot
myplot <- (ggplot(outz)
  + geom_bar(aes(1, y=w, fill=outz$model), stat="identity")
  + scale_fill_manual( values=wes_palette("Zissou", 12, "continuous")))

myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))


#######################################################################
###### Summarize and Compare ALL MODEL Parameter Estimates ############ 
######           to get the BEST MODEL Parameter Estimates ############
#######################################################################
####################################################
# open the file in excel ("_____.Shift.Model.PARAM.ESTIMATES.csv"):
# change first column name to "parameter", & remove numbers
# copy "model" column to end, rename "model.timing", remove numbers from "model" column
####################################################
est.total <- read.csv("/YOUR_DIRECTORY/Meliphagides.Shift.Model.FINALPARAM.ESTIMATES.csv", header=T)
#param.estimates <- summarySE(est.total, measurevar="w", groupvars="parameter")

####################################################
# open the file in excel ("_____.Model.Comparison.csv"):
# change first column name to "model", & copy it to last column, rename it "timing" 
# remove numbers from "model" and remove letters from "timing" (replace no-shift with NAs)
# make sure to remove extra numbers from the "timing" column, ie remove duplicates
# change names of "delta" and "w" columns
####################################################
best.models <- read.csv("/YOUR_DIRECTORY/Meliphagides.Model.Comparison.BEST.csv", header=T)

# (this might not be necessary) To make this all work, we have to combine levels of the two data frames
levels1 <- levels(best.models$model)
levels2 <- append(levels1, levels(est.total$model))
best.models$model <- factor(best.models$model, levels=levels2) # ignore duplicated levels warning
est.total$model <- factor(est.total$model, levels=levels2) # ignore duplicated levels warning

best.parameter.estimates <- NULL
for (i in 1:max(best.models$tree.no)){
  match.tree <- NULL; match.model <- NULL; best.batch.models <- NULL # make our empty data frames
  match.tree <- subset(est.total, i == est.total$tree.no) # choose only models applied to the 'i'th tree
  best.batch.models <- subset(best.models, best.models$tree.no == i)
  #match.model <- subset(match.tree, best.batch.models$model == est.total$model) # choose only the best fitting model from the above subset
  for (h in 1:nrow(best.batch.models)) {
    best.try <- NULL
    best.try <- subset(match.tree, best.batch.models$model[h] == match.tree$model) # choose only the best fitting model from the above subset
    match.model <- rbind(match.model, best.try)
  }
  #match.model <- match.model[complete.cases(match.model[,5]),] # remove any NAs introduced
  best.parameter.estimates <- rbind(best.parameter.estimates, match.model) # write only the model estimates of the BEST models
} # ignore the 50 duplicated levels warnings

write.csv(best.parameter.estimates,
          file="/YOUR_DIRECTORY/Meliphagides.BEST.Model.PARAM.ESTIMATES.csv",
          quote=F, row.names = F)



####################################################
## Plot the model parameter estimates
####################################################
####################################################
# open the file in excel ("____.BEST.Model.PARAM.ESTIMATES.csv"):
# take "model" column & copy it to last column, rename it "timing" #don't think you need to do this anymore
# remove numbers from "model" and remove letters from "timing" #don't think you need to do this anymore
####################################################
input <- read.csv("/YOUR_DIRECTORY/Meliphagides.BEST.Model.PARAM.ESTIMATES.csv", header=T)
#before <- summarySE(input, measurevar="before.shift", groupvars="parameter")
#after  <- summarySE(input, measurevar="after.shift",  groupvars="parameter")

## Subset each parameter into a data frame for plotting
trc.subset <- subset(input, input$model == "TRC")
trc.alpha <- subset(trc.subset, trc.subset$parameter == "alpha")
trc.sigma  <- subset(trc.subset, trc.subset$parameter == "sigma.sq")
trc.theta <- subset(trc.subset, trc.subset$parameter == "theta")

## Pull the mean and 95% range of the Post Shift Scalar
t.test(trc.sigma$after.shift)


## (DUMB) 
#######################
## This is my workaround for fixing the sigma estimates, it must be implemented!
## TRC model gives sigma.sq, and a post-shift scalar (PSS)
## we have to multiply the sigma.sq before shift by the PSS to get the rate after the shift
#######################
actual.after <- NULL
for (i in 1:nrow(trc.sigma)){
  new.after <- NULL
  new.after <- (trc.sigma$before.shift[i])*(trc.sigma$after.shift[i])
  actual.after <- rbind(actual.after, new.after)
}
trc.sigma["after.updated"] <- actual.after

# Compare parameter estimates for Sigma Sq. 
sigma.sq <- (ggplot(trc.sigma)
             + geom_density(aes(x=before.shift),fill="green", alpha=0.5)  # plot param estimates BEFORE SHIFT
             + geom_density(aes(x=after.updated), fill="blue", alpha=0.5, adjust=5)
             + xlim(0,0.02)
             + labs(x="Sigma Squared")) # plot param estimates AFTER  SHIFT

# Compare parameter estimates for Alpha
alpha <- (ggplot(trc.alpha)
          #+ geom_density(aes(x=before.shift),fill="green", alpha=0.2)  # plot param estimates BEFORE SHIFT
          + geom_density(aes(x=after.shift), fill="yellow", alpha=0.5)
          + xlim(0,3)
          + labs(x="Alpha")) # plot param estimates AFTER  SHIFT

# Compare parameter estimates for THETA (change the x-limit if necessary)
size.range <- range(trc.theta$before.shift)
theta<-(ggplot(trc.theta)
        + geom_density(aes(x=before.shift),fill="red", alpha=0.5)  # plot param estimates BEFORE SHIFT
        #+ geom_density(aes(x=after.shift), fill="blue", alpha=0.2)
        + xlim(1.18,1.24) # change this as necessary (limits on the x-axis)
        #+ xlim(size.range[1], size.range[2])
        + labs(x="Theta")) # plot param estimates AFTER  SHIFT

multiplot(sigma.sq, alpha, theta, cols=3)



#######################################################
## Plot Density Distribution plots of Model Shift Times
#######################################################
## Read in your data
model.res <- read.csv("/YOUR_DIRECTORY/Meliphagides.Model.Comparison.BEST.csv")
# you can also use the dating estimates from "input" above

## Plot the timing of shifts as a frequency histogram
(ggplot(model.res, aes(x=timing))
  + geom_histogram(binwidth=1))

## Plot the timing of shifts as a density distribution (with mean estimate)
times <- (ggplot(model.res, aes(x=timing)) 
          + geom_density(fill="orange", alpha=0.5, adjust=2) # the adjust feature will smooth the distribution
          #+ geom_vline(aes(xintercept=mean(timing, na.rm=T)),
          #             color="red", linetype="dashed", size=1)
          + scale_x_reverse(limits=c(20,0))
          + labs(x="Time of Shift")) # use this to reverse then define the limits of the x axis

multiplot(sigma.sq, alpha, theta, times, cols=2)

## Plot the timing of shifts as both a frequency histogram AND a density distribution
(ggplot(model.res, aes(x=timing))
  + geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white")
  + geom_density(alpha=0.2, fill="green")
  + scale_x_reverse(limits=c(20,0))) # use this to reverse then define the limits of the x axis

########################################################################
## Doing these plots for multiple groups (our example is the shift time)
########################################################################
## Read in your data
model.mult.res <- read.csv("/YOUR_DIRECTORY/UPDATED.All.Radiations.Model.Comparison.csv", header=T)
summary(model.mult.res$radiation)
model.mult.res$radiation <- factor(model.mult.res$radiation, levels=c("Pygopodoidea", "Agamidae", "Meliphagoids", "Oz.Marsupials", "Sphenomorphine"))

## Plot the timing of shifts as a density distribution, across groups defined in the "radiation" column
(ggplot(model.mult.res, aes(timing, fill=radiation, colour=radiation))
  + geom_density(alpha=0.5, adjust=1.2) # use alpha to change the opacity
  + scale_x_reverse(limits=c(70,0))
  + scale_fill_manual( values=wes_palette("Zissou", 4, "continuous"))
  + scale_color_manual(values=wes_palette("Zissou", 4, "continuous"))
  + theme_classic() # use this to reverse then define the limits of the x axis
  + geom_vline(xintercept=c(65,64,30,28))) # add a vertical line for the origination of each radiation
lines()

## Plot the timing of shifts as a distribution in a single joyplot using 'ggjoy'
model.mult.res <- subset(model.mult.res, !is.na(model.mult.res$timing)) #make the timing column 'as.numeric' if the plot is jacked up
model.mult.res$timing <- as.numeric(as.character(model.mult.res$timing))
#model.mult.res$radiation <- factor(model.mult.res$radiation, levels=c("Pygopodoidea", "Agamidae", "Meliphagoids", "Oz.Marsupials", "Sphenomorphine"))

(ggplot(model.mult.res, aes(x=timing, y=radiation, fill=radiation, alpha=0.5)) 
  + geom_joy(scale=5)
  + scale_fill_manual(values=wes_palette("Moonrise3", 5, "discrete"))
  + theme_joy()
  + scale_x_reverse())
