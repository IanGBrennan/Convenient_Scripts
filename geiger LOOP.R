library(phytools)
library(geiger)
library(mvtnorm)
library(ggplot2)

#######################################################
# Read in all the data you'll use
######################################################
#tree<-read.nexus("BT.Pygopodoidea.tre") #read in desired tree
trees <- read.nexus("PB.Australian.Marsupials.100.trees")
#trees<-lapply(trees, drop.tip, tip=drip) #drop tips if necessary


#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data<-read.csv("BT.Australian.Marsupials.MlogBL.csv", row.names = 1, header=F) #read in data file

#logBL<-setNames(data[,14], rownames(data)) #choose your data column (here: logBL) and apply rownames
#name.check(tree, data["M.logBL"]) #check to make sure the tips match the data labels
name.check(trees[[1]], data)

total.results <- NULL
tree.model.comparison <- NULL
for (i in 1:100) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  
  ######################################################
  # run Continuous model fitting on the data
  ######################################################
  bmfit<-fitContinuous(trees[[i]], data, SE=NA, model="BM", control=(niter=1000))
  ebfit<-fitContinuous(trees[[i]], data, SE=NA, model="EB", control=(niter=1000))
  deltafit<-fitContinuous(trees[[i]], data, SE=NA, model="delta", control=(niter=1000))
  kappafit<-fitContinuous(trees[[i]], data, SE=NA, model="kappa", control=(niter=1000))
  oufit<-fitContinuous(trees[[i]], data, SE=NA, model="OU", control=(niter=1000))
  lambdafit<-fitContinuous(trees[[i]], data, SE=NA, model="lambda", control=(niter=1000))
  ######################################################
  ######## Slater Paleo Code from 2013 Mammal Paper ####
  ######################################################
  #the fxn 'fitContinuous_paleo' allows timed shifts and release of trait evolution
  source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Slater.2013.Fit Release.Radiate Model.R"); ## now source the function from your WD
  
  res.ts <- NULL
  res.trc <- NULL
  res.src <- NULL
  best.ts <- NULL
  best.trc <- NULL
  best.src <- NULL

  for (h in 3:11) {
    time.shift<-fitContinuous_paleo(trees[[i]], data, model="timeshift", shift.time=h)
    res.ts <- rbind(res.ts, as.data.frame(t(c(time.shift$Trait1$lnl, time.shift$Trait1$aic, time.shift$Trait1$aicc))))
  }
  for (h in 3:11) {
    TRCfit<-fitContinuous_paleo(trees[[i]], data, model="radiate.constrain", shift.time=h)
    res.trc <- rbind(res.trc, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
  }
  for (h in 3:11) {
    SRCfit<-fitContinuous_paleo(trees[[i]], data, model="release.constrain", shift.time=h)
    res.src <- rbind(res.src, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
  }
  
  colnames(res.ts) <- c("lnL", "AIC", "AICc")
  colnames(res.trc) <- c("lnL", "AIC", "AICc")
  colnames(res.src) <- c("lnL", "AIC", "AICc")
  
  rownames(res.ts) <- c(paste("ts",3:11))
  rownames(res.trc) <- c(paste("trc",3:11))
  rownames(res.src) <- c(paste("src",3:11))
  
  best.ts <- subset(res.ts, AICc == min(res.ts$AICc))
  best.trc <- subset(res.trc, AICc == min(res.trc$AICc))
  best.src <- subset(res.src, AICc == min(res.src$AICc))
  
  
  #####################################################
  ###### Summarize and Compare Model Fitting ##########
  #####################################################
  
  results <- NULL
  results.names <- list(bmfit$opt, ebfit$opt, oufit$opt,
                        kappafit$opt, deltafit$opt, lambdafit$opt)
  non.standard <- list(time.shift$Trait1, TRCfit$Trait1, SRCfit$Trait1)
  for (k in 1:length(results.names)) {
    x <- as.data.frame(results.names[k])
    results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
  }
  rownames(results) <- c("BM", "EB", "OU", "kappa", "delta", "lambda")
  colnames(results) <- c("lnL", "AIC", "AICc")
  
  results <- rbind(results, best.ts)
  results <- rbind(results, best.trc)
  results <- rbind(results, best.src)
  
  weight <- aicw(results$AICc)
  results <- cbind(results, weight$delta)
  results <- cbind(results, weight$w)

  total.results <- rbind(total.results, results)
  
  score <- subset(results, AICc == min(results$AICc))
  tree.model.comparison <- rbind(tree.model.comparison, score)
}

write.csv(tree.model.comparison, 
          file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Australian.Marsupials.Model.Comparison.csv",
          quote=F)
write.csv(total.results,
          file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Australian.Marsupials.Model.Comparison.TOTAL.csv",
          quote=F)

#############
# open the files in excel:
# change first column name to "model", & remove numbers
# and last two columns to "delta" and "w"
#############

#output <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Agamidae.Model.Comparison.csv", header=T)
#table(output, output$model)

# Summarize the table with this Standard Error Function
#############################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###################################################
output.total <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Australian.Marsupials.Model.Comparison.TOTAL.csv", header=T)
outz <- summarySE(output.total, measurevar="w", groupvars="model")

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
(ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2))
