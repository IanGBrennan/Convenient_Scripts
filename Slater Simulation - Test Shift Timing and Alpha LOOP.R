library(phytools)
library(mvtnorm)
library(geiger)
library(OUwie)
library(diveRsity)
library(ggplot2)
library(wesanderson)
library(TESS)
library(TreePar)
save.here <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations"
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD

###########################################################
## Step 1. if using empirical trees
###########################################################

## Read in the series of trees:
#vara <- read.nexus("BT.Varanidae.tre") # 28 tips
agam <- read.nexus("BT.Agamids.tre") # 100 tips
mars <- read.nexus("BT.Australian.Marsupials.tre") # 133 tips
bird <- read.nexus("BT.Meliphagides.tre") # 149 tips
pygo <- read.nexus("BT.Pygopodoidea.tre") # 189 tips
skink <- read.nexus("BT.Sphenomorphines.tre") # 240 tips

## Or make your own trees
test1 <- pbtree(n=100, b=1.5, d=1)
test2 <- pbtree(n=200, b=1.5, d=1)
test3 <- pbtree(n=100, b=1, d=0.5)

## Create a multiPhylo object of your empirical trees
emp.trees <- list()
emp.trees <- c(agam, mars, bird, pygo, skink)
class(emp.trees) <- "multiPhylo"

## Create a multiPhylo object of your simulated trees
sim.trees <- list()
sim.trees <- c(test1, test2, test3)
class(sim.trees) <- "multiPhylo"

## Set your parameters relevant to your empirical parameter estimates
# alpha
#alpha = 2 # either a static value
diff.alpha <- seq(from=0.5, to=10, by=0.1) # or set it as a vector of sampled values
diff.alpha <- sample(diff.alpha, size=100, replace=T) # or set it as a vector of sampled values
# pre/post shift sigma
preshift.sigma  = 1
postshift.sigma = 5
# and shift time
sim.shifts <- seq(from=1, to=20, by=1)
sim.shifts <- sample(sim.shifts, size=100, replace=T)

######################################################################################################
## 2. Run through the loop, creating 100 sets of traits (adjusting either time or alpha; Option A + Option C) 
### for each tree. This simulates under a 'TRC' or 'SRC' model, by adjusted Option B (comment out this step)
######################################################################################################
sim.traits.geiger <- list(); sim.traits.ouwie <- list() # make trait lists for all trees in geiger and ouwie data format
# make sure to the outputs depending on if you're simulating different times, or alphas
# this means changing the 'm' and the 'm[[2]]' objects below!
num.sims <- 2 # designate the number of simulated data sets you like to create per tree
input.trees <- output.trees
for (z in 1:length(input.trees)) {
  traits.geiger <- list(); traits.ouwie <- list() # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- input.trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  
  for (i in 1:num.sims) {
    # Option A (comment out top when simulating different shift times, comment out bottom when simulating different alphas)
    m <- split.matrices <- split.vcv(phy, 10) # divide the vcv matrix at a static time
    #m <- split.matrices <- split.vcv(phy, sim.shifts[i]) # or at differing times
    
    # Option B (adjust to change simulating model, comment out both to get the SRC model)
    #m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
    #m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
    
    # Option C (comment out top when simulating different times, comment out bottom when simulating different alphas)
    m[[2]] <- ouMatrix(m[[2]], alpha=diff.alpha[i]) # transform the second half of the vcv matrix according to your alpha
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
  sim.traits.ouwie[[z]] <- traits.ouwie
}

# for transparency, we'll write the data to a file
save(sim.traits.geiger, file="SIMULATED.Trees.and.Data.ALPHA.geigerdata.TRC.RData")
save(sim.traits.ouwie,  file="SIMULATED.Trees.and.Data.ALPHA.ouwiedata.TRC.RData")


####################################################################################
## 3. Now let's create a loop to test the accuracy of time shift parameter estimates
####################################################################################
# remember to do this for each of the 5 tree sizes!
all.timing.res <- NULL
for (j in 1:length(emp.trees)) {
  timing.res <- NULL
  tree <- emp.trees[[j]] # change this to match the tree size you want
  for (i in 1:100) {
    cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
    
    data <- sim.traits.geiger[[j]][[i]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[1]][[i]] # change this to match the proper sized tree
    
    min.shift <- sim.shifts[i]-3; max.shift <- sim.shifts[i]+3
    if(min.shift<=0) {
      min.shift <- 1
    }
    
    ## Fit the Two Rate Constrain (BM/OU) model (this includes saving parameter estimates for each iteration)
    TRC <- NULL
    # create objects to hold the results for each simulation, and just the best fit
    res.trc <- NULL # create a null data frame for the GEIGER Two-Rate Constraint model (BM+OU)
    best.trc <- NULL # create a null data frame for the BEST GEIGER Two-Rate Constraint model (BM+OU)
    
    for (h in min.shift:max.shift) {
      TRCfit<-fitContinuous_paleo(tree, data, model="TRC", shift.time=h)
      res.trc <- rbind(res.trc, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc, h, sim.shifts[i]))))
      #TRC.estimates <- NULL
      #TRC.estimates <- rbind(TRC.estimates, as.data.frame(t(c(TRCfit$Trait1$beta, TRCfit$Trait1$post.shift.scalar))))
      #TRC.estimates <- rbind(TRC.estimates, c(NA, TRCfit$Trait1$alpha))
      #TRC.estimates <- rbind(TRC.estimates, c(time.shift$Trait1$root.state, NA))
      #rownames(TRC.estimates) <- c("sigma.sq", "alpha", "theta")
      #TRC.estimates [,3] <-paste("TRC",h)
      #TRC <- rbind(TRC, TRC.estimates)
    }
    #colnames(TRC) <- c("before.shift", "after.shift", "model")
    
    ## Set the column names of the results data frames
    colnames(res.trc)  <- c("lnL", "AIC", "AICc", "est.time", "sim.time")
    
    ## Set the row names of the results data frames
    #rownames(res.trc) <- c(paste("TRC",min.shift:max.shift))
    
    ## Choose the best (lowest AIC) shift time from each model results data frame
    best.trc <- subset(res.trc, AICc == min(res.trc$AICc))
    timing.res <- rbind(timing.res, best.trc)
  }
  all.timing.res[[j]] <- timing.res
}

write.csv(timing.res, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Agamids.Simulated.Shift.TIMING.csv", quote=F)

head(all.timing.res)
plot(timing.res$est.time, timing.res$sim.time)
lm(timing.res$est.time ~ timing.res$sim.time)
abline(lm(timing.res$est.time ~ timing.res$sim.time))


###########################################################################################
## 4. And now a loop to test the accuracy of alpha parameter estimates
###########################################################################################
# remember to do this for each of the 5 tree sizes!
all.alpha.res <- NULL
#for (i in 1:length(small.traits)) 
for (j in 1:length(emp.trees)){
  alpha.res <- NULL
  tree <- emp.trees[[j]] # change this to match the tree size you want
  #for (i in 1:length(diff.alpha)) {
  for (i in 1:100) {
    cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
    
    data <- sim.traits.geiger[[j]][[i]] # change this to match the proper sized tree
    #data.ouwie <- sim.traits.ouwie[[1]][[i]] # change this to match the proper sized tree
    
    ## Fit the Two Rate Constrain (BM/OU) model (this includes saving parameter estimates for each iteration)
    #TRC <- NULL
    # create objects to hold the results for each simulation, and just the best fit
    res.trc <- NULL # create a null data frame for the GEIGER Two-Rate Constraint model (BM+OU)
    #best.trc <- NULL # create a null data frame for the BEST GEIGER Two-Rate Constraint model (BM+OU)
    
    TRCfit<-fitContinuous_paleo(tree, data, model="TRC", shift.time=10)
    res.trc <- rbind(res.trc, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
    res.trc <- cbind(res.trc, t(c(TRCfit$Trait1$alpha, diff.alpha[i])))
    #for (h in 1:length(diff.alpha)) {
    
    #TRC.estimates <- NULL
    #TRC.estimates <- rbind(TRC.estimates, as.data.frame(t(c(TRCfit$Trait1$beta, TRCfit$Trait1$post.shift.scalar))))
    #TRC.estimates <- rbind(TRC.estimates, c(NA, TRCfit$Trait1$alpha))
    #TRC.estimates <- rbind(TRC.estimates, c(time.shift$Trait1$root.state, NA))
    #rownames(TRC.estimates) <- c("sigma.sq", "alpha", "theta")
    #TRC.estimates [,3] <-paste("TRC",h)
    #TRC <- rbind(TRC, TRC.estimates)
    #}
    #colnames(TRC) <- c("before.shift", "after.shift", "model")
    
    ## Set the column names of the results data frames
    colnames(res.trc)  <- c("lnL", "AIC", "AICc", "est.alpha", "sim.alpha")
    
    ## Set the row names of the results data frames
    #rownames(res.trc) <- c(paste("TRC",min.shift:max.shift))
    
    ## Choose the best (lowest AIC) shift time from each model results data frame
    #best.trc <- subset(res.trc, AICc == min(res.trc$AICc))
    alpha.res <- rbind(alpha.res, res.trc)
  }
  all.alpha.res[[j]] <- alpha.res
}

#write.csv(alpha.res, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Agamids.Simulated.Shift.ALPHA.csv", quote=F)

head(all.alpha.res)
plot(alpha.res$est.alpha, alpha.res$sim.alpha)
lm(alpha.res$est.alpha ~ alpha.res$sim.alpha)
abline(lm(alpha.res$est.alpha ~ alpha.res$sim.alpha))
abline(0,1)



#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
ggplotRegression <- function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.25, color="red") + # change to 0.25 and "red" for time plots
    stat_smooth(method = "lm", col = "black") + # change to "black" for time plots
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
#######################################################################################
# the order: agam, bird, mars, pygo, skink

# write loop to automate this process, otherwise you're dumb
plot.res <- NULL
for (i in 1:length(all.timing.res)) { # change this according to the parameter you simulated
  fit <- lm(est.time ~ sim.time, data=all.timing.res[[i]]) # change this according to the parameter you simulated
  plot.fit <- (ggplotRegression(fit))
  plot.res[[i]] <- plot.fit
}

multiplot(plot.res[[1]], plot.res[[3]],
          plot.res[[2]], plot.res[[4]], 
          plot.res[[5]],
          cols=2)
# plots: (a, b, c, d, e)
# [[a]], [[b]]
# [[c]], [[d]]
# [[e]]




########################################################################################
## 5. Now let's create a loop to comparatively fit a set of models to our simulated data
### this includes the generating model, which we hope provides the best fit!
########################################################################################
# remember to do this for each of the 5 tree sizes!
stree.res <- NULL
input.trees <- output.trees # designate your input trees for this analysis
#for (i in 1:length(small.traits)) 
num.sims <- 2 # designate the number of simulated data sets you like to create per tree
for (i in 1:num.sims) {
  cat("iteration", i, "of", num.sims, "\n") #keep track of what tree/loop# we're on
  
  tree <- input.trees[[1]] # change this to match the tree size you want
  data <- sim.traits.geiger[[1]][[i]] # change this to match the proper sized tree
  data.ouwie <- sim.traits.ouwie[[1]][[i]] # change this to match the proper sized tree
  
  bmfit    <- fitContinuous(tree, data, SE=NA, model="BM")
  ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
  oufit    <- fitContinuous(tree, data, SE=NA, model="OU")
  TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=10)
  SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=10)
  #OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
  #OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
  #OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
  
  
  #####################################################
  ###### Summarize and Compare Model Fitting ##########
  #####################################################
  results.names <- list(bmfit$opt, ebfit$opt, oufit$opt)
  results <- NULL
  for (k in 1:length(results.names)) {
    x <- as.data.frame(results.names[k])
    results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
  }
  
  results.nonstan <- NULL
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
  #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
  #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
  #results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
  
  ## combine both
  results <- rbind.data.frame(results, results.nonstan)
  model <- c("BM", "EB", "OU", "TRC", "SRC") #, "OUM", "OUMA", "OUMV")
  colnames(results) <- c("lnL", "AIC", "AICc")
  results <- cbind(results, model)
  
  ## Use AIC weights to determine best fitting model and model contributions
  weight <- aicw(results$AICc)
  results <- cbind(results, weight$delta)
  results <- cbind(results, weight$w)
  stree.res <- rbind.data.frame(stree.res, results)
}
colnames(stree.res) <- c("lnL", "AIC", "AICc", "model", "delta", "w")
# add in diff.alpha to tell if ones that weren't SRC had low Alpha or something else
diff.alpha
name.of.file <- paste(save.here, "PP.Pygopodoidea.TRC.variable.a.model.selection.csv")
write.csv(stree.res, file=name.of.file, quote=F)





#############
# open the file in excel ("____.Model.Comparison.TOTAL.csv"):
# change first column name to "model", & remove numbers
# and last two columns to "delta" and "w"
#############

#output <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Trimmed.FINAL.Agamidae.Model.Comparison.TOTAL.csv", header=T)
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
#output.total <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Comparison TOTAL/Trimmed.FINAL.All.Radiations.Model.Comparison.TOTAL.csv", header=T)
#output.total <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Comparison TOTAL/Trimmed.FINAL.Meliphagides.Model.Comparison.TOTAL.csv", header=T)
#output.total[,1] <- NULL # drop first column if it's 'X'
outz <- summarySE(stree.res, measurevar="w", groupvars="model")

## Otherwise Just Plot the model results
########################################
#outz$model <- factor(outz$model, 
                     levels=c("delta", "DensDep", "EB", "fit_env", "fit_bgb", 
                              "TS", "OUS", "OUMA", "OUMV", "OUMVA", 
                              "SRC", "TRC")) # re-order the models in the plot
myplot <- (ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + scale_fill_manual( values=wes_palette("Royal2", 8, "continuous"))
  + theme(axis.text.x=element_text(angle=45, hjust=1)))

myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))

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
           + scale_fill_manual( values=wes_palette("Zissou", 8, "continuous")))



########################################################################################
## 6. Finally, if we've been working with simulated trees (those with fossil tips) and 
### data, we need to remove those fossil taxa from our tree and data files and then 
### rerun our model-fitting analyses again (back to Step 5), to compare results.
########################################################################################
drops <- is.extinct(output.trees[[1]], tol=0.0001) # find out which tips are extinct
extant.sim <- drop.extinct(output.trees[[1]], tol=0.00001)
target.trees <- output.trees
extant.trees <- NULL
for (y in 1:length(target.trees)) {
  extant.trees[[y]] <- drop.extinct(output.trees[[y]], tol=0.00001)
  class(extant.trees) <- "multiPhylo"
}
extant.data <- NULL
for (j in 1:length(sim.traits.geiger)) {
  e.data <- NULL
  for (k in 1: length(sim.traits.geiger[[j]])) {
    drops <- is.extinct(target.trees[[j]], tol=0.0001)
    e.data[[k]] <- subset(sim.traits.geiger[[j]][[k]], !rownames(sim.traits.geiger[[j]][[k]]) %in% drops)
  }
  extant.data[[j]] <- e.data
}


### Estimate birth/death rates for the extant and fossil trees
bd.ms(extant.trees[[1]])
bd.km(extant.trees[[1]])


### Estimate speciation/extinction rates via TESS/TreePar
# load in a test data set
tree<-read.tree("Pygopodoidea.tre")
# convert the phylogeny into the branching times
times<-as.numeric(branching.times(extant.trees[[1]]))
ext.times <- as.numeric(branching.times(output.trees[[1]]))

########################################################################
#### Create a likelihood function for constant rate birth-death model
########################################################################
#### Create the prior distributions as functions of the log-transformed probability
prior_delta <- function(times) { dexp(times,rate=10.0,log=TRUE) } 
prior_tau <- function(times) { dexp(times,rate=10.0,log=TRUE) } 
priorsConstBD <- c("diversification"=prior_delta,"turnover"=prior_tau)

#specify the new constant BD likelihood function
#we can make a pure-birth model by setting the extinction rate to 0 (easy peasy)
likelihoodConstBD<-function(params) {
  speciation<-params[1]+params[2]
  extinction<-params[2]
  lnl<-tess.likelihood(ext.times,
                       lambda=speciation,
                       mu=extinction,
                       samplingProbability=1.0,
                       log=TRUE)
  return(lnl)
}
#### Run the MCMC on the likelihood function of constant birthdeath
samplesConstBD <- tess.mcmc(likelihoodFunction = likelihoodConstBD,
                            priors = priorsConstBD,
                            parameters = runif(2,0,1),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 100000,
                            burnin = 10000,
                            thinning = 10,
                            adaptive = TRUE,
                            verbose = TRUE)

summary(samplesConstBD) #summarize the rate estimates
plot(samplesConstBD) #plot the sampling of our ConstantBD model
colMeans(samplesConstBD) #summarize the column means
head(summary(samplesConstBD))



x <- branching.times(output.trees[[1]])
LikConstant(x=x, lambda=0.9, mu=0.8, sampling=1)
start = log(c(0.1, 0.01),1)
lower = log(c(1e-7, 1e-7), 1)
upper = log(c(10, 10), 1)
o <- optim(LikConstant, p = start, lower = lower, upper = upper, 
           method = "L")
bd.shifts.optim(x,sampling=c(1),survival=1)[[2]]
LikConstant()












  
  #####################################################
  ###### Summarize and Compare Model Fitting ##########
  #####################################################
  results.names <- list(bmfit$opt, ebfit$opt, oufit$opt)
  results <- NULL
  for (k in 1:length(results.names)) {
    x <- as.data.frame(results.names[k])
    results <- rbind(results, as.data.frame(t(c(x$lnL, x$aic, x$aicc))))
  }
  
  results.nonstan <- NULL
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(TRCfit$Trait1$lnl, TRCfit$Trait1$aic, TRCfit$Trait1$aicc))))
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(SRCfit$Trait1$lnl, SRCfit$Trait1$aic, SRCfit$Trait1$aicc))))
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMfit$loglik, OUMfit$AIC, OUMfit$AICc))))
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMAfit$loglik, OUMAfit$AIC, OUMAfit$AICc))))
  results.nonstan <- rbind(results.nonstan, as.data.frame(t(c(OUMVfit$loglik, OUMVfit$AIC, OUMVfit$AICc))))
  
  ## combine both
  results <- rbind.data.frame(results, results.nonstan)
  model <- c("BM", "EB", "OU", "TRC", "SRC", "OUM", "OUMA", "OUMV")
  colnames(results) <- c("lnL", "AIC", "AICc")
  results <- cbind(results, model)
  
  ## Use AIC weights to determine best fitting model and model contributions
  weight <- aicw(results$AICc)
  results <- cbind(results, weight$delta)
  results <- cbind(results, weight$w)
  tree.res <- rbind.data.frame(tree.res, results)
}
colnames(tree.res) <- c("lnL", "AIC", "AICc", "model", "delta", "w")
name.of.file <- paste(save.here, "PP.Pygopodoidea.TRC.a2.model.selection.csv")
write.csv(tree.res, file=name.of.file, quote=F)




















## 5. ## Summarize and plot
##############################################################################
# If you're coming back, read in your file ("PP.____tips.model.selection.csv"):
# output <- read.csv("PP.25tips.model.selection.csv", header=T)
#############
output <- stree.res # otherwise use the object you created above
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
#output.total <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/FINAL.Agamidae.Model.Comparison.TOTAL.csv", header=T)
# Use the recent object
outz <- summarySE(output, measurevar="w", groupvars="model")
# Or load and summarize
small.output <- read.csv("PP.25tips.SRC.a0.2.model.selection.csv", header=T)
small.out <- summarySE(small.output, measurevar="w", groupvars="model")
smedium.output <- read.csv("PP.100tips.SRC.a0.2.model.selection.csv", header=T)
smedium.out <- summarySE(smedium.output, measurevar="w", groupvars="model")

## If you want to plot a composite bar graph, subset each radiation, summarize model weights
outz.mars <- summarySE(subset(output.total, group=="marsupials"), measurevar="w", groupvars="model")
outz.agam <- summarySE(subset(output.total, group=="agamidae"), measurevar="w", groupvars="model")
outz.pygo <- summarySE(subset(output.total, group=="pygopodoidea"), measurevar="w", groupvars="model")
outz.skink <- summarySE(subset(output.total, group=="sphenomorphine"), measurevar="w", groupvars="model")
## Then combine them back together
outz.mars["group"] <- "marsupials"
outz.agam["group"] <- "agamidae"
outz.pygo["group"] <- "pygopodoidea"
outz.skink["group"] <- "sphenomorphine"
group.outz <- rbind.data.frame(outz.mars, outz.agam, outz.pygo, outz.skink)

## If you want to summarize each model into a table (you probably don't)
###################################################
TRCs <- subset(outz, model == "TRC")
SRCs <- subset(outz, model == "SRC")
shift.mod <- rbind.data.frame(TRCs, SRCs)
##################################################

## Subset the SRC model from all
SRC.small <- subset(small.out, model == "SRC")
SRC.small["ntips"] <- "25"
SRC.smedium <- subset(smedium.out, model == "SRC")
SRC.smedium["ntips"] <- "50"
all.SRC <- rbind(SRC.small, SRC.smedium)

## Otherwise Just Plot the model results
########################################
(ggplot(outz, aes(x=model, y=w, fill=model))
 + geom_bar(stat="identity")
 + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
 + theme(axis.text.x=element_text(angle=45, hjust=1))
 + theme_classic())
#+ scale_fill_manual(values=wes_palette("Darjeeling")))

(ggplot(shift.mod, aes(x=model, y=w, fill=model))
  + geom_boxplot(stat="identity"))

(ggplot(shift.mod, aes(x=model, y=w))
  + geom_point())

## Plot the composite bar graphs
group.outz$model <- factor(group.outz$model, 
                           levels=c("BM", "delta", "kappa", "lambda", 
                                    "gamma", "EB", "DensDep", "OU",
                                    "BMS", "TS", "OUS", "OUMA", "OUMV",
                                    "OUMVA", "SRC", "TRC")) # this re-orders the models in the legend
(ggplot(group.outz)
  + geom_bar(aes(y=w, x=group, fill=group.outz$model), stat="identity")
  + theme(axis.text.x=element_text(angle=25, hjust=1))
  + scale_fill_manual( values=wes_palette("Cavalcanti", 15, "continuous")))


########## next is fix the loop to write the ouwie format! ################