library(mvtnorm)
library(geiger)
library(OUwie)
library(diveRsity)
library(ggplot2)
library(wesanderson)
save.here <- "/YOUR_DIRECTORY/Trait Simulations"
source("https://github.com/IanGBrennan/MioceneAustralia/blob/master/Scripts/Phenotypic_Evolution/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD

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

## Create a multiPhylo object of your trees
emp.trees <- list()
emp.trees <- c(agam, mars, bird, pygo, skink)
class(emp.trees) <- "multiPhylo"

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
for (z in 1:length(emp.trees)) {
  traits.geiger <- list(); traits.ouwie <- list() # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- emp.trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  
  for (i in 1:100) {
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
save(sim.traits.geiger, file="Slater.SIMULATED.TIMING.geigerdata.TRC.RData")
save(sim.traits.ouwie,  file="Slater.SIMULATED.TIMING.ouwiedata.TRC.RData")


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

write.csv(timing.res, file="/YOUR_DIRECTORY/Trait Simulations/Agamids.Simulated.Shift.TIMING.csv", quote=F)

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
#for (i in 1:length(small.traits)) 
for (i in 1:100) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  
  tree <- emp.trees[[4]] # change this to match the tree size you want
  data <- sim.traits.geiger[[4]][[i]] # change this to match the proper sized tree
  data.ouwie <- sim.traits.ouwie[[4]][[i]] # change this to match the proper sized tree
  
  bmfit    <- fitContinuous(tree, data, SE=NA, model="BM")
  ebfit    <- fitContinuous(tree, data, SE=NA, model="EB")
  oufit    <- fitContinuous(tree, data, SE=NA, model="OU")
  TRCfit   <- fitContinuous_paleo(tree, data, model="TRC", shift.time=10)
  SRCfit   <- fitContinuous_paleo(tree, data, model="SRC", shift.time=10)
  OUMfit  <- OUwie.slice(tree, data.ouwie, model=c("OUM"),  root.station=T, timeslices=c(10))
  OUMAfit <- OUwie.slice(tree, data.ouwie, model=c("OUMA"), root.station=T, timeslices=c(10))
  OUMVfit <- OUwie.slice(tree, data.ouwie, model=c("OUMV"), root.station=T, timeslices=c(10))
  
  
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

#output <- read.csv("/YOUR_DIRECTORY/Trimmed.FINAL.Agamidae.Model.Comparison.TOTAL.csv", header=T)
#table(output, output$model)
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