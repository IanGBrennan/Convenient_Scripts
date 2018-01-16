save.here <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations"

###########################################################
## Step 1. if using empirical trees
###########################################################

## Read in the series of trees:
vara <- read.nexus("BT.Varanidae.tre") # 28 tips
bird <- read.nexus("BT.Honeyeaters.tre") # 58 tips
agam <- read.nexus("BT.Agamids.tre") # 100 tips
mars <- read.nexus("BT.Australian.Marsupials.tre") # 133 tips
pygo <- read.nexus("BT.Pygopodoidea.tre") # 189 tips
skink <- read.nexus("BT.Sphenomorphines.tre") # 240 tips

## Create a multiPhylo object of your trees
emp.trees <- list()
emp.trees <- c(vara, bird, agam, mars, pygo, skink)
class(emp.trees) <- "multiPhylo"

## Set your parameters
alpha = 2
preshift.sigma  = 1
postshift.sigma = 5

## Run through the loop, creating 100 sets of traits for each tree
traits.geiger <- list(); traits.ouwie <- list() # make trait lists for geiger and ouwie data format
for (i in 1:length(emp.trees)) {
  phy <- emp.trees[[i]] # designating the target tree
  m <- split.matrices <- split.vcv(phy, 10) # divide the vcv matrix
  
  # like the 'release radiate' model
  m[[2]] <- ouMatrix(m[[2]], alpha=alpha) # transform the second half of the vcv matrix according to your alpha
  m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
  m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
  m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
  
  # OR, do it like the 'ecological release' model
  # m <- lapply(m, function(x) x*0.1)
  # m.rev.rel <- m[[1]] + m[[2]]
  
  traitz <- list(); traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  traitz.geiger <- NULL; traitz.ouwie <- NULL
  for (t in 1:100){
    # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
    traitz[[t]] <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), 
                                    sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))
    t.ouwie <- NULL
    t.ouwie <- as.data.frame(names(traitz[[t]]))
    t.ouwie[,2] <- as.data.frame(t(traitz[[t]]))
    traitz.ouwie[[t]] <- t.ouwie
    
    t.geiger <- t.ouwie; names <- t.ouwie[,1]
    rownames(t.geiger) <- names; t.geiger[,1] <- NULL
    traitz.geiger[[t]] <- t.geiger
    
  }
  traits.geiger[[i]] <- traitz.geiger
  traits.ouwie[[i]]  <- traitz.ouwie
}

# for transparency, we'll write the data to a file
save(traits.geiger, file="Slater.SIMULATED.geigerdata.TRC.RData")
save(traits.ouwie,  file="Slater.SIMULATED.ouwiedata.TRC.RData")


## 4. Now let's create a loop to comparatively fit the models to our simulated data
################################################################################
# remember to do this for each of the 5 tree sizes!
tree.res <- NULL
#for (i in 1:length(small.traits)) 
for (i in 1:100) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  
  tree <- emp.trees[[5]] # change this to match the tree size you want
  data <- traits.geiger[[5]][[i]] # change this to match the proper sized tree
  data.ouwie <- traits.ouwie[[5]][[i]] # change this to match the proper sized tree
  
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
  tree.res <- rbind.data.frame(tree.res, results)
}
colnames(tree.res) <- c("lnL", "AIC", "AICc", "model", "delta", "w")
name.of.file <- paste(save.here, "PP.Pygopodoidea.TRC.a2.model.selection.csv")
write.csv(tree.res, file=name.of.file, quote=F)

## 4. Now let's create a loop to comparatively fit the models to our simulated data
################################################################################
# remember to do this for each of the 5 tree sizes!
stree.res <- NULL
#for (i in 1:length(small.traits)) 
for (i in 1:100) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  
  tree <- emp.trees[[6]] # change this to match the tree size you want
  data <- traits.geiger[[6]][[i]] # change this to match the proper sized tree
  data.ouwie <- traits.ouwie[[6]][[i]] # change this to match the proper sized tree
  
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
name.of.file <- paste(save.here, "PP.Sphenomorphine.TRC.a2.model.selection.csv")
write.csv(stree.res, file=name.of.file, quote=F)



















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