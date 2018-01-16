library(ape)
library(phytools)
library(OUwie)

## I'm leaving this here to test out against a known phylogeny/result
## Test out the models against the Whales (Slater, 2010)
#######################################
data(whales)
whaletree<- whales$phy ## the whale tree
data<- setNames(whales$dat[,2], rownames(whales$dat)); ## the length data in log(meters)
data.OUwie <- read.csv("Whales.Data.csv", header=F) #read in data file in OUwie format
whaletree <- drop.tip(whaletree, c(	"Balaenoptera_bonaerensis_X75581",
                                    "Balaenoptera_brydei",
                                    "Cephalorhynchus_eutropia_AF084072",
                                    "Delphinus_capensis_AF084087_",
                                    "Delphinus_tropicalis_AF084088",
                                    "Eubalaena_japonica",
                                    "Inia_geoffrensis_humboldtiana_AF521110",
                                    "Mesoplodon_traversii_AY579556",
                                    "Sotalia_guianensis_DQ086827"))
#######################################

save.here <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Simulations/"

##########################################################################
## We're going to:
# 1. simulate phylogenies (I use 100) of varied size (25-400 tips), and keep them in a multiPhylo object
# 2. simmap each tree with two (you can do more) temporal regimes (at shift time 't')
#    using 'make.era.map'
# 3. simulate a trait across each phylogeny, according to a set of OU parameters for each regime
#    and format them for both geiger and OUwie
# 4. estimate the fit (ML) of the models to all of our trees/data in a big ol' loop
# 5. summarize and plot the outputs
##########################################################################

## 1. ## Start by setting up a series of trees, and simulating them with 'pbtree'
##########################################################################
size = c(25, 50, 100, 200, 400) # set your tree sizes
#treesize = c("small", "smedium", "medium", "marge", "large") # if you wanna name em, it's handy!
trees <- list() # set up your list for the trees
for (i in 1:length(size)) {
  trees[[i]] <- pbtree(b=0.1, d=0.05, scale=50, n=size[i], extant.only=T, nsim=100)
  #assign(paste(treesize[i]), tree) # use this if you want to name them
  class(trees[[i]]) <- "multiPhylo" # set the trees as a multiPhylo object!
} # loop through, and add the trees to your list
# note: trees[[1]]=small, trees[[2]]=smedium, trees[[3]]=medium....

## 2. ## Determine where you want to make your shift point, then SIMMAP
##########################################################################
max.height <- max(nodeHeights(trees[[5]][[1]]))
timeslices <- max.height - 10
timeslices <- c(0,timeslices)

era.trees <- list() # make an empty list for your simmap era trees
for (z in 1:length(size)) {
  era.trees[[z]] <- trees[[z]]
  for (i in 1:100) {
    era.trees[[z]][[i]] <- make.era.map(trees[[z]][[i]], timeslices)
  }
} # loop through, simmap (era.map) each tree in the multiPhylo

## If you want to plot your simmap trees for a figure
colorz <- wes_palette("Darjeeling")
names(colorz) <- c(1,2)
par(mfrow=c(1,5))
for (i in 1:length(era.trees)) {
  plotSimmap(era.trees[[i]][[i]], colorz, pts=F, ftype="off", lwd=1)
}

######################################################
## 3. Now we'll simulate the traits onto our phylogenies
######################################################

## Identify your parameter values:
p.alpha <- c(2, 1e-10)
p.sigma.sq <- c(0.1, 0.01)
p.theta0 <- 2
p.theta <- c(2,2)

## Start with the smallest tree and work up to the biggest
small.traits <- list(); small.traits.ouwie <- list()
smedium.traits <- list(); smedium.traits.ouwie <- list()
medium.traits <- list(); medium.traits.ouwie <- list()
marge.traits <- list(); marge.traits.ouwie <- list()
large.traits <- list(); large.traits.ouwie <- list()

for (i in 1:100) {
  # simulate the small traits first
  small.traits[[i]] <- OUwie.sim(era.trees[[1]][[i]], simmap=T, alpha=p.alpha, sigma.sq=p.sigma.sq, theta0=p.theta0, theta=p.theta)
  # set up the small traits for 'OUwie'
  small.traits.ouwie[[i]] <- small.traits[[i]]
  # then get the small data in 'geiger' order
  small.traits[[i]] <- as.data.frame(small.traits[[i]], row.names=era.trees[[1]][[i]]$tip.label)
  small.traits[[i]][,1] <- NULL
  
  # simulate the smedium traits first
  smedium.traits[[i]] <- OUwie.sim(era.trees[[2]][[i]], simmap=T, alpha=p.alpha, sigma.sq=p.sigma.sq, theta0=p.theta0, theta=p.theta)
  # set up the smedium traits for 'OUwie'
  smedium.traits.ouwie[[i]] <- smedium.traits[[i]]
  # then get the smedium data in 'geiger' order
  smedium.traits[[i]] <- as.data.frame(smedium.traits[[i]], row.names=era.trees[[2]][[i]]$tip.label)
  smedium.traits[[i]][,1] <- NULL
  
  # simulate the medium traits first
  medium.traits[[i]] <- OUwie.sim(era.trees[[3]][[i]], simmap=T, alpha=p.alpha, sigma.sq=p.sigma.sq, theta0=p.theta0, theta=p.theta)
  # set up the medium traits for 'OUwie'
  medium.traits.ouwie[[i]] <- smedium.traits[[i]]
  # then get the medium data in 'geiger' order
  medium.traits[[i]] <- as.data.frame(medium.traits[[i]], row.names=era.trees[[3]][[i]]$tip.label)
  medium.traits[[i]][,1] <- NULL
  
  # simulate the marge traits first
  marge.traits[[i]] <- OUwie.sim(era.trees[[4]][[i]], simmap=T, alpha=p.alpha, sigma.sq=p.sigma.sq, theta0=p.theta0, theta=p.theta)
  # set up the marge traits for 'OUwie'
  marge.traits.ouwie[[i]] <- marge.traits[[i]]
  # then get the marge data in 'geiger' order
  marge.traits[[i]] <- as.data.frame(marge.traits[[i]], row.names=era.trees[[4]][[i]]$tip.label)
  marge.traits[[i]][,1] <- NULL
  
  # simulate the large traits first
  large.traits[[i]] <- OUwie.sim(era.trees[[5]][[i]], simmap=T, alpha=p.alpha, sigma.sq=p.sigma.sq, theta0=p.theta0, theta=p.theta)
  # set up the large traits for 'OUwie'
  large.traits.ouwie[[i]] <- large.traits[[i]]
  # then get the large data in 'geiger' order
  large.traits[[i]] <- as.data.frame(large.traits[[i]], row.names=era.trees[[5]][[i]]$tip.label)
  large.traits[[i]][,1] <- NULL
}


## 4. Now let's create a loop to comparatively fit the models to our simulated data
################################################################################
# remember to do this for each of the 5 tree sizes!
tree.res <- NULL
#for (i in 1:length(small.traits)) 
for (i in 1:100) {
  cat("iteration", i, "of 100", "\n") #keep track of what tree/loop# we're on
  
  #tree <- era.trees[[1]][[i]] # change this to match the tree size you want
  tree <- era.tree # change this to match the tree size you want
  data <- small.traits[[i]] # change this to match the proper sized tree
  data.ouwie <- small.traits.ouwie[[i]] # change this to match the proper sized tree
  
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
name.of.file <- paste(save.here, "PP.25tips.SRC.a0.2.model.selection.csv")
write.csv(tree.res, file=name.of.file, quote=F)


## 5. ## Summarize and plot
##############################################################################
# If you're coming back, read in your file ("PP.____tips.model.selection.csv"):
# output <- read.csv("PP.25tips.model.selection.csv", header=T)
#############
output <- tree.res # otherwise use the object you created above
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
 + theme(axis.text.x=element_text(angle=45, hjust=1)))
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








phy <- pbtree(b=0.1, d=0.05, t=50)
phy <- pbtree(b=0.1, scale=50, n=50, extant.only=T)
ephy<-make.era.map(phy, limits = c(0, 25)) # make an era map with a shift half way through the clade’s duration
plotSimmap(ephy)

y <- OUwie.sim(ephy, simmap.tree=TRUE, root.age=50, alpha=c(100,1e-100), sigma.sq=c(0.5,0.5), theta0=0, theta = c(0,1)) 
x <- OUwie.sim(ephy, simmap.tree=TRUE, root.age=50, alpha=c(1, 10^10), sigma.sq=c(0.5,0.5), theta0=0, theta = c(0,1))

plot(diag(vcv(phy)), y[,2])
plot(diag(vcv(phy)), x[,2])
dev.off()

test.tree <- read.nexus("BT.Pygopodoidea.tre")
test.data <- read.csv("BT.Pygopodoidea.logSVL.csv", header=F)
plot(diag(vcv(test.tree)), test.data[,2])
