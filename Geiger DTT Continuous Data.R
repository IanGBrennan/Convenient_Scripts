library(phytools)
library(geiger)
library(mvtnorm)

#######################################################
# Read in all the data you'll use
######################################################
tree<-read.nexus("BT.Australian.Marsupials.tre") #read in desired tree
#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data<-read.csv("BT.Australian.Marsupials.MlogBL.csv", row.names = 1, header=F) #read in data file

#logBL<-setNames(data[,14], rownames(data)) #choose your data column (here: logBL) and apply rownames
#name.check(tree, data["M.logBL"]) #check to make sure the tips match the data labels
name.check(tree, data)

##################################################################
# Quickly view the data by doing an Ancestral State Reconstruction
##################################################################
#quickly view the data by doing an Ancestral State Reconstruction of the focal data
fit<-fastAnc(tree,data,vars=TRUE,CI=TRUE) #run the quick ASR (Anc)
obj<-contMap(tree,logWT,plot=FALSE) #set it to an object
#plot(obj,type="fan",legend=0.7*max(nodeHeights(tree)),
#     fsize=c(0.7,0.9)) #plot it as a fan if you'd like
plot(obj, fsize=0.5) #otherwise just plot it here, and adjust values as you see fit
#################### change the colors so they don't look like doo-doo
## what is the length of the current color ramp?
n<-length(obj$cols)
## change to grey scale
obj$cols[1:n]<-grey(0:(n-1)/(n-1))
plot(obj)## change to blue -> red
obj$cols[1:n]<-colorRampPalette(c("red", "pink", "white", "lightblue", "navy"), space="Lab")(n)
plot(setMap(obj,invert=TRUE), fsize = 0.4, lwd=3, legend = TRUE, type = "phylogram", outline=TRUE, direction="rightwards")


######################################################
# run DTT on the data that you checked above
######################################################
dtt.obj <- dtt(tree, data, nsim=1000, index=c("avg.sq"), plot=TRUE, calculateMDIp=TRUE)
dtt.obj$MDI; dtt.obj$MDIpVal #have a look at the most important values

bmfit$opt
######################################################
# run Continuous model fitting on the data
######################################################
bmfit<-fitContinuous(tree, data, SE=NA, model="BM", control=(niter=1000))
ebfit<-fitContinuous(tree, data, SE=NA, model="EB", control=(niter=1000))
deltafit<-fitContinuous(tree, data, SE=NA, model="delta", control=(niter=1000))
kappafit<-fitContinuous(tree, data, SE=NA, model="kappa", control=(niter=1000))
oufit<-fitContinuous(tree, data, SE=NA, model="OU", control=(niter=1000))
lambdafit<-fitContinuous(tree, data, SE=NA, model="lambda", control=(niter=1000))

######################################################
######## Slater Paleo Code from 2013 Mammal Paper ####
######################################################
#the fxn 'fitContinuous_paleo' allows timed shifts and release of trait evolution
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Slater.2013.Fit Release.Radiate Model.R"); ## now source the function from your WD

# the additional model options are:
# "timeshift", "release", "releaseradiate"
time.shift<-fitContinuous_paleo(tree, data, model="timeshift", shift.time=6)
Radiate.Mio.constrain<-fitContinuous_paleo(tree, data, model="radiate.constrain", shift.time=4)
Release.Mio.constrain<-fitContinuous_paleo(tree, data, model="release.constrain", shift.time=5)

fitContinuous_paleo(tree, data, model="release", shift.time=29)
fitContinuous_paleo(tree, data, model="BM")

#####################################################
###### Summarize and Compare Model Fitting ##########
#####################################################

results<-as.data.frame(t(c(bmfit$opt$lnL, bmfit$opt$aic, bmfit$opt$aicc)))
results<-rbind(results, c(as.data.frame(t(c(ebfit$opt$lnL, ebfit$opt$aic, ebfit$opt$aicc)))))                  
results<-rbind(results, c(as.data.frame(t(c(oufit$opt$lnL, oufit$opt$aic, oufit$opt$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(deltafit$opt$lnL, deltafit$opt$aic, deltafit$opt$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(kappafit$opt$lnL, kappafit$opt$aic, kappafit$opt$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(lambdafit$opt$lnL, lambdafit$opt$aic, lambdafit$opt$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(time.shift$Trait1$lnl, time.shift$Trait1$aic, time.shift$Trait1$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(Radiate.Mio.constrain$Trait1$lnl, Radiate.Mio.constrain$Trait1$aic, Radiate.Mio.constrain$Trait1$aicc)))))
results<-rbind(results, c(as.data.frame(t(c(Release.Mio.constrain$Trait1$lnl, Release.Mio.constrain$Trait1$aic, Release.Mio.constrain$Trait1$aicc)))))

names(results)<-c("lnL", "AIC", "AICc")
weights <- aicw(results$AICc)
results <- cbind(results, weights)
rownames(results)<-c("BM", "EB", "OU", "Delta", 
                     "Kappa", "Lambda", "EOb.TimeShift", 
                     "Radiate.Miocene.Constrain", "Release.Miocene.Constrain")
results
write.csv(results, file="xxx.csv")


fitContinuousMCMC(tree, logBL, model = "BM", Ngens = 100000, sampleFreq =1000,
                  printFreq = 1000, node.priors = NULL, root.prior = NULL,
                  outputName ="BM_M.logBL")
fitContinuousMCMC(tree, logBL, model = "ACDC.lin", Ngens = 10000, sampleFreq =100,
                  printFreq = 1000, node.priors = NULL, root.prior = NULL,
                  outputName ="ACDC.lin_M.logBL")
bm.res <- read.table("BM_M.logBL_model_params.txt", header= TRUE)
head(bm.res)


plot(bm.res$generation, bm.res$logLk, type = "l")

lines(trend.res$generation, trend.res$logLk, col = "red")
legend("bottomleft", c("bm", "trend"), lwd = 3, col = c("black", "red"))


###################################################
######## THE STANLEY ADJUSTED CODE (different) ####
###################################################

##### Creates Ed's MDIbyNode function #####
MDIplot<-function(phy,data,SimNo=100,TITLE=NA,sigfigs=3,SaveFile="Y"){
  length(phy$tip.label)+1->N
  length(phy$tip.label)+phy$Nnode->M
  MDIValues<-data.frame()
  
  for (i in N:M)
  {
    MDI<-0
    extract.clade(phy,i)->phy1
    if(length(phy1$tip.label)>2){
      dtt(phy1,data,nsim=100)->DTT
      DTT$MDI->MDI
      signif(MDI,digits = sigfigs)->MDI
    }
    append(MDIValues,MDI)->MDIValues
  }
  makeNodeLabel(phy,"n")->tr
  tr$node.label<-as.list(MDIValues)
  plot(tr,show.node.label = TRUE)
  write.tree(tr,"MDI.tre")
  return(tr)}


MDItree<-MDIplot(tree,data["M.logBL"]) #plot the resulting tree with node values as MDI
write.tree(MDItree, file="Marsupials.MDItree.tre") #write the tree to view in FigTree

#alternatively, you can view this info with the plotBranchbyTrait function
plotBranchbyTrait(MDItree, MDItree$node.label, palette="rainbow")


###################################################
######## THE SLATER ADJUSTED CODE (better) ########
#*use with caution as it does something funny with CIs
###################################################
## this script will repeat the analysis of disparity through time conducted in Slater et al. 2010 and reproduce figure 4
pdf("~/Google.Drive/R.Analyses/Marsupials.DTT.pdf") #setup your output file
svl <- read.csv("Marsupials.data.csv", row.names=1, header=TRUE)
data<-setNames(svl[,15], rownames(svl)) #choose your data file (here: svl) and column ([,5]), apply rownames
source("/Users/Ian/Google.Drive/R.Analyses/GEIGER/dtt.full with confidence limits and p value.R"); ## now source the function from your WD
# and run the analysis. i recommend using ~10,000 sims to get a stable P value
##if you've used the regular dtt.full, you'll notice that this function doesn't immediately plot up the empirical dtt curve. don't worry - it's just that the way I get the nice shaded 95% region requires a blank plot at first with the empirical and median curves overlain on the shaded region.
dttOut<-dttFullCIs(tree, data, nsims=1000);
title(main="Weight (log mass) - disparity through time, avg.sq.") #apply title to your plot
dev.off()
dttOut$MDI; dttOut$Pvalue


