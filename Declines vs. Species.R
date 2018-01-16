library(laser)
library(geiger)
library(ape)
library(BAMMtools)
library(phytools)
library(TESS)#fire it up!
library(TreePar)


plot(tree, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="yellow");
plot(tree, cex=0.2); tiplabels(cex=0.2, frame="circle", bg="yellow");
oz.mar<-drop.tip(tree, c(1:47)) #remove non-Australian marsupials
plot(oz.mar, cex=0.2); nodelabels(cex=0.2, frame="circle", bg="yellow");
write.tree(oz.mar, file="oz.marsupials.tre")

tree<-drop.tip(tree, c("Chelydra_serpentina", "Podocnemis_expansa", "Gallus_gallus", "Dromaius_novaehollandiae", "Crocodylus_porosus", "Alligator_mississippiensis", "Sphenodon_punctatus"))
write.tree(tree, file="squamates.tre")

## Summarize info from inter-node intervals 
tree<-read.tree("chiroptera.tre")
gammaStat(tree)
2*(1 - pnorm(abs(gammaStat(tree))))
gamStat(branching.times(tree), return.list=TRUE)

#### Diversity dependence estimate from LASER ####
#### Marsupials
mar<-getBtimes(file="oz.marsupials.tre") #using the 10% tree to test for artefactual shifts (<10 Ma)
result<-fitdAICrc(mar, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=170)
y2r<-yule2rate(mar, ints=120)
y3r<-yule3rate(mar, ints=120) ## best model for Oz marsupial rates, as decided by LASER
y4r<-yule4rate(mar, ints=120)
y5r<-yule5rate(mar, ints=120)
lh.comp.2.3<-pchisq(2*(y2r[1]-y3r[1]),3)
lh.comp.3.4<-pchisq(2*(y3r[1]-y4r[1]),3)


###############################################
#### TESS run of Pygopodoidea ######
################################################
tree<-read.tree("old.radiation.tre")
times<-branching.times(tree)

#### shortest possible version with less parameters, and most estimated by the package
#I've commented out a number of the parameters, because they overlap (i.e. "estimateNumberRateChanges" and "numExpectedRateChanges"), just uncomment and switch if you want to specify values
#this version will estimate most parameters, and as such, might be a good place to start
tess.analysis(tree=tree,
              samplingProbability = 1, #extant sampling probability at current time (sampling fraction)
              empiricalHyperPriors = TRUE, #should the hyperpriors be estimated empirically?
              estimateNumberRateChanges = TRUE, #should the number of rate shifts be estimated empirically?
              #numExpectedRateChanges = log(2), #expected number of rate changes following Poisson proces. at default 0 changes = 0.5 (50%) probability
              #pMassExtinctionPriorShape1 = 2, #alpha (first shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event
              #pMassExtinctionPriorShape2 = 10, #beta (second shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event
              estimateNumberMassExtinctions = TRUE, #should we estimate the number of mass-extinction events?
              #numExpectedMassExtinctions = log(2), #expected number of extinction events following Poisson proces. at default 0 changes = 0.5 (50%) probability
              estimateMassExtinctionTimes = TRUE, #should we estimate the times of mass-extinction events?
              #tInitialMassExtinction = c(30,35), #the intial value of the vector or times of the mass-extinction events, (initial MCMC values)
              MRCA = TRUE, #does the process start with the most recent common ancestor? must have root edge!
              CONDITION = "survival", #should we condition the process on "time", "survival" or "taxa"?
              BURNIN = 10000, #burnin quantity
              MAX_ITERATIONS = 10000000, #number of MCMC iterations
              THINNING = 10000, #recording frequency during the MCMC search
              OPTIMIZATION_FREQUENCY = 500, #how frequently the MCMC moves optimization
              CONVERGENCE_FREQUENCY = 1000, #how frequently we check for convergence
              MAX_TIME = Inf, MIN_ESS = 500, #maximum time the MCMC can run #the minimum ESS to assume convergence
              ADAPTIVE = TRUE, #auto-tuning of MCMC?
              dir = "" , #subdirectory for the output files, default ("") is the present directory
              priorOnly = FALSE, #sample from the prior only?
              verbose = TRUE) #want detailed output?

#### Summarizing the output of a diversification rate estimation including mass-extinction events
yrOutput <- tess.process.output(dir=getwd(),
                                 tree=(tree),
                                 criticalBayesFactors=c(2,6,10),
                                 burnin=.25)
#numExpectedRateChanges=log(3), #uncomment these two lines if you give it an expected number of extinction events
#numExpectedMassExtinctions=log(3))

#### Plotting the output of a diversification rate estimation including mass-extinction events
figures<-tess.plot.output(yrOutput,
                          fig.types=c("speciation rates",
                                      "speciation shift times",
                                      "speciation Bayes factors",
                                      "extinction rates",
                                      "extinction shift times",
                                      "extinction Bayes factors",
                                      "net-diversification rates",
                                      "relative-extinction rates",
                                      "mass extinction times",
                                      "mass extinction Bayes factors"),
                          xlab="million years ago",
                          col=NULL,
                          col.alpha=50,
                          xaxt="n",
                          yaxt="s",
                          pch=19,
                          plot.tree=TRUE)


###############################################
#### Condensed version of TreePar for OZ Marsupials
###############################################
tree<-read.tree("Squamates.tre")
y<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
p.noshift<-bd.shifts.optim(y, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
#use miniall=(previous) to give it the sampling from your previous run and save time!
p.oneshift<-bd.shifts.optim(y, sampling=c(1,1), miniall=p.noshift, survival=1, grid=1, start=0, end=205)[[2]]
p.twoshift<-bd.shifts.optim(y, sampling=c(1,1,1), miniall=p.oneshift, survival=1, grid=1, start=0, end=205)[[2]]
p.threeshift<-bd.shifts.optim(y, sampling=c(1,1,1,1), miniall=p.twoshift, survival=1, grid=1, start=0, end=205)[[2]]
p.fourshift<-bd.shifts.optim(y, sampling=c(1,1,1,1,1), miniall=p.threeshift, survival=1, grid=1, start=0, end=205)[[2]]
p.fiveshift<-bd.shifts.optim(y, sampling=c(1,1,1,1,1,1), miniall=p.fourshift, survival=1, grid=1, start=0, end=205)[[2]]
p.sixshift<-bd.shifts.optim(y, sampling=c(1,1,1,1,1,1,1), miniall=p.fiveshift, survival=1, grid=1, start=0, end=205)[[2]]
i<-1
comparison.0.1<-pchisq(2*(p.fourshift[[i]][1] - p.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(p.fourshift[[i+1]][1]-p.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(p.fourshift[[i+2]][1]-p.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(p.fourshift[[i+3]][1]-p.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(p.fourshift[[i+4]][1]-p.fourshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(p.fourshift[[i+5]][1]-p.fourshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
p.fourshift
comparison.0.1<-pchisq((p.sixshift[[i]][1]-p.sixshift[[i+1]][1]),1) # trying to determine if I should be using a single degree of freedom , it's just a comparison of two models
comparison.1.2<-pchisq((p.sixshift[[i+1]][1]-p.sixshift[[i+2]][1]),1) #but it should still stay as 3 DF because there are 4 parameters (2 lambda, 2 mu)

# Plot parameter estimates
###### you have to do all the "p.x.shift" options above if you want to plot them all together
###### if you just want to compare the values to get the best configuration, you can go directly to the highest number of shifts
plot0<-bd.shifts.plot(list(p.noshift),0,67,0,0.15)
par(new=T)
plot1<-bd.shifts.plot(list(p.oneshift),1,67,0,0.15)
par(new=T)
plot2<-bd.shifts.plot(list(p.twoshift),2,67,0,0.15)
par(new=T)
plot3<-bd.shifts.plot(list(p.threeshift),3,67,0,0.15)
par(new=T)
plot4<-bd.shifts.plot(list(p.fourshift),4,67,0,0.15)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(p.oneshift),1,25,0,0.2, plotturnover=TRUE)


#######################################
#### LTT plots and Null PB models
#######################################
# Create and Plot the Lineage Through Time
chi<-read.tree("chiroptera.tre")
ltt.plot(chi, log="y") #plot the LTT, check tree to find out # of tips
lines(c(-57.6,0), c(2,812), lty=2) #give the age (depth of tree), and the number of taxa

# Create and Plot the Null Pure Birth model simulations, including a 95% confidence interval
chipbtree<-pbtree(n=812, scale=57.6, nsim=1000, extant.only=TRUE) #5000 sims, based on age and number of taxa
ltt(chipbtree) #plot the 1000 simulations
ltt95(pbtree, log.lineages=TRUE)
plot(ltt95(chipbtree, log=TRUE), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped
par(new=T) #keep our simulation plot in place, while we drop on the 
ltt.plot(chi, log="y") #drop the empirical LTT on top

#####################################################################
#### Plot Distribution of Speciation Events (not all branching events)
old<-read.tree("old.radiation.tre")
n<-length(old$tip.label)
ors<-setNames(old$edge.length[sapply(1:n,function(x,y)   which(y==x),y=old$edge[,2])],old$tip.label)
hist<-hist(ors, breaks=30, xlab="Branching Times", col="lightpink", ylim=c(0,50), xlim=c(30,0)) #breaks determines the # of bins distributed across the whole xlim=
multiplier<-hist$counts/hist$density
mydensity<-density(ors) #pull the speciation frequencies out
mydensity$y<-mydensity$y*multiplier[1]
lines(mydensity) #plot the smoothed-out line of best fit across our histogram
abline(v=mean(ms), col="blue", lwd=2) #add a line for the mean
abline(v=median(ors), col="red", lwd=2) #add a line for the median


