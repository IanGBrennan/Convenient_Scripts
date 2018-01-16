install.packages('rwty')
library(rwty)
x<-citation('rwty')
toBibtex(x)
####################################
# Check for Convergence using RWTY #
####################################

sptrees<-load.trees("MSY.spheno.extended.trees", type="nexus", format="beast", logfile="MSY.3clade.extended.log", skip=0) #read in your *beast trees
#older beast files have fewer skipped lines, so designate the header as the first line with "skip=0"
approx.ess<-topological.approx.ess(sptrees, burnin=52) #check the ESS of the topology, we want at least >200
allplots<-analyze.rwty(sptrees, burnin=2000, window.size=100, treedist='RF') #analyze the trees and designate burnin
#the 'analyze.rwty' function above will tell you all the plots it's creating
#but we've put em into the 'allplots' object, because we don't want it to 
#print every single one of them (maybe you do, but I don't)


#if you get an error: "All MCMC chains must be the same length as their associated p tables"
#do a quick check to see where the discrepancy is:
length(xtrees$ptable[,1])
length(xtrees$trees)
#write.csv(sptrees$ptable, file="ptable.200genes")

# So, instead let's just have a look through the plots that might 
#actually mean something, most will have two plots (a trace and a distribution)
#for all parameters, make sure ESS>200 to determine sufficient MCMC mixing
##########################################################################
allplots$posterior.trace #plot the estimate of the posterior
allplots$likelihood.trace #plot the estimate of the likelihood
allplots$prior.trace #plot the estimate of the prior
allplots$speciescoalescent.trace #check the estimate of gene coalescence with "species"
allplots$topology.trace.plot #check the estimate of topological convergence
allplots$Chain.1.correlations #you've worked hard, have a pretty graph of the parameter estimates through your chain
allplots$BirthDeath.t.Species.trace #get a handle on the Birth/Death ratios
allplots$birthRate2.t.Species.trace #look a the Birth (speciation) rate estimate
allplots$relativeDeathRate2.t.Species.trace #try to get an idea of Death (extinction) rate
allplots$treespace.heatmap


###################################################################
# Explanations of these plots can be found in Dan Warren's tutorial. 
# You can find this information available here: https://github.com/danlwarren/RWTY
#################################################################################


