library(lattice)
library(coda)
library(btw)

################################################################
# Check convergence of MCMC in Coda
################################################################

agamids<-read.csv("Agamids.logSVL.csv", header=T) #this is a BayesTraits VarRates output
res <- mcmc(agamids,
            start=1, #define your starting row
            end=1050) #define your ending row, you can also thin the sampling
traceplot(res[,7]) #trace the variables to visualize mixing
ess(res) #check to make sure all ESS > 200
densityplot(res[,7])
summary(res)
HPDinterval(res)

