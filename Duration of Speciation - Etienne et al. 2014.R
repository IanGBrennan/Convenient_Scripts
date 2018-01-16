library(PBD)

tree<-read.tree("pygopodoidea.tre")

##############################################
#### Estimating the Duration of Speciation using Etienne et al. 2014 (updated from E&R 2012)
##############################################
btimes<-branching.times(tree) #pull out the branch lengths
brts = .1*rev(sort(as.numeric(branching.times(tree))))
missnumspec = 0#input the estimated number of missing species
pars2 = c(1,1,0) #apply pars2 parameters
endmc = 1000 #direct the number of iterations
btimesML<-pbd_ML(btimes) #get a quick and dirty estimate of the starting parameters
brtsML<-pbd_ML(brts) #get a quick and dirty estimate of the starting parameters based on the btimes*10
MLout = pbd_ML(btimes,initparsopt = c(0.0744,0.0255,1.72),idparsopt = 1:3,exteq = 1,missnumspec = missnumspec)
MLout = pbd_ML(brts,initparsopt = c(0.0744,0.0255,1.72),idparsopt = 1:3,exteq = 1,missnumspec = missnumspec)
MLpars = as.numeric(unlist(MLout[1:4]))
exp_durspec = pbd_durspec_mean(c(MLpars[1],MLpars[4],MLpars[3]))
median_durspec = pbd_durspec_quantile(c(MLpars[1],MLpars[3],MLpars[4]),0.5)
mode_durspec = pbd_durspec_mode(c(MLpars[1],MLpars[3],MLpars[4]))
sum(btimes)/36


####################################################
#compare the loglikelihood value against LASER estimates
####################################################
calc_AIC_vals(btimesML[5], 4) #calculate the AIC value from the loglikelihood of the btimesML from PBD, number at end refers to # of parameters
calc_AIC_vals(-519.54,4)
result<-fitdAICrc(btimes, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=570)

tree<-read.tree("primates.tre")
brts = 100*rev(sort(as.numeric(branching.times(tree))))
btimes<-branching.times(tree)
btimesML<-pbd_ML(btimes)
brtsML<-pbd_ML(brts) #get a quick and dirty estimate of the starting parameters based on the btimes*10
result<-fitdAICrc(btimes, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=100)
calc_AIC_vals(btimesML[5],4)
calc_AIC_vals(brtsML[5],4)
y4r<-yule4rate(btimes, ints=100)
y5r<-yule5rate(btimes, ints=100)
