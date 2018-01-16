install.packages('TreeSimGM')
library(PBD)
library(BAMMtools)
library(laser)
library(BioGeoBEARS)


data(warblers)
plotLtt(warblers)
ltt.plot(warblers)

################################################
#### Assessing Protracted Speciation via the Pure Birth Death Model
################################################

tree<-read.tree("pygopodoidea.tre")
btimes<-branching.times(tree)
pbd_ML(btimes)
pbd_loglik(pars1=c(0.133, 0.0, 0.081, 0.00), brts=btimes) #returns the log likelihood based on the estimates from our ML search
pbd_durspec_mean(pars=c(0.133, 0.081, 0.00)) #compute the mean duration to speciation value, value of 0.53 means 5.3 million years
pbd_durspec_quantile(pars=c(0.133, 0.081, 0.00), 0.05) #can give you the confidence intervals for the mean, must do twice, at 0.95 and 0.05
pbd_durspec_quantile(pars=c(0.133, 0.081, 0.00), 0.95) #can give you the confidence intervals for the mean, must do twice, at 0.95 and 0.05
fitdAICrc(btimes, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=600) #100 times tree depth
yule4rate(btimes, ints=600)
yule5rate(btimes, ints=600)

tree<-read.tree("carphodactylidae.tre")
btimes<-branching.times(tree)
pbd_ML(btimes)

tree<-read.tree("core.diplodactylidae.tre")
btimes<-branching.times(tree)
pbd_ML(btimes)
tree<-read.tree("Pygopodidae.tre")
btimes<-branching.times(tree)
pbd_ML(btimes)

tree<-read.tree("oz.marsupials.tre")
btimes<-branching.times(tree)
pbd_ML(btimes) #maximum likelihood estimate based on the branching times
pbd_loglik(pars1=c(0.156, 0.109, 1.916, 0.10), brts=btimes) #returns the log likelihood based on the estimates from our ML search
pbd_durspec_mean(pars=c(0.156, 0.109, 1.916, 0.10)) #compute the mean duration to speciation value, value of 0.53 means 5.3 million years
pbd_durspec_quantile(pars=c(0.156, 0.109, 1.916, 0.10), 0.05) #can give you the confidence intervals for the mean, must do twice, at 0.95 and 0.05
pbd_sim_cpp(pars=c(0.156, 0.109, 1.916, 0.10), age=66.6, soc=1, plotltt=1)
par(new=T)
ltt.plot(tree, log="y", col="red")
marpbtree<-pbtree(n=146, scale=66.7, nsim=100, extant.only=TRUE) #1000 sims, based on age and number of taxa
par(new=T)
plot(ltt95(marpbtree, log=TRUE, col="red"), xaxis="flipped") #plot the 95% confidence intervals with the y axis log transformed and the x axis flipped


tree<-read.tree("oz.crayfish.tre")
btimes<-branching.times(tree)
pbd_ML(btimes, btorph=0) #maximum likelihood estimate based on the branching times, returns the values necessary below
x<-pbd_loglik(pars1=c(0.133, 0.000, 0.081, 0.000), pars2=c(0, 0, 2, 1, "lsoda", 0, 0), brts=btimes, missnumspec=0) #returns the log likelihood based on the estimates from our ML search
#logLik = 105.5848
fitdAICrc(btimes, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=570)
pbd_durspec_mean(pars=c(0.133, 0.000, 0.081, 0.000)) #compute the mean duration to speciation value, value of 0.53 means 5.3 million years
pbd_durspec_quantile(pars=c(0.133, 0.000, 0.081, 0.000), 0.05) #can give you the confidence intervals for the mean, must do twice, at 0.95 and 0.05
?logLik()

pbd_sim_cpp(pars=c(0.133, 0.000, 0.081, 0.000), age=57, soc=1, plotltt=1)
par(new=T) #keep our simulation plot in place, while we drop on the 
ltt.plot(tree, log="y", col="red") #drop the empirical LTT on top

tree<-read.tree("primates.tre")
btimes<-branching.times(tree)
missnumspec=9
pars2=c(1,1,0)
endmc=1000
MLout=pbd_ML(btimes, initparsopt=c(0.66,0.63,0.99), idparsopt=1:3,exteq=1, missnumspec=missnumspec)
MLpars=as.numeric(unlist(MLout[1:4]))
exp_durspec=pbd_durspec_mean(c(MLpars[1],MLpars[4],MLpars[3]))
median_durspec = pbd_durspec_quantile(c(MLpars[1],MLpars[3],MLpars[4]),0.5)

pbd_ML(btimes, btorph=0)
pbd_durspec_mean(c(15.019, 0.1207, 14.996))
pbd_durspec_quantile(c(15.019, 0.1207, 14.996), 0.05)
x<-pbd_loglik(pars1=c(0.114841, 0.00000, 12.893589, 0.0000), pars2=c(0, 0, 2, 1, "lsoda", 0, 0), brts=btimes, missnumspec=0)
#logLik = 6499.911
fitdAICrc(btimes, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=278) #100 times tree depth

sum(btimes/125)
calc_AIC_vals(105.588444, 4)

#####################################################################
#### Plot Distribution of Speciation Events (not all branching events)
#####################################################################
beet<-read.tree("oz.diving.beetles.tre") #126 tips, median TTS 1.58
n<-length(beet$tip.label)
bs<-setNames(beet$edge.length[sapply(1:n,function(x,y)   which(y==x),y=beet$edge[,2])],beet$tip.label)
hist<-hist(bs, breaks=25, xlab="Branching Times", col="lightpink", ylim=c(0,20), xlim=c(25,0)) #breaks determines the # of bins distributed across the whole xlim=
multiplier<-hist$counts/hist$density
mydensity<-density(bs) #pull the speciation frequencies out
mydensity$y<-mydensity$y*multiplier[1]
lines(mydensity) #plot the smoothed-out line of best fit across our histogram
abline(v=mean(bs), col="blue", lwd=2) #add a line for the mean
abline(v=median(bs), col="red", lwd=2) #add a line for the median
median(bs)
mean(bs) #gives the time to speciation of just the tips, no deeper branches
btimes<-branching.times(beet)
sum(btimes)/116
