library(devtools)
library(githubinstall)
library(phytools)

source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/get.landscape.BBMV_no_bounds.r',chdir=F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/ML_functions.r',chdir=F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/utils_BBMV.r',chdir=F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/Simulate\ BBM+V.r',chdir = F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/charac_time.r',chdir=F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/Uncertainty_BBMV.r',chdir=F)
source('/Users/Ian/Google.Drive/R.Analyses/BBMV-master/R/ACE_FPK.R',chdir=F)


library(geiger)
tree=sim.bdtree(stop='taxa',n=50)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x))
V6=10*(x^4-0.5*(x^2)+0.*x)
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size)
par(mfrow=c(1,1))
plot(V6_norm)
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=1,bounds=bounds)
hist(TRAIT,breaks=20)


#######################################################
# Read in all the data you'll use
######################################################
tree<-read.nexus("BT.Agamids.tre") #read in desired tree
#trees <- read.nexus("PB.Meliphagides.100.trees")
#drip <- c("S_Papuascincus_sp","S_Prasinohaema_virens","S_Sphenomorphus_jobiensis", "S_Sphenomorphus_muelleri","S_Sphenomorphus_solomonis")
#trees<-lapply(trees, drop.tip, tip=drip) #drop tips if necessary

#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data       <- read.csv("BT.Agamids.logSVL.csv", row.names = 1, header=F) #read in data file in GEIGER format
data.OUwie <- read.csv("BT.Agamids.logSVL.csv", header=F) #read in data file in OUwie format
data.fitenv<- data.OUwie[,2]; names(data.fitenv) <- data.OUwie[,1] #read in data file in RPANDA format

hist(data.fitenv, breaks=50)
bounds=c(min(data.fitenv),max(data.fitenv))


#logBL<-setNames(data[,14], rownames(data)) #choose your data column (here: logBL) and apply rownames
name.check(tree, data); name.check(tree, data.fitenv) #check to make sure the tips match the data labels

ll_FPK4=lnL_FPK(tree,data.fitenv,Npts=50,a=NULL,b=NULL,c=NULL) # the full model
ll_FPK2=lnL_FPK(tree,data.fitenv,Npts=25,a=0,b=NULL,c=NULL) # OU
ll_FPK0=lnL_FPK(tree,data.fitenv,Npts=25,a=0,b=0,c=0) # BM

fit4=find.mle_FPK(model=ll_FPK4)
get.landscape.FPK(fit=fit4)
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit2=find.mle_FPK(model=ll_FPK2)
get.landscape.FPK(fit=fit2) # this shape of the landscape cannot have 2 peaks
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit0=find.mle_FPK(model=ll_FPK0) 
get.landscape.FPK(fit=fit0) # this one is forced to be flat
lines(V6_norm~seq(from=min(bounds),to=max(bounds),length.out=length(V6_norm)))

fit4$aic 
fit2$aic
fit0$aic

charac_time(fit=fit4)
charac_time(fit=fit2)
charac_time(fit=fit0)
max(branching.times(tree))

Uncertainty_FPK(fit=fit4, tree, trait=data.fitenv, Npts=25, effort_uncertainty=100,
                scope_a=c(-1,10), scope_b=c(-5,5),scope_c=c(-2,2))

oufit    <- fitContinuous(tree, data, SE=NA, model="OU", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))
bmfit    <- fitContinuous(tree, data, SE=NA, model="BM", bounds=list(alpha=c((log(2)/max(nodeHeights(tree))), 10)))

ACE_nodes=ACE_FPK(fit4,specific.point=NULL)
plot(ACE_nodes[[90]],type='l')


ll_BBMV4=lnL_BBMV(tree,data.fitenv,Npts=25,bounds=bounds,a=NULL,b=NULL,c=NULL)
ll_BBMV2=lnL_BBMV(tree,data.fitenv,Npts=25,bounds=bounds,a=0,b=NULL,c=NULL)
ll_BBMV1=lnL_BBMV(tree,data.fitenv,Npts=25,bounds=bounds,a=0,b=0,c=NULL)
ll_BBMV0=lnL_BBMV(tree,data.fitenv,Npts=25,bounds=bounds,a=0,b=0,c=0) # this is the BBM model


fit4b=find.mle_FPK(model=ll_BBMV4)
get.landscape.FPK(fit=fit4b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit2b=find.mle_FPK(model=ll_BBMV2)
get.landscape.FPK(fit=fit2b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit1b=find.mle_FPK(model=ll_BBMV1)
get.landscape.FPK(fit=fit1b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))

fit0b=find.mle_FPK(model=ll_BBMV0)
get.landscape.FPK(fit=fit0b)
lines(Vb_norm~seq(from=min(bounds),to=max(bounds),length.out=length(Vb_norm)))










