library(laser)
library(DDD)
library(TreePar)

tree<-read.nexus("YOURTREE") #read in your tree
bt<-branching.times(tree) #pull out the branching times (some programs want the phylo object, others want b-times)

### Check the Gamma statistic for the tree(s)
##################################
# x>0 suggests more branching towards the tips
# x<0 suggests more branching deeper in the tree, slowing down towards present
gamStat(bt)


### Rate Estimates Using LASER
##################################
LASER.ddl<-DDL(bt) # you can individually check the model estimates one by one (this is linear diversity-dependent)
LASER.ddx<-DDX(bt) # here again, this time it's exponential diversity-dependence
# or you can bunch all the "constant rate" models together
fitdAICrc(bt, modelset=c("pureBirth", "bd", "DDL", "DDX"), ints=NULL)
# if you're curious, you can also investigate temporal shifts in the yule model speciation rate
# you could use these to compare against multi-rate models from TreePar
yule2rate(bt, ints=NULL)
yule3rate(bt, ints=NULL) 
yule4rate(bt, ints=NULL)
yule5rate(bt, ints=NULL)


### Birth Death as measured by APE
##################################
APE.bd<-birthdeath(tree)
APE.bd$para[2] #technically this is the diversification rate (lambda-mu), but mu is so small it's negligible

#compare the speciation rate estimate of LASER to APE
LASER.bd<-bd(bt)
LASER.bd$r1; APE.bd$para[2] #samesies?


### Birth Death as measured by DDD
##################################
DDD.bd<-bd_ML(bt)

LASER.bd$r1; APE.bd$para[2]; DDD.bd$lambda0 #all similar still?


### Density Dependence as measured by DDD
##################################
DDD.dd<-dd_ML(bt) # what about getting the density dependent measurements?

# are the rates the same as measured by LASER?
LASER.ddl$r1; DDD.dd$lambda #samesies? close at least?


### BirthDeath specation/extinction rates as estimated by TreePar
##################################
TreePar.bd<-bd.shifts.optim(bt, sampling=1, survival=1) #this will give you the rates with 0 shifts

#compare these to the estimates from LASER, APE, & DDD
LASER.bd$r1; APE.bd$para[2]; DDD.bd$lambda0;  TreePar.bd[[2]] #all similar still?


#alternatively, you can try the bd.densdep.optim which will estimate the rates under a DD model in TreePar
#this might require some further tweaking/reading
bd.densdep.optim(x=bt, continuous=T, lambdainit=0.16, muinit=0.0002, Kinit=209, model=-1)


