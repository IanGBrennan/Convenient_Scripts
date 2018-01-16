library(expoTree)
library(TreePar)
library(laser)



comparison.2.3<-pchisq(2*(y2r[5]-y3r[7]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.3.4<-pchisq(2*(y3r[7]-y4r[9]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.4.5<-pchisq(2*(y4r[9]-y5r[11]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)

tree.three<--(threeshift[[4]][1])+sum(log(2:length(x)))

?fitdAICrc

tree<-read.tree("Pygopodoidea.tre")
plot(tree)

x<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
x
shift<-LikShifts()

#### Condensed version for the 10% Pygopodoidea tree
tree<-read.tree("Pygopodoidea.ten.percent.tre")
plot(tree)
l<-sort(branching.times(tree),decreasing=TRUE)
l.noshift<-bd.shifts.optim(l, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
#use miniall=(previous) to give it the sampling from your previous run and save time!
l.oneshift<-bd.shifts.optim(l, sampling=c(1,1), survival=1, grid=1, start=0, end=60)[[2]]
l.twoshift<-bd.shifts.optim(l, sampling=c(1,1,1), survival=1, grid=1, start=0, end=60)[[2]]
l.threeshift<-bd.shifts.optim(l, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=60)[[2]]
l.fourshift<-bd.shifts.optim(l, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=60)[[2]]
l.fiveshift<-bd.shifts.optim(l, sampling=c(1,1,1,1,1,1), miniall=l.fourshift,survival=1, grid=1, start=0, end=60)[[2]]
l.sixshift<-bd.shifts.optim(l, sampling=c(1,1,1,1,1,1,1), miniall=l.fiveshift, survival=1, grid=1, start=0, end=60)[[2]]
i<-1
comparison.0.1<-pchisq(2*(p.fourshift[[i]][1]-p.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(p.fourshift[[i+1]][1]-p.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(p.fourshift[[i+2]][1]-p.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(p.fourshift[[i+3]][1]-p.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(p.sixshift[[i+4]][1]-p.sixshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(p.sixshift[[i+5]][1]-p.sixshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)

comparison.0.1<-pchisq((p.sixshift[[i]][1]-p.sixshift[[i+1]][1]),1) # trying to determine if I should be using a single degree of freedom , it's just a comparison of two models
comparison.1.2<-pchisq((p.sixshift[[i+1]][1]-p.sixshift[[i+2]][1]),1) #but it should still stay as 3 DF because there are 4 parameters (2 lambda, 2 mu)


# Read in functions for estimating the paramters
source("MEPP_Homework5_functions.R")

#### Likelihood score of a tree with constant speciation/extinction rates (conditioned on time only)
LikConstant(lambda, mu, sampling=1, x, root=0, survival=1) #lik=592.1059
LikConstant(lambda=0.08, mu=0, sampling=1, x, root=0, survival=1) #or without extinction (lik=591.2128)

#### Likelihood score of a tree with constant speciation/extinction rates (conditioned on time+taxa) ####
lambda=0.08
mu=0.003
LikConstantn(lambda, mu, sampling=1, x, root=0) #likelihood based on constant rates (lik=585.6097)
res<-LikConstantn(lambda=0.07, mu=0.003, sampling=1, x, root=1) #can designate lambda/mu specifically this way too, though this isn't better (lik=586.1533)
LikConstantn(lambda=0.08, mu=0, sampling=1, x, root=0) #or without extinction (lik=584.8835, better!)

#### Trying to get bd.shifts.optim to work... ####
bd.shifts.optim(x, sampling=c(1,1), posdiv=FALSE, survival=1, maxitk=5, start=0, grid=5, end=56)
respygo<-bd.shifts.optim(x, sampling=1, grid=1, start=30, end=56, survival=1, maxitk=5, posdiv=FALSE) [[2]]
i<-1
test<-pchisq(2*(res[[i]][1]-res[[i+1]][1]),3)
test
noshift<-bd.shifts.optim(x, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
oneshift<-bd.shifts.optim(x, sampling=c(1,1), miniall=noshift, survival=1, grid=1, start=0, end=60)[[2]]
twoshift<-bd.shifts.optim(x, sampling=c(1,1,1), miniall=oneshift, survival, grid=1, start=0, end=60)[[2]]
threeshift<-bd.shifts.optim(x, sampling=c(1,1,1,1), miniall=twoshift, survival, grid=1, start=0, end=60)[[2]]
fourshift<-bd.shifts.optim(x, sampling=c(1,1,1,1,1), miniall=threeshift, survival=1, grid=1, start=0, end=60)[[2]]
fiveshift<-bd.shifts.optim(x, sampling=c(1,1,1,1,1,1), miniall=fourshift, survival=1, grid=1, start=0, end=60)[[2]]
sixshift<-bd.shifts.optim(x, sampling=c(1,1,1,1,1,1,1), miniall=fiveshift, survival=1, grid=1, start=0, end=60)[[2]]


MEshift<-bd.shifts.optim(x, sampling=c(1, 0.2), survival=1, grid=1, start=0, end=60, ME=TRUE)[[2]]
MEshift2<-bd.shifts.optim(x, sampling=c(1, 0.2, 0.2), miniall=MEshift, survival=1, grid=1, start=0, end=60, ME=TRUE)[[2]]
MEshift3<-bd.shifts.optim(x, sampling=c(1, 0.2, 0.2, 0.2), miniall=MEshift2, survival=1, grid=1, start=0, end=60, ME=TRUE)[[2]]
MEshift4<-bd.shifts.optim(x, sampling=c(1, 0.2, 0.2, 0.2, 0.2), miniall=MEshift3, survival=1, grid=1, start=0, end=60, ME=TRUE)[[2]]

comparison.ME.0.1<-pchisq(2*(MEshift4[[i]][1]-MEshift4[[i+2]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.ME.1.2<-pchisq(2*(MEshift4[[i+1]][1]-MEshift4[[i+2]][1]),3)
comparison.ME.2.3<-pchisq(2*(MEshift4[[i+2]][1]-MEshift4[[i+3]][1]),3)
comparison.ME.3.4<-pchisq(2*(MEshift4[[i+3]][1]-MEshift4[[i+4]][1]),3)

plot0<-bd.shifts.plot(list(MEshift),1,60,0,1)
par(new=T)
plot1<-bd.shifts.plot(list(MEshift2),2,60,0,1) # (num.shifts ,x.range, , y.range)
par(new=T)
plot2<-bd.shifts.plot(list(MEshift3),3,60,0,1) # ( ,x.range, , y.range)
par(new=T)
plot2<-bd.shifts.plot(list(MEshift4),4,60,0,1) # ( ,x.range, , y.range)
par(new=T)



#you can skip the "noshift" option, because by doing the above "threeshift", you include 0 [[1]], 1 [[2]], 2 [[3]], and 3 [[4]] shifts
#or as many others as you want to test, this is controlled by the sampling fraction c(1,...) one unit per shift to estimate
#now let's compare the loglikelihoods to evaluate which model is best
comparison.0.1<-pchisq(2*(sixshift[[i]][1]-sixshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(sixshift[[i+1]][1]-sixshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(sixshift[[i+2]][1]-sixshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(sixshift[[i+3]][1]-sixshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(sixshift[[i+4]][1]-sixshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)

#above, note that [[i+...]] indicates the shift model to compare, and [1] refers to the first column of that model (loglikelihood)

AICnoshift<-2*5+2*noshift[[i]][1]
AIConeshift<-2*5+2*oneshift[[i+1]][1]
AICtwoshift<-2*5+2*twoshift[[i+2]][1]
AICthreeshift<-2*5+2*threeshift[[i+3]][1]
AICfourshift<-2*5+2*fourshift[[i+4]][1]
AICfiveshift<-2*5+2*fiveshift[[i+5]][1]
test0<-pchisq(2*(AIConeshift-AICnoshift),3)
test0<-pchisq(2*(AICnoshift-AIConeshift),3)
test1<-pchisq(2*(AICtwoshift-AIConeshift),3)
test1<-pchisq(2*(AIConeshift-AICtwoshift),3)
test2<-pchisq(2*(AICthreeshift-AICtwoshift),3)
test2<-pchisq(2*(AICtwoshift-AICthreeshift),3)
test3<-pchisq(2*(AICfourshift-AICthreeshift),3)
test3<-pchisq(2*(AICthreeshift-AICfourshift),3)
test4<-pchisq(2*(AICfiveshift-AICfourshift),3)
test4<-pchisq(2*(AICfourshift-AICfiveshift),3)
test<-pchisq(2*(AICthreeshift-AICfourshift),3)
test<-pchisq(2*(AICthreeshift-AICnoshift),3)


# Plot parameter estimates
plot4<-bd.shifts.plot(list(fourshift),4,60,0,0.15)
par(new=T)
plot3<-bd.shifts.plot(list(threeshift),3,60,0,0.15) # (num.shifts ,x.range, , y.range)
par(new=T)
plot2<-bd.shifts.plot(list(twoshift),2,60,0,0.15) # ( ,x.range, , y.range)
par(new=T)
plot1<-bd.shifts.plot(list(oneshift),1,60,0,0.15)
par(new=T)
plot0<-bd.shifts.plot(list(noshift),0,60,0,0.15)
par(new=T)
plot5<-bd.shifts.plot(list(fiveshift),5,60,0,0.15)
###### Most supported shift scheme with turnover plotted too
plot6<-bd.shifts.plot(list(threeshift),3,60,0,0.15, plotturnover=TRUE)

# Plot parameters used for simulation
lines(c(-max(x),-time[2]),c((lambda-mu)[2],(lambda-mu)[2]),col="red")
lines(c(-time[2],0),c((lambda-mu)[1],(lambda-mu)[1]),col="red")

#### MLE based on a Density Dependent Model ####
v<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
PygoDD<-bd.densdep.optim(v, discrete=TRUE, rho=1) [[2]]
resDD[[1]]

par<-c(0.1, 0.03, 0.1)
PygolikDD<-LikDD(par,model=-1, v, sampling=1, root=1)

#### Condensed version for Pygopodidae
tree<-read.tree("pygopodidae.tre")
y<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
p.noshift<-bd.shifts.optim(y, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
p.oneshift<-bd.shifts.optim(y, sampling=c(1,1), survival=1, grid=1, start=0, end=25)[[2]]
p.twoshift<-bd.shifts.optim(y, sampling=c(1,1,1), survival=1, grid=1, start=0, end=25)[[2]]
p.threeshift<-bd.shifts.optim(y, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=25)[[2]]
p.fourshift<-bd.shifts.optim(y, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=25)[[2]]
p.sixshift<-bd.shifts.optim(y, sampling=c(1,1,1,1,1,1,1), survival=1, grid=1, start=0, end=25)[[2]]
i<-1
comparison.0.1<-pchisq(2*(p.sixshift[[i]][1]-p.sixshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(p.sixshift[[i+1]][1]-p.sixshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(p.sixshift[[i+2]][1]-p.sixshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(p.sixshift[[i+3]][1]-p.sixshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(p.sixshift[[i+4]][1]-p.sixshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(p.sixshift[[i+5]][1]-p.sixshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)

comparison.0.1<-pchisq((p.sixshift[[i]][1]-p.sixshift[[i+1]][1]),1) # trying to determine if I should be using a single degree of freedom , it's just a comparison of two models
comparison.1.2<-pchisq((p.sixshift[[i+1]][1]-p.sixshift[[i+2]][1]),1) #but it should still stay as 3 DF because there are 4 parameters (2 lambda, 2 mu)

# Plot parameter estimates
plot0<-bd.shifts.plot(list(p.noshift),0,25,0,0.4)
par(new=T)
plot1<-bd.shifts.plot(list(p.oneshift),1,25,0,0.4)
par(new=T)
plot2<-bd.shifts.plot(list(p.twoshift),2,25,0,0.4)
par(new=T)
plot3<-bd.shifts.plot(list(p.threeshift),3,25,0,0.4)
par(new=T)
plot4<-bd.shifts.plot(list(p.fourshift),4,25,0,0.4)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(p.oneshift),1,25,0,0.2, plotturnover=TRUE)

#### Condensed version for Carphodactylidae
tree<-read.tree("carphodactylidae.tre")
c<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
c.noshift<-bd.shifts.optim(c, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
c.oneshift<-bd.shifts.optim(c, sampling=c(1,1), survival=1, grid=1, start=0, end=30)[[2]]
c.twoshift<-bd.shifts.optim(c, sampling=c(1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
c.threeshift<-bd.shifts.optim(c, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
c.fourshift<-bd.shifts.optim(c, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
c.sixshift<-bd.shifts.optim(c, sampling=c(1,1,1,1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
i<-1
comparison.0.1<-pchisq(2*(c.fourshift[[i]][1]-c.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(c.fourshift[[i+1]][1]-c.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(c.fourshift[[i+2]][1]-c.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(c.fourshift[[i+3]][1]-c.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(c.fourshift[[i+4]][1]-c.fourshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(c.fourshift[[i+5]][1]-c.fourshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)


# Plot parameter estimates
plot0<-bd.shifts.plot(list(c.noshift),0,25,0,0.2)
par(new=T)
plot1<-bd.shifts.plot(list(c.oneshift),1,25,0,0.2)
par(new=T)
plot2<-bd.shifts.plot(list(c.twoshift),2,25,0,0.2)
par(new=T)
plot3<-bd.shifts.plot(list(c.threeshift),3,25,0,0.2)
par(new=T)
plot4<-bd.shifts.plot(list(c.fourshift),4,25,0,0.2)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(c.oneshift),1,25,0,0.1, plotturnover=TRUE)

#### Condensed version for Core Diplodactylidae
tree<-read.tree("core.diplodactylidae.tre")
cd<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
cd.noshift<-bd.shifts.optim(cd, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
cd.oneshift<-bd.shifts.optim(cd, sampling=c(1,1), survival=1, grid=1, start=0, end=30)[[2]]
cd.twoshift<-bd.shifts.optim(cd, sampling=c(1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
cd.threeshift<-bd.shifts.optim(cd, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
cd.fourshift<-bd.shifts.optim(cd, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
cd.fiveshift<-bd.shifts.optim(cd, sampling=c(1,1,1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
cd.sixshift<-bd.shifts.optim(cd, sampling=c(1,1,1,1,1,1,1), survival=1, grid=1, start=0, end=30)[[2]]
i<-1
comparison.0.1<-pchisq(2*(cd.fourshift[[i]][1]-cd.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(cd.fourshift[[i+1]][1]-cd.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(cd.fourshift[[i+2]][1]-cd.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(cd.fourshift[[i+3]][1]-cd.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(cd.fourshift[[i+4]][1]-cd.fourshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(cd.fourshift[[i+5]][1]-cd.fourshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)


# Plot parameter estimates
plot0<-bd.shifts.plot(list(cd.noshift),0,30,-.05,0.15)
par(new=T)
plot1<-bd.shifts.plot(list(cd.oneshift),1,30,-.05,0.15)
par(new=T)
plot2<-bd.shifts.plot(list(cd.twoshift),2,30,-.05,0.15)
par(new=T)
plot3<-bd.shifts.plot(list(cd.threeshift),3,30,-.05,0.15)
par(new=T)
plot4<-bd.shifts.plot(list(cd.fourshift),4,30,-.05,0.15)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(cd.twoshift),2,30,0,0.8, plotturnover=TRUE)

#### Condensed version for Carpho/Pygo
tree<-read.tree("Carphodactylidae.Pygopodidae.tre")
cp<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
cp.noshift<-bd.shifts.optim(cp, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
cp.oneshift<-bd.shifts.optim(cp, sampling=c(1,1), survival=1, grid=1, start=0, end=52)[[2]]
cp.twoshift<-bd.shifts.optim(cp, sampling=c(1,1,1), survival=1, grid=1, start=0, end=52)[[2]]
cp.threeshift<-bd.shifts.optim(cp, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=52)[[2]]
cp.fourshift<-bd.shifts.optim(cp, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=52)[[2]]
cp.sixshift<-bd.shifts.optim(cp, sampling=c(1,1,1,1,1,1,1), survival=1, grid=1, start=0, end=52)[[2]]
i<-1
comparison.0.1<-pchisq(2*(cp.fourshift[[i]][1]-cp.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(cp.fourshift[[i+1]][1]-cp.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(cp.fourshift[[i+2]][1]-cp.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(cp.fourshift[[i+3]][1]-cp.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(cp.fourshift[[i+4]][1]-cp.fourshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(cp.fourshift[[i+5]][1]-cp.fourshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)


# Plot parameter estimates
plot0<-bd.shifts.plot(list(cp.noshift),0,55,0,0.2)
par(new=T)
plot1<-bd.shifts.plot(list(cp.oneshift),1,55,0,0.2)
par(new=T)
plot2<-bd.shifts.plot(list(cp.twoshift),2,55,0,0.2)
par(new=T)
plot3<-bd.shifts.plot(list(cp.threeshift),3,55,0,0.2)
par(new=T)
plot4<-bd.shifts.plot(list(cp.fourshift),4,55,0,0.2)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(cp.twoshift),2,55,0,0.2, plotturnover=TRUE)

#### Condensed version for Diplodactylidae
tree<-read.tree("diplodactylidae.tre")
d<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
d.noshift<-bd.shifts.optim(d, sampling=c(1), survival=1)[[2]] 
#above estimates loglikelihood, turnover (extinction), and diversification (speciation-extinction)
d.oneshift<-bd.shifts.optim(d, sampling=c(1,1), survival=1, grid=1, start=0, end=46)[[2]]
d.twoshift<-bd.shifts.optim(d, sampling=c(1,1,1), survival=1, grid=1, start=0, end=46)[[2]]
d.threeshift<-bd.shifts.optim(d, sampling=c(1,1,1,1), survival=1, grid=1, start=0, end=46)[[2]]
d.fourshift<-bd.shifts.optim(d, sampling=c(1,1,1,1,1), survival=1, grid=1, start=0, end=46)[[2]]
d.sixshift<-bd.shifts.optim(d, sampling=c(1,1,1,1,1,1,1), survival=1, grid=1, start=0, end=46)[[2]]
i<-1
comparison.0.1<-pchisq(2*(d.fourshift[[i]][1]-d.fourshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(d.fourshift[[i+1]][1]-d.fourshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(d.fourshift[[i+2]][1]-d.fourshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.3.4<-pchisq(2*(d.fourshift[[i+3]][1]-d.fourshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.4.5<-pchisq(2*(d.fourshift[[i+4]][1]-d.fourshift[[i+5]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.5.6<-pchisq(2*(d.fourshift[[i+5]][1]-d.fourshift[[i+6]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)


# Plot parameter estimates
plot0<-bd.shifts.plot(list(d.noshift),0,46,0,0.12)
par(new=T)
plot1<-bd.shifts.plot(list(d.oneshift),1,46,0,0.12)
par(new=T)
plot2<-bd.shifts.plot(list(d.twoshift),2,46,0,0.12)
par(new=T)
plot3<-bd.shifts.plot(list(d.threeshift),3,46,0,0.12)
par(new=T)
plot4<-bd.shifts.plot(list(d.fourshift),4,46,0,0.12)
###### Most supported shift scheme with turnover plotted too
plot5<-bd.shifts.plot(list(d.threeshift),3,46,0,1, plotturnover=TRUE)

d.threeshift
k<-d.threeshift[[4]][4] #turnover (extinction/speciation)
d<-d.threeshift[[4]][5] #diversification (speciation-extinction)
d.lam<-d/(1-k) #should give use the speciation rate
d.mu<-(d*k)/(1-k) #should give us the extinction rate



#### Trying to get DD shifts to work...####
#[i] or [[1]] will tell you the log likelihood of the param
#[[1]] will tell  you the 
x<-sort(branching.times(tree),decreasing=TRUE) #pulls and sorts the branching times in the tree
resDD<-bd.densdep.optim(y, discrete=TRUE) [[2]]
resDD[[1]]


tree<-sim.rateshift.taxa(10,1,c(2,0.1),c(0,0.05),frac=c(1,1),times=time,complete=FALSE)
z<-sort(getx(tree[[1]]),decreasing=TRUE)
resTest<-bd.densdep.optim(z, discrete=FALSE, continuous=TRUE) [[2]]
resTest
resTestShifts<-bd.shifts.optim(z,c(rho,1),0.1,0.1,2.1)[[2]][[2]]
z
x
# Best model where AIC smallest
AICDD <- 2*3+2*resTest$value
AICShifts <- 2*5+2*resTestShifts[1]
AICDD
AICShifts



#### Plot of parameter estimates of best models across multiple clades ####
plot0<-bd.shifts.plot(list(p.oneshift),1,60,0,0.2)
par(new=T)
plot0<-bd.shifts.plot(list(c.oneshift),1,60,0,0.2)
par(new=T)
plot0<-bd.shifts.plot(list(cd.twoshift),2,60,0,0.2)
par(new=T)
plot0<-bd.shifts.plot(list(threeshift),3,60,0,0.2)
par(new=T)
plot0<-bd.shifts.plot(list(cp.twoshift),2,60,0,0.2)
par(new=T)
plot0<-bd.shifts.plot(list(d.threeshift),3,60,0,0.2)
par(new=T)












###################################################################
# Simulation example - time-dependent rates
#
# First we simulate a tree:
# Number of species
nspecies <- 100

# At time 1 in the past, we have a rate shift:
time <- c(0,1)
#=(t0,t1=t2) from slide 11 in Macroevolution lecture

# Half of the present day species are sampled (rho[1]=0.5):
rho <- c(0.5,1)

# speciation rates (between t[i],t[i+1] we have speciation rate lambda[i]):
lambda <- c(2,5)  
#=(lambda0,lambda1=lambda2) from slide 11 in Macroevolution lecture

# extinction rates (between t[i],t[i+1] we have extinction rate mu[i]):
mu <- c(1.5,0)
#=(mu0,mu1=mu2) from slide 11 in Macroevolution lecture

# Simulation of a tree:
tree<-sim.rateshift.taxa(nspecies,1,lambda,mu,frac=rho,times=time,complete=FALSE)

# Extracting the speciation times x:
q<-sort(getx(tree[[1]]),decreasing=TRUE)
q

# When estimating the rate shift times t based on branching times x, 
# we allow the shift times to be 0.6, 0.8, 1, 1.2, .. ,2.4:
start <- 0.4
end <- 2
grid <- 0.02

# We fix rho and estimate time, lambda, mu:
res <- bd.shifts.optim(q,c(1,1,1,1),grid,start,end)[[2]]
res
res[[2]][1]
res[[1]][1]
res[[3]][1]
res[[4]][1]
compare<-pchisq(2*(res[[i]][1]-res[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
compare2<-pchisq(2*(res[[i]][1]-res[[i+2]][1]),3)
compare3<-pchisq(2*(res[[i]][1]-res[[i+3]][1]),3)

# res[[2]] tells us about the maximum likelihood estimate given one rate shift:
# - log lik = 60.6940763.
# rate shift at time 1.0.
# turnover (extinction/speciation) = 0.68 more recent than 1.0,
#     and = 0.19 more ancestral than 1.0.
# net diversification (speciation-extinction) rate = 0.81 more recent than 1.0, 
#     and = 3.66 more ancestral than 1.0.

#values used for simulation:
mu/lambda
lambda-mu

#test if 1 shift explain the tree significantly better than 0 shifts:
#if test>0.95 then 1 shift is significantly better than 0 shifts at a 5% error
i<-1
test<-pchisq(2*(res[[i]][1]-res[[i+1]][1]),3)
test

#test if 2 shifts explain the tree significantly better than 1 shift:
i<-2
test<-pchisq(2*(res[[i]][1]-res[[i+1]][1]),3)
test

# Plot parameter estimates
bd.shifts.plot(list(res),1,2.1,0,5)
# Plot parameters used for simulation
lines(c(-max(x),-time[2]),c((lambda-mu)[2],(lambda-mu)[2]),col="red")
lines(c(-time[2],0),c((lambda-mu)[1],(lambda-mu)[1]),col="red")









sampling <- c(0.1, 1)
backsamp <- c(1,0.1)
me.res.T <- bd.shifts.optim(x, sampling=c(1,1), ME=T, all=F, start=27, grid=0.5, end=35)[[2]]
me.res.TA <- bd.shifts.optim(x, sampling=c(1,1), ME=T, all=T, start=27, grid=0.5, end=35)[[2]]

me.res.F <- bd.shifts.optim(x, sampling=c(1,0.1), ME=F, all=T, start=27, grid=0.5, end=35)[[2]]
alt.me <- bd.shifts.optim(x, sampling=c(1,1), ME=F, start=27, grid=0.5, end=35)[[2]]

l.noshift<-bd.shifts.optim(l, sampling=c(1), survival=1)[[2]] 
l.oneshift<-bd.shifts.optim(l, sampling=c(1,1), miniall=l.noshift, survival=1, grid=1, start=0, end=60)[[2]]
l.twoshift<-bd.shifts.optim(l, sampling=c(1,1,1), miniall=l.oneshift, survival=1, grid=1, start=0, end=60)[[2]]
l.threeshift<-bd.shifts.optim(l, sampling=c(1,1,1,1), miniall=l.twoshift, survival=1, grid=1, start=0, end=60)[[2]]

calc_AIC_vals((-me.res.F[[1]][1]), 5)
calc_AIC_vals((-alt.me[[1]][1]), 5)


i <- 1
comparison.me.0<-pchisq(2*(me.res.F[[2]][1]-alt.me[[2]][1]),5)
comparison.0.me<-pchisq(2*(alt.me[[2]][1]-me.res.F[[2]][1]),5)

comparison.0.1<-pchisq(2*(l.threeshift[[i]][1]-l.threeshift[[i+1]][1]),3) #compare noshift vs. oneshift via chi square (we can reject oneshift)
comparison.1.2<-pchisq(2*(l.threeshift[[i+1]][1]-l.threeshift[[i+2]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.3<-pchisq(2*(l.threeshift[[i+2]][1]-l.threeshift[[i+3]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)
comparison.2.me<-pchisq(2*(l.threeshift[[i+2]][1]-me.res[[2]][1]),3)
comparison.3.4<-pchisq(2*(l.threeshift[[i+3]][1]-l.threeshift[[i+4]][1]),3) #oneshift vs. twoshift (meaningless considering noshift is the best, but we can still reject twoshift)

plot0<-bd.shifts.plot(list(me.res),1,60,0,1, plotturnover = T)

bd.shifts.plot(list(me.res.T),1,60,0,1, plotturnover = T)
bd.shifts.plot(list(me.res.TA),1,60,0,1, plotturnover = T)
bd.shifts.plot(list(me.res.F),1,60,0,1, plotturnover = T)


bd.shifts.plot(list(backwards.me),1,60,0,1, plotturnover = T)
par(new=T)

#construct a vector of the likelihoods, named by the models
candidateModels <- c("MassEx.F"=-(me.res.F[[2]][1]),
                     "Shift30"=-(alt.me[[2]][1]))
                     #"NoShift"=-(l.noshift[[1]][1]),
                     #"OneShift"=-(l.oneshift[[2]][1]),
                     #"TwoShift"=-(l.twoshift[[3]][1]),
                     #"ThreeShift"=-(l.threeshift[[4]][1]))
                     
                     #"shiftBDME"=lik.me.shift)
#make a grid with all possible combinations of the models
LikelihoodGrid <- expand.grid(M0=names(candidateModels),
                              M1=names(candidateModels))
#add a column that is the 2lnBF for each pair of models
LikelihoodGrid$BF <- 2 * (candidateModels[LikelihoodGrid$M0] -
                            candidateModels[LikelihoodGrid$M1])
#sort the comparisons by their 2lnBF in descending order
LikelihoodGrid <- LikelihoodGrid[order(LikelihoodGrid$BF,
                                       decreasing=TRUE),]
LikelihoodGrid #look at the grid, determine which is the preferred model


-me.res[[1]][1]+sum(log(2:length(x)))
-l.noshift[[1]][1]+sum(log(2:length(x)))



tess.plot



trees <- read.tree("Pygopodoidea.TESS.loop.trees")
me.vs.shift <- NULL
for (i in 1:100) {
  x <- sort(branching.times(trees[[i]]),decreasing=TRUE) #pulls and sorts the branching times in the tree

  me.res <- bd.shifts.optim(x, sampling=c(1,0.1), ME=F, all=T, start=27, grid=0.5, end=34)[[2]]
  alt.me <- bd.shifts.optim(x, sampling=c(1,1), ME=F, start=27, grid=0.5, end=34)[[2]]
  
  #construct a vector of the likelihoods, named by the models
  candidateModels <- c("MassEx"=-(me.res[[2]][1]),
                       "Shift30"=-(alt.me[[2]][1]))
  
  #make a grid with all possible combinations of the models
  LikelihoodGrid <- expand.grid(M0=names(candidateModels),
                                M1=names(candidateModels))
  #add a column that is the 2lnBF for each pair of models
  LikelihoodGrid$BF <- 2 * (candidateModels[LikelihoodGrid$M0] -
                              candidateModels[LikelihoodGrid$M1])
  #sort the comparisons by their 2lnBF in descending order
  LikelihoodGrid <- LikelihoodGrid[order(LikelihoodGrid$BF,
                                         decreasing=TRUE),]
  me.vs.shift <- rbind(me.vs.shift, LikelihoodGrid[1,])
}

write.table(me.vs.shift, file="/Users/Ian/Google.Drive/R.Analyses/TESS/Analysis Loop/TreePar.Model.Comparison.table.xls", 
            sep="\t", quote=F, row.names=F) # write it all to a file
table(me.vs.shift[,1])
