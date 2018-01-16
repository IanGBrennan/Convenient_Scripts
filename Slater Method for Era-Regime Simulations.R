library(mvtnorm)

#######################################################################################
# This is Graham's original script, transforms the second half of the VCV matrix by the
# alpha parameter FIRST, then applies the rate scalars (sigma.sq) AFTERWARDS. I'm not
# sure this works correctly, because it estimates nonsensical values of alpha.
#######################################################################################

# Using my approach, we just split the phylogenetic variance co-variance matrix 
## into two and transform the elements of the OU process one before adding them 
## back together. Use the function split.vcv() to split the matrix at your 
## desired time before present.
#phy <- read.nexus("Zheng.Wiens.squamates.tre") # test a really big tree (>2k tips)
phy <- read.nexus("BT.Sphenomorphines.tre")
#phy <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")
m <- split.matrices <- split.vcv(phy, 10)

# then transform the elements of the matrix evolving under an OU process with 
# the appropriate alpha and sigmasq. There’s a function in here called ouMatrix 
# that was written by Luke Harmon and used to be part of geiger. You can use that 
# to transform it. So if the second half of the tree evolves under an OU 
# process with alpha = 1, then
m[[2]] <- ouMatrix(m[[2]], alpha=2)

# you now have two time slice matrices, one with variances and covariances 
# transformed according to a stationary peak OU process. To get a process 
# analogous to a reversed “release” model from my 2013 paper, multiply both 
# by the same sigma.sq and add them back together
m <- lapply(m, function(x) x*0.1)
m.rev.rel <- m[[1]] + m[[2]]

# if you want a reverse release radiate, just multiply each by a 
# different rate scalar and add them
m[[1]] <- m[[1]] * 01
m[[2]] <- m[[2]] * 5
m.rev.rel.rad <- m[[1]] + m[[2]]

# then draw your simulated data from a multivariate normal distribution, 
# choosing an appropriate root state and using the desired vcv.
x <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))

x.ouwie <- NULL
x.ouwie <- as.data.frame(names(x))
x.ouwie[,2] <- as.data.frame(t(x))
x.geiger <- x.ouwie; names <- x.ouwie[,1]
rownames(x.geiger) <- names; x.geiger[,1] <- NULL


trc <- fitContinuous_paleo(phy, x.geiger, model="TRC", shift.time=10)
r.c <- fitContinuous_paleo(phy, x.geiger, model="radiate.constrain", shift.time=10)

OUwie.slice(phy, x.ouwie, model="OUMVA", root.station=T, timeslices=NA)


diag(vcv(phy))
diag(vcv(x))














#######################################################################################
# This is MY adjusted script, which applies the rate scalars (sigma.sq) to the two split
# VCV matrices FIRST, then transforms the second matrix under an OU process by applying
# the alpha parameter. I think this makes more sense.
#######################################################################################

# Using my approach, we just split the phylogenetic variance co-variance matrix 
## into two and transform the elements of the OU process one before adding them 
## back together. Use the function split.vcv() to split the matrix at your 
## desired time before present.
#phy <- read.nexus("Zheng.Wiens.squamates.tre") # test a really big tree (>2k tips)
phy <- read.nexus("BT.Sphenomorphines.tre")
#phy <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")
m <- split.matrices <- split.vcv(phy, 10)

# if you want a reverse release radiate, just multiply each by a 
# different rate scalar
m[[1]] <- m[[1]] * 01
m[[2]] <- m[[2]] * 3

# then transform the elements of the matrix evolving under an OU process with 
# the appropriate alpha and sigmasq. There’s a function in here called ouMatrix 
# that was written by Luke Harmon and used to be part of geiger. You can use that 
# to transform it. So if the second half of the tree evolves under an OU 
# process with alpha = 1, then
m[[2]] <- ouMatrix(m[[2]], alpha=0.2)

# add them back together
m.rev.rel.rad <- m[[1]] + m[[2]]

# then draw your simulated data from a multivariate normal distribution, 
# choosing an appropriate root state and using the desired vcv.
x <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))

x.ouwie <- NULL
x.ouwie <- as.data.frame(names(x))
x.ouwie[,2] <- as.data.frame(t(x))
x.geiger <- x.ouwie; names <- x.ouwie[,1]
rownames(x.geiger) <- names; x.geiger[,1] <- NULL


trc <- fitContinuous_paleo(phy, x.geiger, model="TRC", shift.time=10)
r.c <- fitContinuous_paleo(phy, x.geiger, model="radiate.constrain", shift.time=10)

