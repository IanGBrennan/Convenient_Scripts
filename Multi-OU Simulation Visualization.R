library(OUwie)
library(randomcoloR)

# ==How to simulate Ornstein-Uhlenbeck ==

# Look at the code for the Ornstein-Uhlenbeck model in Functions_ContTraitSim.R. You can see it simply converts the equation into R code but also includes how many generations you want to run (time) and the starting value for the trait (x0). 

# time is the length of the simulation, theta is the optima, 
# alpha is the pull towards the optima, sigma is the rate and
# X0 is the starting trait value, you could combine the optima 
# and X0 but I find it helpful to separate them as you can also
# investigate how changes in parameter values change how quickly 
# the optima is reached
ornstein_uhlenbecksim <- function(time,theta,alpha,sigma2,x0){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma2*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  return(x);
}

OU.sim <- function(n, theta, alpha, sigma,x0){
  dw  <- rnorm(n, 0)
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma*dw[i-1]
  }
  return(x);
}

multiOU.sim <- function(time, theta, alpha, sigma, x0, theta1, alpha1, sigma1, shift){
  dw  <- rnorm(time, 0)
  x <- c(x0)
  for (i in 2:(time-shift)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma*dw[i-1]
  }
  for (i in ((time-shift)+1):(time+1)) {
    x[i]  <-  x[i-1] + alpha1*(theta1-x[i-1]) + sigma1*dw[i-1] 
  }
  x <- rev(x)
  return(x);
}



a.OU.sims <- replicate(10, ornstein_uhlenbecksim(time=100, theta=0.75, alpha=0.5, sigma=0.03, x0=0.25))
b.OU.sims <- replicate(10, OU.sim(n=100, theta=0.75, alpha=0.5, sigma=0.03, x0=0.25), simplify=FALSE)
plot(b.OU.sims[[1]], type="n", ylim=c(0,1), xlab="Time", ylab="Trait", main="alpha=0.5, sigma=0.03")
obj <- lapply(b.OU.sims, lines, col="blue")
c.OU.sims <- replicate(10, multiOU.sim(time=100,theta=0,alpha=0,sigma=0.03,x0=0,theta1=0,alpha1=0.5,sigma1=0.03,shift=25))
plot(c.OU.sims[[1]], type="n", ylim=c(0,1), xlab="Time", ylab="Trait", main="alpha=0.5, sigma=0.03")
obj <- lapply(c.OU.sims, lines, col="blue")
plot(c.OU.sims[,1], xlim=c(100,0),ylim=c(-1, 1), type="l", col="red", xlab="Time", ylab="Trait")
for(i in 2:100) lines(c.OU.sims[,i], col=randomColor(length(empiricalTRC[,1])))




#### New Multi-OU Simulator
multi.ornstein_uhlenbecksim <- function(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time-shift)) {
    x[i]  <-  x[i-1] + alpha0*(theta0-x[i-1]) + sigma0*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  for (i in ((time-shift)+1):(time+1)) {
    x[i]  <-  x[i-1] + alpha1*(theta1-x[i-1]) + sigma1*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  x <- rev(x)
  return(x);
  print(c(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift))
}


time=65; theta0=1.2; alpha0=0; sigma0=0.05; x0=1.2; theta1=1.2; alpha1=0.5; sigma1=0.0075; shift=7
birdsTRC <- replicate(100, multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift))
plot(birdsTRC[,1], xlim=c(time,0),ylim=c((x0-.3*x0), (x0+.3*x0)), type="l", col="red", xlab="Time", ylab="Trait")
for(i in 2:100) lines(birdsTRC[,i], col="red")





##### TRY FIXING THIS TOMORROW
multi.ornstein_uhlenbecksim <- function(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time-shift)) {
    x[i]  <-  x[i-1] + alpha0*(theta0-x[i-1]) + sigma0*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  for (i in ((time-shift)+1):(time+1)) {
    x[i]  <-  x[i-1] + alpha1*(theta1-x[i-1]) + sigma1*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  x <- rev(x)
  return(x);
}


theta0=1.2; theta1=1.2
alpha0=0; alpha1=0.5
sigma0=0.005; sigma1=0.0075
time=65; shift=7; x0=1.2
#TRC birds:  time=65; theta0=1.2; alpha0=0; sigma0=0.00578; x0=1.2; theta1=1.2; alpha1=0.458; sigma1=0.0773; shift=7
#SRC birds:  time=65; theta0=1.2; alpha0=0; sigma0=0.005; x0=1.2; theta1=1.2; alpha1=0.5; sigma1=0.005; shift=7
#null birds: time=65; theta0=1.2; alpha0=0; sigma0=0.005; x0=1.2; theta1=1.2; alpha1=0; sigma1=0.005; shift=7

#TRC agamid: time=35; theta0=2.1; alpha0=0; sigma0=0.025; x0=2.1; theta1=2.1; alpha1=1; sigma1=0.0025; shift=8
#TRC skinks: time=31; theta0=1.8; alpha0=0; sigma0=0.0025; x0=1.8; theta1=1.8; alpha1=0.5; sigma1=0.00125; shift=8

#test <- multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift)
#plot(test, xlim=c(100,0),ylim=c(-1, 7), type="l", col="red", xlab="Time", ylab="Trait")
empiricalTRC <- replicate(100, multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift))
plot(empiricalTRC[,1], xlim=c(time,0),ylim=c(0, 2), type="l", col="red", xlab="Time", ylab="Trait")
for(i in 2:100) lines(empiricalTRC[,i], col=randomColor(length(empiricalTRC[,1])))

empiricalSRC <- replicate(100, multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift))
plot(empiricalSRC[,1], xlim=c(65,0),ylim=c(1.5, 1.5), type="l", col="black", xlab="Time", ylab="Trait")
for(i in 2:1000) lines(empiricalSRC[,i], col="black")

nullBM <- replicate(100, multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift))
plot(nullBM[,1], xlim=c(65,0),ylim=c(1, 1.5), type="l", col="green", xlab="Time", ylab="Trait")
for(i in 2:100) lines(nullBM[,i], col="green")



#### Now with Empirical Data!
#############################
data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Parameter Estimates/Meliphagides.BEST.Model.PARAM.ESTIMATES.csv", header=T)
trc.data <- subset(data, model=="TRC")
trc.data$model.timing <- as.numeric(trc.data$model.timing)

nums <- unique(trc.data$tree.no)

trc.parameters <- NULL
params.used <- NULL
for (i in 1:length(nums)) {
  focus.data <- subset(trc.data, tree.no==nums[i])
  time = 31
  theta0 = focus.data$before.shift[3]
  alpha0 = 0
  sigma0 = (focus.data$before.shift[1])
  x0 = focus.data$before.shift[3]
  theta1 = focus.data$before.shift[3]
  alpha1 = focus.data$after.shift[2]
  #alpha1 = 0.79
  sigma1 = ((focus.data$after.shift[1])*(focus.data$before.shift[1])) # this is the BM sigma * the post-shift-scalar
  shift = focus.data$model.timing[1]
  params <- as.data.frame(t(c(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift)))
  params.used <- rbind(params.used, params)
  empiricalTRC <- multi.ornstein_uhlenbecksim(time,theta0,alpha0,sigma0,x0,theta1,alpha1,sigma1,shift)
  trc.parameters <- cbind(trc.parameters, empiricalTRC)
}
colnames(params.used) <- c("time","theta0","alpha0","sigma0","x0","theta1","alpha1","sigma1","shift")
plot(trc.parameters[,1], xlim=c(time,0),ylim=c(1.1,1.4), type="l", col="red", xlab="Time", ylab="Trait")
for(i in 2:length(nums)) lines(trc.parameters[,i], col=randomColor(length(trc.parameters[,1])))


focus.data <- subset(trc.data, tree.no==1)
class(focus.data$model.timing[3])

abline(a=0, b=0.005)


#######################################
# http://schmitzlab.info/BMandOU.html #
#######################################
#### Let's start with Brownian Motion
changes <- rnorm(1000) # sd=1
plot(changes, type="l", col="blue", ylab="Trait Change")
# the cumulative sum represents the path of the trait through time
sum.changes <- cumsum(changes)
plot(sum.changes, type="l", col="blue", xlab="Time", ylab="Trait Value")
# now let's compare against a change in the standard deviation (jump size)
heated.up <- rnorm(1000, sd=2) #new sd-value
plot(changes, type="l", col="blue", ylim=c(-6,6), ylab="Trait Change")
lines(heated.up, col="red")
# and the cumulative sum of this new pattern:
plot(cumsum(changes), type="l", col="blue", ylim=c(-80,80), xlab="Time", ylab="Trait Value")
lines(cumsum(heated.up), col="red")


# BM function where all we can change is the Standard Deviation
bm.sim0 <- function(sd){
  x <- rnorm(1000, sd=sd)
  sim <- cumsum(x)
  plot(sim, type="l", xlab="Time", ylab="Trait")
}
bm.sim0(2)

# BM function where we can change the SD, time, and mean
bm.sim <- function(sd, t=1000, x0=0, plot=FALSE, col="black"){
  sim <- cumsum(rnorm(t, mean=x0, sd=sd))
  if (plot==TRUE) {plot(sim, type="l", xlab="Time", ylab="Trait", col=col)}
  if (plot=="ADD") {lines(sim, col=col)}
  return(sim);
}
# test it out
bm.sim(2, 100, 0, plot="ADD", col="green")
par(new=F); bm.sim(2, 100, 0, plot=T, col="red")


# We therefore end up with four parameters: theta, the optimal 
# trait value, alpha, the tendency to stay close to the optimum, 
# sigma, stochastic motion parameter, and x0, the starting value 
# for the trait.
  
t <- 0:60
sig2 <- 0.03
nsim <- 100
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-5, 5), type = "l", col="gray")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

# working off the BM set-up above, this simulates a rate declining model (EB-esque)
Y <- NULL
#decline <- seq(0.0001,1,0.016)
#decline <- rev(decline)
#dec1 <- seq(0.8,0.999,0.015); dec1 <- rev(dec1)
#dec2 <- seq(0.1,0.75, 0.005); dec2 <- rev(dec2)
#decline <- append(dec1, dec2)
inc1 <- seq(1.2,3,0.1); inc1 <- rev(inc1)
inc2 <- seq(1,1.2,0.0045); inc2 <- rev(inc2)
incline <- append(inc1,inc2)

for (i in 1:length(X[1,])) {
 rate <- X[,i]*incline[i]
 Y <- cbind(Y, rate)
}
plot(t, Y[1, ], xlab = "time", ylab = "phenotype", ylim = c(-5, 5), type = "l")
apply(Y[2:nsim, ], 1, function(x, t) lines(t, x), t = t)



  