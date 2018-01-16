# Functions for simulating various continuous models of trait evolution without phylogenies 

## Written by Samantha Price - July 2014, for the Paleobiological and phylogenetic approaches to macroevoution NESCent Academy

# == Brownian motion

bmsimple<-function(time, stdev){ # time the length of the simulation and stdev the standard deviation of the normal distribution (in this simple case it can be thought of as sigma^2 as the time is not an issue), assumes the mean is 0 .
  displace<-rnorm(time, mean=0, sd=stdev) 
  traj<-cumsum(displace) # calculates the trajectory over time		
  return(traj) 
}


# == Brownian with a trend in the mean trait value

## a function that requires mu which is the trend parameter showing how the mean of the normal distribution changes over time, time the length of the simulation and stdev the standard deviation of the normal distribution, it assumes the starting ancestral trait value is 0


bmtrend<-function(time, mu, stdev){	displace<-vector(,length=time) # set up an empty vector of the same length as time
for(i in 1:time){
  displace[i]<-rnorm(1, mean=mu, sd=stdev) # sample from rnorm with a mean of mu multiplied by the iteration (i.e. time)
}
traj<-cumsum(displace) # calculates the trajectory over time		
return(traj) 
}


# == Ornstein-Uhlenbeck 

## time is the length of the simulation, theta is the optima, alpha is the pull towards the optima, sigma is the rate and X0 is the starting trait value, you could combine the optima and X0 but I find it helpful to separate them as you can also investigate how changes in parameter values change how quickly the optima is reached

ornstein_uhlenbecksim <- function(time,theta,alpha,sigma2,x0){
  dw  <- rnorm(time, 0)
  x <- c(x0) # the first trait value is the ancestral value given by x0
  for (i in 2:(time+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma2*dw[i-1] # this then calculates the trait value for each subsequent generation 
  }
  return(x);
}


# == Early or Late Burst

### a function that requires r which determines the exponential change in the rate over time, time the length of the simulation and the initial rate of Brownian motion, it assumes the starting ancestral trait value is 0

eblb<-function(time, r, sigma0){ 
  displace<-vector(,length=time) # set up an empty vector of the same length as time
  for(i in 1:time){
    e<-exp(1)
    displace[i]<-rnorm(1, mean=0, sd=sigma0*e^(-r*i)) # sample from normal distribution with mean of zero and standard deviation determined by sigma0*e^(-r*i)			
  }
  eblbtraj<-cumsum(displace) # calculates the trajectory over time
  return(eblbtraj) 
}