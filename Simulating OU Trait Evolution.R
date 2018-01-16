# ===================How to simulate Ornstein-Uhlenbeck===================


# For an OU model a change in trait X over an increment in time t (dXi(t))= alpha[theta-Xi(t)]dt + sigma*dBi(t). So to simulate this we need 4 parameters, theta (the optimal trait value), alpha (the pull towards the optima), sigma (sigma^2 the stochastic (Brownian) motion parameter) and the starting value for the trait(x0). To keep things simple we are assuming that the species share no common ancestry so dt is ignored (it would also appear in the BM portion where variance of the normal distribution is given as sigma^2*dt).

# So if we start the trait at the optima which is 1

alpha<-0.4
theta<-1
sigma<-0.05
x0<-1

OUalpha0.1_sample1<- x0 + alpha*(theta-x0)+sigma*(rnorm(1,mean=0))

# the next sample would be

OUalpha0.1_sample2<- OU_sample1 + alpha*(theta-OU_sample1)+sigma*(rnorm(1,mean=0))

OUalpha0.1_sample3<- OU_sample2 + alpha*(theta-OU_sample2)+sigma*(rnorm(1,mean=0))

plot(x=0:3, y=c(1,OUalpha0.1_sample1, OUalpha0.1_sample2, OUalpha0.1_sample3), type="l", main="Alpha=0.4, Theta=1, Sigma=0.05, x0=1")


# It will take a while to keep writing these out so we can see the actual pattern, so lets create a simple function to do it for us, it will run the simulation for n time units

ornstein_uhlenbecksim <- function(n,theta,alpha,sigma,x0){# n is the number of time units to simulate for, theta is the optima, alpha is the pull towards the optima, sigma is the rate and X0 is the starting trait value
  dw  <- rnorm(n, 0)
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma*dw[i-1]
  }
  return(x);
}

# we can now use this function to explore the effect that changing various parameters will have on the distribution of the trait value

alpha1sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 1, 0.01, 0)) # replicate just runs it 10 times for us to illustrate the variance between runs


colr<-topo.colors(9)

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )

# play with different parameters.

alpha0.5sigma0.01theta1  <-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.01, 0))
alpha0.1sigma0.01theta1  <-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.01, 0))
alpha0sigma0.01theta1    <-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.01, 0))
alpha1sigma0.1theta1     <-replicate(10, ornstein_uhlenbecksim(100, 1, 1, 0.1, 0))
alpha0.5sigma0.1theta1   <-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.1, 0))
alpha0.1sigma0.1theta1   <-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.1, 0))
alpha0sigma0.1theta1     <-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.1, 0))

par(mfrow=c(4,2))

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha1sigma0.1theta1[,1], type="l", main="alpha1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.1theta1)) lines(alpha1sigma0.1theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.01theta1[,1], type="l", main="alpha0.5 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.01theta1)) lines(alpha0.5sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.1theta1[,1], type="l", main="alpha0.5 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.1theta1)) lines(alpha0.5sigma0.1theta1[,i], type="l", col=colr[i-1] )
plot(alpha0.1sigma0.01theta1[,1], type="l", main="alpha0.1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.01theta1)) lines(alpha0.1sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0.1sigma0.1theta1[,1], type="l", main="alpha0.1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.1theta1)) lines(alpha0.1sigma0.1theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.01theta1[,1], type="l", main="alpha0 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.01theta1)) lines(alpha0sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.1theta1[,1], type="l", main="alpha0 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.1theta1)) lines(alpha0sigma0.1theta1[,i], type="l", col=colr[i-1] )


# Once again the above simulations were assuming the species were independent but we know that they are not, they share evolutionary history and this affects the evolution of their traits. So when simulating along a phylogeny time is also an important component, to be able to account for this you need to include time and dt back into the equations. Fortunately, you don't need to be able to do this yourself as the R package OUwie has a simulation function OUwie.sim. 






library(ape)
library(phytools)
library(OUwie)
tree1		<- rtree(40)
tree1 		<- chronopl(tree1, 1)
regime1		<- rep(c(0,1),each=20)
names(regime1) <- tree1$tip.label
sim1 		<- make.simmap(tree1, regime1, model="ER", nsim=1, pi="estimated")
traits1 <- OUwie.sim(sim1, simmap.tree=TRUE, alpha=c(1.5,1.5), 
sigma.sq=c(0.1,0.1), theta0=2.5, theta=c(1.0,6.0))
traits1 <- data.frame(traits1[,1],regime1,traits1[,2])
OUwie_OUM <- OUwie(sim1, traits1, model="OUM", simmap.tree=TRUE)
par(mfrow=c(1,2))
plotTraitgram(fitted.model=OUwie_OUM, anc.model="ER", reps=5, plot.grey.only = "TRUE")
plotTraitgram(fitted.model=OUwie_OUM, anc.model="ER", reps=5, plot.grey.only = "FALSE")




#####################################
## Create a fxn to simulate OU
## 'ornstein_uhlenbecksim'
#####################################
ornstein_uhlenbecksim <- function(n,theta,alpha,sigma,x0){# n is the number of time units to simulate for, theta is the optima, alpha is the pull towards the optima, sigma is the rate and X0 is the starting trait value
  dw  <- rnorm(n, 0)
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + alpha*(theta-x[i-1]) + sigma*dw[i-1]
  }
  return(x);
}
#####################################
# ornstein_uhlenbecksim(n, theta, alpha, sigma, x0)

dub <- replicate(10, ornstein_uhlenbecksim(1000, 1, c(0.1,0.5), 2, 0))
plot(dub[,1], type="l", main="dub")
for(i in 2:ncol(dub)) lines(dub[,i], type="l",col=colr[i-1] )


alpha1sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(1000, 1, 1, 0.01, 0)) # replicate just runs it 10 times for us to illustrate the variance between runs
colr<-topo.colors(9)

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )
# play with different parameters.

alpha0.5sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.01, 0))
alpha0.1sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.01, 0))
alpha0sigma0.01theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.01, 0))
alpha1sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 1, 0.1, 0))
alpha0.5sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.5, 0.1, 0))
alpha0.1sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0.1, 0.1, 0))
alpha0sigma0.1theta1<-replicate(10, ornstein_uhlenbecksim(100, 1, 0, 0.1, 0))



par(mfrow=c(4,2))

plot(alpha1sigma0.01theta1[,1], type="l", main="alpha1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.01theta1)) lines(alpha1sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha1sigma0.1theta1[,1], type="l", main="alpha1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha1sigma0.1theta1)) lines(alpha1sigma0.1theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.01theta1[,1], type="l", main="alpha0.5 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.01theta1)) lines(alpha0.5sigma0.01theta1[,i], type="l",col=colr[i-1] )
plot(alpha0.5sigma0.1theta1[,1], type="l", main="alpha0.5 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.5sigma0.1theta1)) lines(alpha0.5sigma0.1theta1[,i], type="l", col=colr[i-1] )
plot(alpha0.1sigma0.01theta1[,1], type="l", main="alpha0.1 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.01theta1)) lines(alpha0.1sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0.1sigma0.1theta1[,1], type="l", main="alpha0.1 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0.1sigma0.1theta1)) lines(alpha0.1sigma0.1theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.01theta1[,1], type="l", main="alpha0 sigma=0.01", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.01theta1)) lines(alpha0sigma0.01theta1[,i], type="l", col=colr[i-1] )
plot(alpha0sigma0.1theta1[,1], type="l", main="alpha0 sigma=0.1", ylim=c(-2, 2))
for(i in 2:ncol(alpha0sigma0.1theta1)) lines(alpha0sigma0.1theta1[,i], type="l", col=colr[i-1] )

