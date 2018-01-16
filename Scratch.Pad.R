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
plotTraitgram(fitted.model=OUwie_OUM, anc.model="ER", reps=2, plot.grey.only = "FALSE")



tree <- read.nexus('BT.Pygopodoidea.tre')
data.OUwie <- read.csv("BT.Pygopodoidea.logSVL.csv", header=F) #read in data file in OUwie format

max.height <- max(nodeHeights(tree))
timeslices <- max.height - 10
timeslices <- c(0,55)
phy.sliced <- make.era.map(tree, timeslices)
leg <- c("blue3", "red3")
names(leg) <- c(1,2)
plotSimmap(phy.sliced, leg, pts=F, ftype="off", lwd=1)

alpha=c(1.0,0.5)
sigma.sq=c(0.45,0.9)
theta0=1.0
theta=c(1.0,2.0)
sim.time <- OUwie.sim(phy.sliced,simmap.tree=T,scaleHeight=FALSE,
                    alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)

regime <- 1
regime1		<- rep(1,100)
regime0 <- rep(0, 89)
regime <- append(regime1, regime0)
names(regime) <- tree$tip.label
sim 		<- make.simmap(tree, regime, model="ER", nsim=1, pi="estimated")

traits <- OUwie.sim(sim, simmap.tree=TRUE, alpha=c(1.5,1.5), 
                     sigma.sq=c(0.1,0.1), theta0=2.5, theta=c(1.0,6.0))
traits <- data.frame(traits[,1],regime,traits[,2])

emp.OUM <- OUwie(sim, traits, model="OUM", simmap.tree=T)

plotTraitgram(fitted.model=emp.OUM, anc.model="ER", reps=2, plot.grey.only = "FALSE")







tree <- read.nexus("BT.Varanidae.tre")
sig2 <- 0.01
x <- fastBM(tree, sig2=sig2, internal=T)
data <- read("BT.Varanidae.logTL.csv", header=F, row.names=1)
name.check(tree, data)
dataz<-setNames(data[,1], rownames(data)) #choose your data file (here: svl) and column ([,5]), apply rownames
phenogram(tree, dataz, spread.labels=T, spread.cost=c(1,0))
contMap(tree, dataz, lwd=7)
xx <- bmPlot(tree, sig2 = 5/500, type="threshold", thresholds = c(-2, 2))
bmPlot(tree, sig2 = c(0.05, 0.0000005))

t <- 0:100
nsim <- 100
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 1)
X <- cbind(rep(0, nsim), t(apply(X, 1, cumsum)))
plot(t, X[1, ], xlab = "time", ylab = "phenotype", ylim = c(-2, 2), type = "l")
apply(X[2:nsim, ], 1, function(x, t) lines(t, x), t = t)

BM = fastBM(tree, a=2, sig2=0.75, alpha=1e-10, theta=2.9, internal = TRUE,)
phenogram(tree, BM, spread.labels = TRUE, spread.cost=c(1,0))

highOU = fastBM(tree, a=2, sig2=0.75, alpha=10, theta=2.9, internal = TRUE,)
phenogram(tree, highOU, spread.labels = TRUE, spread.cost=c(1,0))

# Alpha acts as a 'rubber band', pulling the trait to the optimum:
x = fastBM(tree, a=0, sig2=1.0, alpha=2.0, theta=4.0, internal = TRUE)
phenogram(tree, x, spread.labels = TRUE)

# Alpha acts as a 'rubber band', pulling the trait to the optimum:
x = fastBM(tree, a=0, sig2=c(2, 0.000000005), alpha=c(0.0001, 3), theta=4.0, internal = TRUE)
phenogram(tree, x, spread.labels = TRUE)












