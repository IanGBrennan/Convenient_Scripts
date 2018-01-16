library(wesanderson)
library(strap)
library(nlme)
library(mgcv)

tree <- read.nexus('BT.Australian.Marsupials.tre')

max.height <- max(nodeHeights(tree))
tree$root.time <- max.height

timeslice1 <- max.height - 5
timeslice2 <- max.height - 10
timeslices <- c(0,timeslice2, timeslice1)
phy.sliced <- make.era.map(tree, timeslices)
leg <- c("blue3", "red3", "green")
colorz <- wes_palette("Darjeeling", 3, "discrete")
names(leg) <- c(1,2,3)
names(colorz) <- c(1,2,3)
plotSimmap(phy.sliced, colorz, pts=F, ftype="off", lwd=1)

phy.sliced$root.time <- max.height

geoscalePhylo(tree=ladderize(phy.sliced, right=FALSE), label.offset=0, cex.age=0.4, cex.ts=0.4, cex.tip=0.2)
geoscalePhylo(tree=ladderize(tree, right=FALSE), label.offset=0, cex.age=0.4, cex.ts=0.4, cex.tip=0.2)






##########################
## http://www.fromthebottomoftheheap.net/2011/06/12/additive-modelling-and-the-hadcrut3v-global-mean-temperature-series/
##########################

## load custom functions
tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/derivFun.R",
              tmp, method = "wget")
source(tmp)

tmp <- tempfile()
download.file("https://github.com/gavinsimpson/random_code/raw/master/tsDiagGamm.R",
              tmp, method = "wget")
source(tmp)


data
plot(meanCI ~ timing, data=data)
#p1 <- predict(m1$gam, newdata = pdat)
m1 <- gamm(meanCI ~ s(timing, k = 8), data = data)
m2 <- gamm(meanCI ~ s(timing, k = 10), data = data)
mx <- gamm(meanCI ~ s(timing, k = 10), data = sim.data)

m3 <- gamm(meanCI ~ s(timing, k=50), data=data)
m4 <- gamm(meanCI ~ s(timing, k=100), data=data)

plot(mx$gam)
plot(m2$gam, add=T)

curve(m2$gam, add=TRUE)
lines(m2, col="red")
anova(m1$lme, m2$lme, m3$lme, m4$lme)

plot(mx$gam, col="red", ylim=range(c(-0.2, 0.2)), xlim=range(c(0,20))); 
par(new=T); 
plot(m2$gam, col="green", ylim=range(c(-0.2, 0.2)), xlim=range(c(0,20))); 
par(new=T); plot(m3$gam, col="blue")
summary(m1$gam)

acf(resid(m1$lme, type = "normalized"))
pacf(resid(m1$lme, type = "normalized"))

m2 <- gamm(meanCI ~ s(timing, k = 30), data = data, correlation = corARMA(form = ~ timing, p = 2))
