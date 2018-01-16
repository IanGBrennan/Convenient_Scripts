library(RPANDA)
library(geiger)
library(mvtnorm)
library(ape)
library(ggplot2)

#####
# This needs to get included in the geiger/ouwie loop
#####


data("Cetacea")
data(InfTemp)

trait <- sim_t_env(Cetacea, param=c(0.1,-0.2), env_data=InfTemp, model="EnvExp", 
                   root.value=0, step=0.001, plot=TRUE)
plot(trait)
res.exp <- fit_t_env(Cetacea, trait, env_data=InfTemp, df=10, scale=T)
lines(res.exp, col="red")
res.exp$LH; res.exp$aic; res.exp$aicc

tree <- read.nexus("BT.Australian.Marsupials.tre")
data <- read.csv("BT.Australian.Marsupials.mlogBL.csv", header=F) #, row.names=1)
trait <- data[,2]
names(trait) <- data[,1]
trait

res.exp.50 <- fit_t_env(tree, trait, env_data=InfTemp, df=50, scale=T, plot=T)
res.exp.10 <- fit_t_env(tree, trait, env_data=InfTemp, df=10, scale=T, plot=T)

res.exp.50$LH; res.exp.50$aic; res.exp.50$aicc
res.exp.10$LH; res.exp.10$aic; res.exp.10$aicc

plot(InfTemp, col="black")
plot(res.exp.10, col="red")
lines(res.exp.50, col="blue")

head(InfTemp)
(ggplot(data=InfTemp, mapping=aes(x=Age, y=Temperature)) 
  + geom_point()
  + geom_smooth()
  + scale_x_reverse())
par()

#the fxn 'fitContinuous_paleo' allows timed shifts and release of trait evolution
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Slater.2013.Fit Release.Radiate Model.R"); ## now source the function from your WD
gfit <- fitContinuous_paleo(tree, data, model="TRC", shift.time=8)
gfit$Trait1$lnl; gfit$Trait1$aic; gfit$Trait1$aicc


############################

tree <- read.nexus("BT.Pygopodoidea.tre")
emp.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BioGeoBEARS/Data.EMPIRICAL/BGB.Australian.Marsupials.UPDATED.EMPIRICAL.VJ.test.txt", header=T, sep="\t") #, row.names=1)
sim.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BioGeoBEARS/Data.SIMULATED/BGB.Australian.Marsupials.UPDATED.SIMULATIONS.VJ.test.txt", header=T, sep="\t")

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Slater.2013.Fit Release.Radiate Model.R"); ## now source the function from your WD

trait.data.geiger  <- read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=F) #read in data file in GEIGER format
trait.data <- read.csv("BT.Pygopodoidea.logSVL.csv", header=F) #, row.names=1)
trait <- trait.data[,2]
names(trait) <- trait.data[,1]

res.vj <- fit_t_env(tree, trait, env_data=emp.data, df=8, scale=T)
plot(res.vj, col="red")
res.vj$LH; res.vj$aic; res.vj$aicc

res.env <- fit_t_env(tree, trait, env_data=InfTemp, df=10, scale=T)
plot(res.env, col="blue")
res.env$LH; res.env$aic; res.env$aicc

trc <- fitContinuous_paleo(tree, trait.data.geiger, model="TRC", shift.time=11)
trc$Trait1$lnl; trc$Trait1$aic; trc$Trait1$aicc

anova()



fifteen <- smooth.spline(InfTemp, df=15)
all <- smooth.spline(InfTemp)
plot(fifteen)
plot(all)
