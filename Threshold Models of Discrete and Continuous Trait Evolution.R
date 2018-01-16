library(phytools)

tree <- pbtree(n=50, scale=1)
r <- 0.75 # simulate using high correlation
V <- matrix(c(1,r,r,1),2,2) # build your matrix
X <- sim.corrs(tree, V)
sample <- 500
ngen <- 2e+05
burnin <- 0.2*ngen
AA <- threshBayes(tree, X, types=c("cont", "cont"), ngen=ngen,
                  control=list(sample=sample))

x <- fastBM(tree, sig2=1, a=0.5, internal=T)
th <- sapply(x, threshState, thresholds=setNames(c(0,0.5,2,Inf), letters[1:4]))

data <- read.csv("/Users/Ian/Desktop/Pygo.test.csv", row.names=1)
data.thresh <- data[,2:3]
eco <- data$Ecology
names(eco) <- rownames(data)
hl <- data$HL.Trunk
names(hl) <- rownames(data)
eco.hl <- NULL
eco.hl <- cbind(eco.hl, eco)
eco.hl <- cbind(eco.hl, hl)
mcmc <- ancThresh(new.tree, eco.hl, ngen=100000,
                  control=list(sample=sample,plot=F, print=T))

plotTree(new.tree,ftype="off")
#tiplabels(pie=to.matrix(eco.hl[,1], ),piecol=palette()[1:3],cex=0.3)
nodelabels(pie=mcmc$ace,piecol=palette()[1:3],cex=2)

colMeans(mcmc$par[(0.2*ngen/sample):(ngen/sample)+1,])
plot(density(mcmc$par[(0.2*ngen/sample):(ngen/sample)+1,
                      letters[2]],bw=0.5),xlab="threshold a->b",main="")
lines(c(0.5,0.5),c(0,1000),lty="dashed")

n <- 50
ngen <- 100000
sample <- 500
tree <- pbtree(n=n, scale=1)
x <- fastBM(tree, sig2=1, a=0.5, internal=T)

th<-sapply(hl,threshState,thresholds=setNames(c(0,0.05,0.09,1), c("no", "under", "mid", "over")))
mcmc <- ancThresh(new.tree, th[1:n], ngen=ngen, sequence=c("no", "under", "mid", "over"), 
                  control=list(sample=sample, plot=F))

plotTree(tree,ftype="off")
tiplabels(pie=to.matrix(th[1:length(new.tree$tip)],c("no", "under", "mid", "over")),piecol=palette()[1:3],cex=0.8)
nodelabels(pie=mcmc$ace,piecol=palette()[1:4],cex=1)

colMeans(mcmc$par[(0.2*ngen/sample):(ngen/sample)+1,c("no", "under", "mid", "over")])

plot(density(mcmc$par[(0.2*ngen/sample):(ngen/sample)+1, "under"],bw=0.5),xlab="threshold under->mid",main="")
lines(c(4,0.5),c(0,1000),lty="dashed")

plot(density(mcmc$par[(0.2*ngen/sample):(ngen/sample)+1,letters[3]],bw=0.5),xlab="threshold b->c",main="")
lines(c(12,0.1),c(0,1000),lty="dashed")







plotThresh(new.tree, eco.hl, mcmc, piecol=colors, cex=0.2)

rm(colnames(data.thresh))
mcmcT <- threshBayes(new.tree, hl, types=c("continuous", "continuous"), 
                    ngen=100000, control=list(sample=sample))

V.mle<-phyl.vcv(eco.hl,vcv(new.tree),lambda=1)$R
V.mle[1,2]/sqrt(V.mle[1,1]*V.mle[2,2])
mean(mcmc$par[(burnin/sample+1):nrow(mcmc$par),"r"])

tree <- read.nexus("BT.Pygopodoidea.tre")

test <- rep("blue", length(data))



lio.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Liolaemus/Liolaemidae.data.csv")
lio.tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Liolaemus/Liolaemidae_MCC_June.tre")
ldata <- lio.data[,c("Parity", "Max_altitude")]
rownames(ldata) <- lio.data[,"Name_in_tree"]
ldata$Parity <- as.character(ldata$Parity)
ldata$Parity[ldata$Parity == "O"] <- "0"
ldata$Parity[ldata$Parity == "V"] <- "1"
ldata$Parity <- as.numeric(ldata$Parity)
ldata$Max_altitude <- as.numeric(ldata$Max_altitude)
ldata <- as.matrix(ldata)
## set some parameters for the MCMC
sample<-10000
ngen<-100000 ## in 'real' studies this should be larger
burnin<-0.2*ngen

mcmc_liolaemidae<-threshBayes(lio.tree,ldata,types=c("disc","cont"),
                     ngen=ngen,control=list(sample=sample))
mean.r <- mean(mcmc_liolaemidae$par[(burnin/sample+1):nrow(mcmc_liolaemidae$par),"r"])
plot(mcmc_liolaemidae$par[,"gen"],mcmc_liolaemidae$par[,"logL"],type="l",
     xlab="generation",ylab="logL")
histo <- plot(density(mcmc_liolaemidae$par[(burnin/sample+1):nrow(mcmc_liolaemidae$par),
                         "r"],bw=0.1),xlab="r",main="posterior density for r")
lines(c(mean.r,mean.r),c(0,3), lty="dashed")

bmPlot(lio.tree, type="threshold", thresholds=c(0,1,2))



ngen=1000000
sample=10000
x<-setNames(as.character(lio.data[,"Viviparous.Oviparous"]),rownames(ldata))
mcmc.o<-ancThresh(lio.tree,x,ngen=ngen,sequence=c("O","V"),
                control=list(sample=sample,plot=FALSE))
mcmc.v<-ancThresh(lio.tree,x,ngen=ngen,sequence=c("V","O"),
                control=list(sample=sample,plot=FALSE))

plotTree(lio.tree,ftype="off")
tiplabels(pie=to.matrix(x,c("O","V")), # make sure the order of the traits matches the ancThresh input order!
          piecol=palette()[1:2],cex=0.2)
nodelabels(pie=mcmc.o$ace,piecol=palette()[1:2],cex=0.4)
colMeans(mcmc.o$par)
colMeans(mcmc.o$par[(0.2*ngen/sample):(ngen/sample)+1,c("O","V")])




library(caper)
library(phytools)

all.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Liolaemus/Liolaemus.data.csv", header=T)
lio.tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Liolaemus/Liolaemidae_MCC_June.tre")
lio.data <- lio.data[,c("Viviparous.Oviparous", "Max_altitude")]

lio <- comparative.data(phy=lio.tree, data=all.data, names.col=Name_in_tree, vcv=F, na.omit=F)
model.pgls <- pgls(Max_altitude ~ Viviparous.Oviparous, data=lio, lambda="ML")
summary(model.pgls)
test <- plot(Max_altitude ~ Viviparous.Oviparous, data=all.data)
abline(model.pgls)

violin <- (ggplot(all.data, aes(Viviparous.Oviparous, fill=Viviparous.Oviparous, Max_altitude))
      + geom_violin(scale="count", draw_quantiles=c(0.25,0.5,0.75)) 
      + geom_jitter(height = 0, width = 0.05) + theme_classic())
boxes <- (ggplot(all.data, aes(Viviparous.Oviparous, fill=Viviparous.Oviparous, Max_altitude))
           + geom_box(scale="count") 
           + geom_jitter(height = 0, width = 0.05) + theme_classic())


multiplot(p,test)






sim.lio <- make.simmap(lio.tree, ldata[,1], model="SYM")
phenogram(sim.lio, ldata[,2], fsize=0.6,spread.costs=c(1,0), colors=setNames(c("blue","red"),c(0,1)))


(ggplot(test, aes(x=Max_altitude, y=Max_altitude, fill=Viviparous.Oviparous, color=Viviparous.Oviparous))
  + geom_point())




