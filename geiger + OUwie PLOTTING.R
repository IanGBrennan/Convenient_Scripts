library(Rmisc)
library(phytools)
library(geiger)
library(mvtnorm)
library(ggplot2)
library(OUwie)
library(plyr)
library(dplyr)
library(wesanderson)
library(Rmisc)
library(phylobase)
library(RPANDA)
library(scales)


total.results <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Comparison TOTAL/ALLFINAL.Australian.Marsupials.Model.Comparison.TOTAL.csv", header=T)
#total.results <- subset(total.results, total.results$tree.no < 26)
total.results[,1] <- NULL

## If you want to drop some models before you determine AICWt [check with 'levels(total.results$model)']
#total.results <- total.results[!(total.results$model=="kappa"),]
#total.results <- total.results[!(total.results$model=="BM"),]
#total.results <- total.results[!(total.results$model=="lambda"),]
#total.results <- total.results[!(total.results$model=="BMS"),] # check it isn't "BMS " or "bms " or something similar
#total.results <- total.results[!(total.results$model=="OU"),]
total.results <- total.results[!(total.results$model=="BMJN"),]
total.results <- total.results[!(total.results$model=="BMNIG"),]


## Use AIC weights to determine best fitting model and model contributions
total.results$tree.no <- as.numeric(total.results$tree.no) # set the tree numbers as actual numbers
sorted.first <- subset(total.results, total.results$tree.no < 10) # subset 1-9
sorted.first <- arrange(sorted.first, tree.no) # sort this subset
sorted.second <- subset(total.results, total.results$tree.no > 9) # subset 10-100
sorted.second <- arrange(sorted.second, tree.no) # sort this subset, this is done to avoid counting errors
sorted.results <- rbind.data.frame(sorted.first, sorted.second) # now bind them together
weight.delta <- NULL
for (j in 1:25){
  targetdata <- subset(sorted.results, sorted.results$tree.no == j)
  res <- aicw(as.numeric(targetdata$AICc))
  weight.delta <- rbind.data.frame(weight.delta, res)
}
weight.delta <- within(weight.delta, rm(fit))
sorted.results <- cbind.data.frame(sorted.results, weight.delta)

outz <- summarySE(total.results, measurevar="w", groupvars="model")

mars <- (ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
  + geom_text(data=outz, aes(x = model, y = ((w+se) + .02), # add percents to the top of each bar!
                             label = percent(w)), colour="black", size = 2)
)
multiplot(agam, bird, 
          pygo, skink, 
          mars, cols=2)

multiplot(mars, mars, 
          mars, mars, 
          mars, cols=2)
agam

# sort the models into the order you'd like
outz$model <- factor(outz$model, 
                     levels=c("JN", "NIG", "BMJN", "BMNIG",
                              "BM", "BMS", "delta", "kappa", "lambda", "TS", "DensDep",
                              "EB", "OU", "OUMA", "OUMV", "OUMVA", "OUS",
                              "SRC", "TRC")) # re-order the models in the plot
# if you just want to do a single plot
myplot <- (ggplot(outz)
           + geom_bar(aes(1, y=w, fill=outz$model), stat="identity")
           + scale_fill_manual( values=wes_palette("Zissou", 19, "continuous")))



#######################################################
# Read in all the data you'll use
######################################################
#tree<-read.nexus("BT.Pygopodoidea.tre") #read in desired tree
trees <- read.nexus("PB.Skinks.100.trees")
#drip <- c("S_Papuascincus_sp","S_Prasinohaema_virens","S_Sphenomorphus_jobiensis", "S_Sphenomorphus_muelleri","S_Sphenomorphus_solomonis")
#trees<-lapply(trees, drop.tip, tip=drip) #drop tips if necessary


#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data       <- read.csv("BT.Skinks.logSVL.csv", row.names = 1, header=F) #read in data file in GEIGER format
data.OUwie <- read.csv("BT.Skinks.logSVL.csv", header=F) #read in data file in OUwie format
data.fitenv<- data.OUwie[,2]; names(data.fitenv) <- data.OUwie[,1] #read in data file in RPANDA format

name.check(trees[[3]], data); name.check(trees[[3]], data.fitenv) #check to make sure the tips match the data labels

taxa <- NULL
for (j in 1:10) {
  pic.out <- pic(data.fitenv, trees[[j]], var.contrasts=T)
  zeroes <- subset(pic.out, pic.out[,1] == 0)
  nodez <- rownames(zeroes)
  for (k in 1:length(nodez)) {
    interim <- getDescendants(trees[[j]], nodez[k])
    taxon <- NULL
    for (p in 1:length(interim)) {
      namez <- trees[[j]]$tip.label[interim[p]]
      taxon[p] <- namez
    }
    taxa[[j]] <- append(taxa, taxon)
  }
}
pic.out <- pic(data.fitenv, trees[[10]], var.contrasts = T)
as.data.frame(pic.out[,1] == 0)
getDescendants(trees[[1]], 239)
trees[[1]]$tip.label[c(4,10)]
plot (trees[[10]]); nodelabels(round(pic.out[,1], 5), adj = c(1, -0.5), frame="n", cex=0.3)



test <- subset(pic.out, pic.out[,1] == 0)
nodez <- rownames(test)
getDescendants(trees[[4]], nodez[1])
trees[[4]]$tip.label[53]

sort(pic.out[,1])
plot(trees[[4]]); nodelabels(frame="c", cex=0.3)
