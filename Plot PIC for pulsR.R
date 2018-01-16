#### pulsR can be a bit annoying because nodes with 0 value contrasts (PIC) 
### cause a singularity in the pure jump models (JN, NIG) resulting in 
## erroneously high likelihoods and favored models
# to get around this, we need to adjust our data just enough to have no 0 value contrasts

#######################################################
# Read in all the data you'll use
######################################################
#tree<-read.nexus("BT.Pygopodoidea.tre") #read in desired tree
trees <- read.nexus("PB.Australian.Marsupials.100.trees")
#drip <- c("S_Papuascincus_sp","S_Prasinohaema_virens","S_Sphenomorphus_jobiensis", "S_Sphenomorphus_muelleri","S_Sphenomorphus_solomonis")
#trees<-lapply(trees, drop.tip, tip=drip) #drop tips if necessary


#data<-read.csv("BT.Pygopodoidea.logSVL.csv", row.names = 1, header=TRUE) #read in data file
data       <- read.csv("BT.Australian.Marsupials.MlogBL.csv", row.names = 1, header=F) #read in data file in GEIGER format
data.OUwie <- read.csv("BT.Australian.Marsupials.MlogBL.csv", header=F) #read in data file in OUwie format
data.fitenv<- data.OUwie[,2]; names(data.fitenv) <- data.OUwie[,1] #read in data file in RPANDA format

name.check(trees[[3]], data); name.check(trees[[3]], data.fitenv) #check to make sure the tips match the data labels


### Now perform your Phylogenetically Independent Contrasts
pic.out <- pic(data.fitenv, trees[[1]], var.contrasts = T)
#as.data.frame(pic.out[,1] == 0) # view nodes have 0 value contrasts
as.data.frame(pic.out[,1] < 0.005 & pic.out[,1] > -0.005)
getDescendants(trees[[1]], 277) # get the descendant tips/nodes from those nodes
trees[[1]]$tip.label[c(225,201)] # and check which tips it actually is
# or do it the easier way, plot the contrasts straight onto the tree and look for 0s
plot(trees[[1]], cex=0.3); nodelabels(round(pic.out[,1], 5), adj = c(1, -0.5), frame="n", cex=0.3)
plot(trees[[1]], cex=0.3); nodelabels(frame="c", cex=0.3)


## Below is an attempt at a loop to pull out taxa with 0 value contrasts, 
### it doesn't quite work yet, and I doubt I'll bother
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

load("/Users/Ian/Google.Drive/R.Analyses/pulsR/vertebrate.body_size.raw.rds")
dataz <-readRDS("/Users/Ian/Google.Drive/R.Analyses/pulsR/vertebrate.body_size.raw.rds")
