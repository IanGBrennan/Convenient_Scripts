#library(devtools)
#install_github("hmorlon/PANDA")
library(phytools)
library(RPANDA)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/RPANDA_extras.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/CreateGeoObject_fromSP.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Calculate_AICs.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/plot.distmaps.R")


aprasia <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.tre")
  aprasia.trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.data.csv", header=T)
    aprasia.trait <- aprasia.trait[,c("Name_in_Tree", "SVL")]; aprasia.trait$SVL <- log(aprasia.trait$SVL)
      pygo.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Aprasia.csv", header=T)
      pygo.dist <- pygo.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
    ### trim down distributional data to match the tree
    aprasia.dist <- filter(pygo.dist, Name_in_Tree %in% aprasia$tip.label)
    ### and from the trait data
    #aprasia.trait <- filter(aprasia.traits, V1 %in% aprasia$tip.label)
    aprasia.traits <- aprasia.trait[,2]; names(aprasia.traits) <- aprasia.trait[,1] #read in data file in RPANDA format
    
    
lerista <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Lerista.tre")
  skink.traits <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Skinks.logSVL.csv", header=F)
    skink.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Skinks.csv", header=T)
      skink.dist <- skink.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
      ### trim down distributional data to match the tree
      lerista.dist <- filter(skink.dist, Name_in_Tree %in% lerista$tip.label)
      ### and from the trait data
      lerista.trait <- filter(skink.traits, V1 %in% lerista$tip.label)
      lerista.traits <- lerista.trait[,2]; names(lerista.traits) <- lerista.trait[,1] #read in data file in RPANDA format

anilios <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Anilios.tre")  
  anilios.trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Anilios.data.csv", header=T)
    anilios.trait <- anilios.trait[,c("Name_in_Data", "SVL")]; anilios.trait$SVL <- log(anilios.trait$SVL)
      anilios.traits <- anilios.trait[,2]; names(anilios.traits) <- anilios.trait[,1] #read in data file in RPANDA format
  anilios.dist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Anilios.csv", header=T)
    anilios.dist <- anilios.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
      ### trim down distributional data to match the tree
      anilios.dist <- filter(anilios.dist, Name_in_Tree %in% anilios$tip.label) # or anilios$tip.label
        #anilios <- drop.tip(anilios, setdiff(anilios$tip.label, anilios.dist$Name_in_Tree)) # drop tips not present in distribution data
      
      
traits <- append(aprasia.traits, lerista.traits)
co.distribution <- rbind(aprasia.dist, lerista.dist)

joint.geo.object <- CreateCoEvoGeoObject_SP(phy1=aprasia, phy2=lerista, map=co.distribution,
                                            rase.obj1=aprasia.rase, rase.obj2=lerista.rase)
saveRDS(joint.geo.object, file="/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.Lerista.geo.object.RDS")

# joint.geo.object <- CreateCoEvoGeoObject_SP(aprasia, anilios, distribution)
#     joint.geo.object <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.Anilios.geo.object.RDS")
# saveRDS(joint.geo.object, file="/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Aprasia.Anilios.geo.object.RDS")

modelGMM <- createModelCoevolution(aprasia, lerista, keyword="GMM")
fitGMM0 <- fitTipData(modelGMM, traits, GLSstyle=T)
  fitGMM <- fitTipData(modelGMM, traits, GLSstyle=T, params0=fitGMM0$inferredParams)
    simulateTipData(modelGMM, fitGMM$inferredParams, method=2)
modelGMM_geo <- createModelCoevolution(aprasia, lerista, geo.object=joint.geo.object, keyword="GMM+geo")
fitGMM_geo0 <- fitTipData(modelGMM_geo, traits, GLSstyle=T)
  fitGMM_geo <- fitTipData(modelGMM_geo, traits, GLSstyle=T, params0=fitGMM_geo0$inferredParams)
both.phy <- c(aprasia, lerista); class(both.phy) <- "multiPhylo"
  multiphy.AIC("fit", both.phy, models=c("GMM", "GMM0", "GMM_geo", "GMM_geo0"))
simulateTipData(modelGMM, fitGMM_geo$inferredParams, method=2)
    sim.par <- fitGMM$inferredParams; sim.par[3] = -sim.par[4]
        simulateTipData(modelGMM, sim.par, method=2)

modelPM <- createModel(lerista, keyword="PM")
  fitPM0 <- fitTipData(modelPM, lerista.traits, GLSstyle=T)
    fitPM <- fitTipData(modelPM, lerista.traits, GLSstyle=T, params0=fitPM0$inferredParams)
modelMC <- createModel_MC(lerista)
  fitMC <- fitTipData(modelMC, lerista.traits, GLSstyle=T)
lerista.geo <- CreateGeoObject_SP(lerista, lerista.dist)
modelMC_geo <- createModel_MC_geo(lerista, lerista.geo)
  fitMC_geo <- fitTipData(modelMC_geo, lerista.traits, GLSstyle=T)
multiphy.AIC("fit", lerista, models=c("MC", "MC_geo"))

modelPMOUless <- createModel(lerista, keyword="PM_OUless")
  fitPMOUless <- fitTipData(modelPMOUless, lerista.traits, GLSstyle=T)
  simulateTipData(modelPMOUless, c(1.8, 2.8e-8, -2.2e-8, 0.02), method=2)

## I'm trying to think of a way to speed up filling in the geography matrix
### worth coming back to, as I think I can fill in the lower off-diagonals, then
#### just transfer them to the upper off-diagonals
mas.matrix <- matrix(nrow=length(namm), ncol=length(namm))
mas.matrix <- as.data.frame(mas.matrix)

rownames(mas.matrix) <- namm
colnames(mas.matrix) <- namm

for (j in 1:length(namm)) {
  for (k in 1:length(namm)) {
    if (j == k) {
      mas.matrix[j,k] <- 1
    } else {
      mas.matrix[j,k] <- filter(test.mat, (species1==namm[j] & species2==namm[k]) |
                                     (species2==namm[j]) & species1==namm[k])$range_overlap
    }
  }
}

test.mat <- as.data.frame(test.mat)
for (k in 1:length(namm)) {
  inter <- filter(test.mat, species1==namm[k])
  for (p in length(inter[,1])) {
    mas.matrix[inter[p,1], inter[p,2]] <- 
    mas.matrix[which(mas.matrix[inter[p,1],1])]
  }
}

# if you (again) forget that RPANDA's loglik is negative (-)
# you can quickly prove it to yourself here
test <- createModel(aprasia, "BM")
fitTipData(test, data=aprasia.traits, GLSstyle=T)
geiger::fitContinuous(aprasia, aprasia.traits, model="BM")
