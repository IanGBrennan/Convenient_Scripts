#library(devtools)
#install_github("hmorlon/PANDA")
library(phytools)
library(RPANDA)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/RPANDA_extras.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/CreateGeoObject_fromSP.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Calculate_AICs.R")


newick <- "((((A:1,B:0.5):2,(C:3,D:2.5):1):6,E:10.25):2,(F:6.5,G:8.25):3):1;"
tree <- read.tree(text=newick)
plot(tree)
    distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/TEST_DATA.csv", header=T)
        distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
            distribution <- distribution[complete.cases(distribution),]


goanna <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/TXXX_Varanus/Trees/Lin&Wiens.Varanus.tre")
trait.data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/TXXX_Varanus/Data/Lin & Wiens - 2016 - Materials/Goanna.Data.csv",
                       row.names=1, header=T)
    trait.data <- trait.data[complete.cases(trait.data),]
          goanna.data <- as.matrix(log(trait.data[,3]))
              rownames(goanna.data) <- rownames(trait.data); 
              colnames(goanna.data) <- "SVL"
              # adjust tree to match data
              all <- goanna$tip.label
              keep <- rownames(goanna.data)
              drops <- setdiff(all, keep)
              goanna <- drop.tip(goanna, drops)
    name.check(goanna, goanna.data)
              

modelBM <- createModel(tree, "BM")
modelOU <- createModel(tree, "OU")
modelPM <- createModel(tree, "PM")
modelMC <- createModel_MC(tree)
    print(modelMC)
geobj <- CreateGeoObject_SP(tree, distribution)
    modelMC_geo <- createModel_MC_geo(tree, geobj)

## Model Attributes:
    # ['name'] - a name,
    # ['paramsNames'] the names of all parameters,
    # ['comment'] a short description,
    # ['period'] the vector of times at which successive branching and death of lineages occur,
    # ['numbersCopy'] vector containing the lineage number which branches or dies at the end of each period,
    # ['numbersPaste'] vector containing the lineage number in which a daughter lineage is placed at the end of each period (zero if the end of the period corresponds to a death),
    # ['initialCondition'] a function of the parameters giving the initial mean and variance of the gaussian process at the root of the tree,
    # ['aAGamma'] the functions corresponding to ai(t), Ai, and Γi(t) that define the evolution of the process on each period, depending on parameters,
    # ['constraints'] a function of the parameters giving the definition range, 
    # ['params0'] a vector of defaut parameter values.


gmodelBM <- createModel(goanna, "BM")
  #gmodelBM['params0'] <- c(1, 1, 1, 0.05)
gmodelOU <- createModel(goanna, "OU")
gmodelPM <- createModel(goanna, "PM")
  gmodelPM['params0'] <- c(5.5, 200, 100, 100, 100, 0.05)
gmodelGM <- createModel(goanna, "GMM")

dataBM <- simulateTipData(modelBM, c(0,0,0,1))
fitTipData(modelBM, dataBM, GLSstyle=T)
fitTipData(modelOU, dataBM, GLSstyle=T)
fitTipData(modelPM, dataBM, GLSstyle=T)
fitTipData(modelMC, dataBM, GLSstyle=T)
fitTipData(modelMC_geo, dataBM, GLSstyle=T)


goannaBM <- simulateTipData(gmodelBM, c(5,0,0,0.05), method=1)
fitTipData(gmodelBM, goannaBM, GLSstyle=T)
fitTipData(gmodelBM, goanna.data, GLSstyle=T)
fitTipData(gmodelOU, goanna.data, GLSstyle=T)
fitTipData(gmodelPM, goanna.data, GLSstyle=T)
getDataLikelihood(gmodelPM, goanna.data, c(-100, -100, -100, 151, -151, 0.06))


## Building the GMM
newick1 <- "(((A:1,B:1):3,(C:3,D:3):1):2,E:6);"
tree_1 <- read.tree(text=newick1)
plot(tree_1)
newick2 <- "((X:1.5,Y:1.5):3,Z:4.5);"
tree_2 <- read.tree(text=newick2)
plot(tree_2)
    endOfPeriodsGMM(tree_1, tree_2)
    modelGMMbis <- createModelCoevolution(tree_1, tree_2, keyword="GMMbis")
    dataGMM <- simulateTipData(modelGMMbis, c(0,0,5,-5, -0.05, 1), method=2)
    getTipDistribution(modelGMMbis, c(0,0,5,-5,-0.5,1))
    fitTipData(modelGMMbis, dataGMM, c(0,0,5,-5,-0.5,1))
    test0 <- fitTipData(modelGMMbis, dataGMM, GLSstyle=T)
      fitTipData(modelGMMbis)
  
modelGMM <- createModelCoevolution(tree_1, tree_2, keyword="GMM")
geo_obj <- CreateCoEvoGeoObject_SP(tree_1, tree_2, distribution)
  modelGMM_pheno <- createGeoModelCoevolution(tree_1, tree_2, geo_obj, keyword="GMM")
modelGMM <- createModelCoevolution(tree_1, tree_2, keyword="GMM")
  modelGMM_geo <- createModelCoevolution(tree_1, tree_2, geo.object=geo_obj, keyword="GMM+geo")
  
  fitTipData(modelGMM, dataGMM, GLSstyle=T)
  fitTipData(modelGMM_geo, dataGMM, GLSstyle=T)
  

dataGMM <- simulateTipData(modelGMM, c(0,0.04,0.01,-0.01, -0.2, 0.04), method=2)

testo <- fitTipData(modelGMM, dataGMM, GLSstyle=T)
    fitTipData(modelGMM, dataGMM, GLSstyle=T, params0=testo$inferredParams)
testa <- fitTipData(modelGMM_pheno, dataGMM, GLSstyle=T)
    fitTipData(modelGMM_pheno, dataGMM, GLSstyle=T, params0=testo$inferredParams)
# we can change the 'PhenotypicModel' class, adjusting how the parameters are estimated
modelGMM.o <- modelGMM; class(modelGMM.o)[1] <- "PhenotypicModel"
    fitGMM.o <- fitTipData(modelGMM.o, dataGMM, GLSstyle=T)
        fitGMM.or <- fitTipData(modelGMM.o, dataGMM, GLSstyle=T, params0=fitGMM.o$inferredParams)
modelGMM.geo <- modelGMM_pheno; class(modelGMM.geo)[1] <- "PhenotypicModel"
    fitGMM.ge <- fitTipData(modelGMM.geo, dataGMM, GLSstyle=T)
        fitGMM.geo <- fitTipData(modelGMM.geo, dataGMM, GLSstyle=T, params0=fitGMM.ge$inferredParams)
multiphy.AIC("fit", tree_1, c("GMM.geo", "GMM.or"))
        
    
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/GMM_TEST_DATA.csv", header=T)
    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]






data(Anolis.data)
#Create a geography.object with a modified edge matrix
#First, specify which region each branch belonged to:
Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
Anolis.map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
CreateGeoObject(Anolis.data$phylo,map=Anolis.map)

#Create a geography.object with a make.simmap object
#First, specify which region each branch belonged to:
require(phytools)
geo<-c(rep("cuba",7),rep("hispaniola",9),"puerto_rico")
names(geo)<-Anolis.data$phylo$tip.label
stochastic.map<-make.simmap(Anolis.data$phylo,geo, model="ER", nsim=1)
map.object <- CreateGeoObject(Anolis.data$phylo,map=stochastic.map)
CreateGeoObject_SP(Anolis.data$phylo, map=distribution)


#######################################################
# Let's run through an example using the anolis tree
#######################################################
data(Anolis.data)
phy <- Anolis.data$phylo
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/ANOLIS_TEST_DATA.csv", header=T)
    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
traits <- Anolis.data$data

# start by building the models (of class 'PhenotypicModel')
modelBM <- createModel(phy, "BM") # Brownian Motion - random, dictated by time and rate
modelOU <- createModel(phy, "OU") # Ornstein-Uhlenbeck - random, dictated by time and rate, but constrained around optimum
modelPM <- createModel(phy, "PM") # Phenotypic Matching - a given lineage's value is influenced by those around it
model_simpPM <- createModel(phy, "PM_OUless")
modelMC <- createModel_MC(phy) # Matching Competition - same as PM model, but with different code
geo_object <- CreateGeoObject_SP(phy, distribution)
  resgeo_object <- resortGeoObject(phy, geo_object)
    modelMC_geo <- createModel_MC_geo(phy, resgeo_object)
    modelPM_geo <- createModel_PM_geo(phy, resgeo_object, "PM")
    modelMC_geo_MINE <- createGeoModel(phy, resgeo_object, "MC")
    modelPM_geo_MINE <- createGeoModel(phy, resgeo_object, "PM")
map.object <- CreateGeoObject(phy,map=stochastic.map)
    modelMC_geo_isl <- createModel_MC_geo(phy, map.object)
    modelPM_geo_isl <- createModel_PM_geo(phy, map.object, "PM")
testMC <- createGeoModel(phy, geo_object, "MC")
testPM <- createGeoModel(phy, geo_object, "PM")
modelDC <- createModel(phy, "ACDC")
    

fitTipData(modelPM, traits, GLSstyle=T)
fitTipData(modelMC_geo, traits, GLSstyle=T)
MC_test <- fitTipData(modelMC_geo_MINE, traits, GLSstyle=T)
PM_test <- fitTipData(modelPM_geo_MINE, traits, GLSstyle=T, 
                      params0=c(mean(traits), 0, mean(traits), 0.01, -0.01, 0.01))
out <- fitTipData(modelPM_geo, traits, GLSstyle=T)
    fitTipData(modelPM_geo, traits, GLSstyle=T, params0=out$inferredParams)
test <- fitTipData(modelPM_geo_MINE, traits, GLSstyle=T)
fitTipData(modelPM_geo_MINE, traits, GLSstyle=T, params0=c(mean(traits), -mean(traits), -0.1, 0.001, -0.1, 0.02))
    test2 <- fitTipData(modelPM_geo_MINE, traits, GLSstyle=T, params0=test$inferredParams)
        test3 <- fitTipData(modelPM_geo_MINE, traits, GLSstyle=T, params0=test$inferredParams)
    

# now fit the model to the trait data
fitBM <- fitTipData(modelBM, traits, GLSstyle=T)
fitOU <- fitTipData(modelOU, traits, GLSstyle=T)
fitPM <- fitTipData(modelPM, traits, GLSstyle=T)
    fitPM2 <- fitTipData(modelPM, traits, GLSstyle=T, params0=fitPM$inferredParams)
    fitPM3 <- fitTipData(modelPM, traits, GLSstyle=T, params0=c(10,10,10,10,10,10))
        fitPM4 <- fitTipData(modelPM, traits, GLSstyle=T, params0=fitPM3$inferredParams)
fitPM_simp <- fitTipData(model_simpPM, traits, GLSstyle=T)
    
fitMC <- fitTipData(modelMC, traits, GLSstyle=T)
fitMC_geo <- fitTipData(modelMC_geo, traits, GLSstyle=T)
    fitMC_geo_test <- fitTipData(testMC, traits, GLSstyle=T)
    getTipDistribution(testMC, fitMC_geo_test$inferredParams)
    simulateTipData(testMC, fitMC_geo_test$inferredParams, method=2)
fitComp <- fit_t_comp(phy, traits, model="MC", geography.object = geo_object, method=2)

fitPM_geo <- fitTipData(modelPM_geo, traits, GLSstyle=T)
    fitPM_geo2 <- fitTipData(modelPM_geo, traits, GLSstyle=T, params0=fitPM_geo$inferredParams)
    fitPM_geo3 <- fitTipData(modelPM_geo, traits, GLSstyle=T, params0=c(10,10,10,10,10,10))
        fitPM_geo_test <- fitTipData(testPM, traits, GLSstyle=T)
            fitPM_geo_test <- fitTipData(testPM, traits, GLSstyle=T, params0=fitPM_geo_test$inferredParams)
fitDC <- fitTipData(modelDC, traits, GLSstyle=T, params0=c(0,0,1,1))
    
res <- as.data.frame(c(fitBM$value, fitOU$value, fitPM2$value, 
                       fitPM_geo2$value, fitMC$value, fitMC_geo$value,
                       fitMC_geo_isl$value, fitPM_geo_isl2$value))
    rownames(res) <- c("BM", "OU", "PM", "PM+geo", "MC", "MC+geo",
                       "MC+islands", "PM+islands"); colnames(res) ="logLikelihood"; res

fitMC_geo_isl <- fitTipData(modelMC_geo_isl, traits, GLSstyle=T)
fitPM_geo_isl <- fitTipData(modelPM_geo_isl, traits, GLSstyle=T)
    fitPM_geo_isl2 <- fitTipData(modelPM_geo_isl, traits, GLSstyle=T, params0=fitPM_geo_isl$inferredParams)
#

#######################################################
# Let's try with the Agamids
#######################################################   
phy <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.tre")

distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Agamids.csv", header=T)
    distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Agamids.logSVL.csv", header=F)

### We don't have distributional data for everything, so we need to drop a few from the tree
drop <- setdiff(phy$tip.label, unique(distribution$Name_in_Tree))
phy <- drop.tip(phy, tip=drop);
### and from the trait data
trait <- trait[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
    traits <- trait[,2]; names(traits) <- trait[,1] #read in data file in RPANDA format

# build a geography object for the tree/distribution data, or load one you've already built
geo_object <- CreateGeoObject_SP(phy, distribution)
    # this probably took a while, so let's save it externally, so we don't have to do it again
    saveRDS(geo_object, file="/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Agamids.geo.object.RDS")
geo_object <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Agamids.geo.object.RDS")
    
# start by building the models (of class 'PhenotypicModel')
modelBM <- createModel(phy, "BM") # Brownian Motion - random, dictated by time and rate
modelOU <- createModel(phy, "OU") # Ornstein-Uhlenbeck - random, dictated by time and rate, but constrained around optimum
modelMC <- createModel_MC(phy) # Matching Competition - same as PM model, but with different code
    modelMC_geo <- createGeoModel(phy, geo_object, "MC")
modelPM <- createModel(phy, "PM") # Phenotypic Matching - a given lineage's value is influenced by those around it
    modelPM_geo <- createGeoModel(phy, geo_object, "PM")
    modelPM_OUless <- createModel(phy, "PM_OUless")

# now fit the model to the trait data
fitBM <- fitTipData(modelBM, traits, GLSstyle=T)
    simulateTipData(modelBM, fitBM$inferredParams, method=2)
    simulateTipData(modelBM, c(2.048, 0.039, 0, 0.039), method=1)
fitOU <- fitTipData(modelOU, traits, GLSstyle=T)
    fitOU <- fitTipData(modelOU, traits, GLSstyle=T, params0=fitOU$inferredParams)
fitMC <- fitTipData(modelMC, traits, GLSstyle=T)
    fitMC_geo <- fitTipData(modelMC_geo, traits, GLSstyle=T)
    fitMC_geo_t_comp <- fit_t_comp(phy, traits, model="MC", geography.object=geo_object) # this is the official RPANDA version
fitPM <- fitTipData(modelPM, traits, GLSstyle=T)
    fitPM <- fitTipData(modelPM, traits, GLSstyle=T, params0=fitPM$inferredParams)
    fitPM_geo <- fitTipData(modelPM_geo, traits, GLSstyle=T)
        fitPM_geo <- fitTipData(modelPM_geo, traits, GLSstyle=T, params0=fitPM_geo$inferredParams)
            simulateTipData(modelPM_geo, fitPM_geo$inferredParams, method=2)
            simPM_geo <- simulateTipData(modelPM_geo, c(1.95, 2.061, 1.51, 0.099, -0.070, 0.048), method=1)
            sim_fitPM_geo <- fitTipData(modelPM_geo, simPM_geo, GLSstyle=T)
            sim_fitPM <- fitTipData(modelPM, simPM_geo, GLSstyle=T)
            sim_fitPM_OUless <- fitTipData(modelPM_OUless, simPM_geo, GLSstyle=T)
            sim_fitBM <- fitTipData(modelBM, simPM_geo, GLSstyle=T)
            sim_fitOU <- fitTipData(modelOU, simPM_geo, GLSstyle=T)
            sim_fitMC <- fitTipData(modelMC, simPM_geo, GLSstyle=T)
            sim_fitMC_geo <- fitTipData(modelMC_geo, simPM_geo, GLSstyle=T)
        multiphy.AIC("sim_fit", phy, model.names)
            
model.names <- c("BM","OU","MC","MC_geo","PM","PM_geo","PM_OUless")

multiphy.AIC("fit", phy, model.names)


    
    