library(geiger)

###########################################
# START HERE: Fit various types of models to your tree/traits
###########################################
tree<- read.tree("pygopodoidea.tre")
data <- read.csv("body.biome.csv", row.names=1, header=TRUE)

# You have to fit both a MODEL and a TRANSFORMATION to get at the questions!
##########################################

# Fit Size
##########################################
#trait<-data["logSVL"] #don't do this, fitting a discrete model to cont data will freeze R
none.size<-fitDiscrete(tree, trait, transform="none", niter=1000)
eb.size<-fitDiscrete(tree, trait, transform="EB", niter=100)
delta.size<-fitDiscrete(tree, trait, transform="delta", niter=1000)
kappa.size<-fitDiscrete(tree, trait, transform="kappa", niter=1000)
results<-list(none.size$opt$aic, eb.size$opt$aic, delta.size$opt$aic, kappa.size$opt$aic)
#a discrete fitting of body size fits a null (BM) model
##########################################
# Fit Ecology
##########################################
trait<-data["ecology"]
none.ecology<-fitDiscrete(tree, trait, transform="none", niter=1000)
none.er.ecology<-fitDiscrete(tree, trait, transform="none", model="ER", niter=1000)
none.ard.ecology<-fitDiscrete(tree, trait, transform="none", model="ARD", niter=1000)

eb.er.ecology<-fitDiscrete(tree, trait, transform="EB", model="ER", niter=1000)
eb.ard.ecology<-fitDiscrete(tree, trait, transform="EB", model="ARD", niter=1000)
eb.er.ecology$opt$aicc
eb.ard.ecology$opt$aicc


eb.ecology<-fitDiscrete(tree, trait, transform="EB", niter=1000)
delta.ecology<-fitDiscrete(tree, trait, transform="delta", niter=1000)
kappa.ecology<-fitDiscrete(tree, trait, transform="kappa", niter=1000)
results<-list(none.ecology$opt$aic, eb.ecology$opt$aic, delta.ecology$opt$aic, kappa.ecology$opt$aic)
#an eb model best ecological evolution

# Fit Toe Morphology
##########################################
trait<-data["toes"]
none.toes<-fitDiscrete(tree, trait, transform="none", niter=1000)
eb.toes<-fitDiscrete(tree, trait, transform="EB", niter=1000)
delta.toes<-fitDiscrete(tree, trait, transform="delta", niter=1000)
kappa.toes<-fitDiscrete(tree, trait, transform="kappa", niter=1000)
results<-list(none.toes$opt$aic, eb.toes$opt$aic, delta.toes$opt$aic, kappa.toes$opt$aic)
#an eb model best fits evolution of toe morphology

# Fit Biome Distribution
##########################################
trait<-data["biome"]
none.biome<-fitDiscrete(tree, trait, transform="none", niter=1000)
eb.biome<-fitDiscrete(tree, trait, transform="EB", niter=1000)
delta.biome<-fitDiscrete(tree, trait, transform="delta", niter=1000)
kappa.biome<-fitDiscrete(tree, trait, transform="kappa", niter=1000)
results<-list(none.biome$opt$aic, eb.biome$opt$aic, delta.biome$opt$aic, kappa.biome$opt$aic)
#a kappa model marginally (non-significantly) best fits evolution of biome distributions
##########################################

# Fit Multiple Ecological Characters (here: strata and diel)
##########################################
trait<-data["fitDis"] #set your data
none.er.all<-fitDiscrete(tree, trait, transform="none", model="ER", niter=1000, ncores=NULL) #model with no transformation and equal rates
#none.ard.all<-fitDiscrete(tree, trait, transform="none", model="ARD", niter=100, ncores=NULL) #model with no transformation and all rates different

eb.er.all<-fitDiscrete(tree, trait, transform="EB", model="ER", niter=1000, ncores=NULL) #early burst transformation, equal rates
eb.ard.all<-fitDiscrete(tree, trait, transform="EB", model="ARD", niter=1000, ncores=NULL) # early burst transformation, all rates different
#eb.sym.all<-fitDiscrete(tree, trait, transform="EB", model="SYM", niter=1000, ncores=NULL) #early burst transformation, symmetrical rates
#eb.mer.all<-fitDiscrete(tree, trait, transform="EB", model="meristic", ncores=NULL) #consider creating a matrix for the transitions

delta.er.all<-fitDiscrete(tree, trait, transform="delta", model="ER", niter=1000, ncores=NULL) #delta transform, equal rates
delta.ard.all<-fitDiscrete(tree, trait, transform="delta", model="ARD", niter=1000, ncores=NULL) #delta transform, all rates different

kappa.er.all<-fitDiscrete(tree, trait, transform="kappa", model="ER", niter=1000, ncores=NULL) #kappa transform, equal rates
kappa.ard.all<-fitDiscrete(tree, trait, transform="kappa", model="ARD", niter=1000, ncores=NULL) #kappa transform, all rates different

lambda.er.all<-fitDiscrete(tree, trait, transform="lambda", model="ER", niter=1000, ncores=NULL) #lambda transform, equal rates
lambda.er.all<-fitDiscrete(tree, trait, transform="lambda", model="ARD", niter=1000, ncores=NULL) #lambda transform, all rates different

# Create a Table to hold the results across models
res<-as.matrix(eb.er.all$opt) #start with first
res<-cbind(res, delta.er.all$opt) #add another column of results
res<-cbind(res, kappa.er.all$opt) #etc
res<-cbind(res, lambda.er.all$opt) #etc
#res<-as.data.frame(res)
#colnames(res) <- c("early burst","delta", "kappa", "lambda")
write.csv(res, file="Pygo.fitDiscrete results.csv") #write results to file

none.er.all$opt #this is left out, because it's a different length from the other results

################################################################################
######## Using the Slater Code to DTT discrete characters as continuous ########
################################################################################
## this script will repeat the analysis of disparity through time conducted in Slater et al. 2010 and reproduce figure 4
#DTT on SVL of Pygopodoidea

trait<-setNames(data[,7], rownames(data)) #choose your data file (here: svl) and column ([,5]), apply rownames
source("dtt.full with confidence limits and p value.R"); ## now source the function from your WD
# and run the analysis. i recommend using ~10,000 sims to get a stable P value
##if you've used the regular dtt.full, you'll notice that this function doesn't immediately plot up the empirical dtt curve. don't worry - it's just that the way I get the nice shaded 95% region requires a blank plot at first with the empirical and median curves overlain on the shaded region.
dttOut<-dttFullCIs(tree, trait, nsims=1000);
#dttOut<-dttFullCIs(pygopodoidea, data, nsims=1000, mdi.range=c(0,0.5)); #designate a time slice to measure the MDI, here first half of the tree
dttOut.burst<-dttFullCIs(tree, trait, nsims=1000, mdi.range=c(0.25,0.5)); #designate a time slice to measure the MDI, here second half of the tree
dttOut.modern<-dttFullCIs(tree, trait, nsims=1000, mdi.range=c(0,0.25)); #designate a time slice to measure the MDI, here third quarter

dttOut$MDI; dttOut$Pvalue
dttOut.burst$MDI; dttOut.burst$Pvalue
dttOut.modern$MDI; dttOut.modern$Pvalue
