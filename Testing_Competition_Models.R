#library(devtools)
#install_github("hmorlon/PANDA")
library(phytools)
library(RPANDA)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/RPANDA_extras.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/CreateGeoObject_fromSP.R")
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Calculate_AICs.R")

tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Varanus.LW.tre")
  varanus.traits <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Modelling_Competition/Varanus.Data.csv", header=T)
    traits <- filter(varanus.traits, Continent == "Australia"); rownames(traits) <- traits$Name_in_tree
      traits <- traits[,c(4:length(traits))]
        traits <- traits[complete.cases(traits),]
          #trait <- traits[,2]; names(trait) <- trait[,1] #read in data file in RPANDA format
        
      drop <- setdiff(tree$tip.label, traits$Name_in_tree)
        tree <- drop.tip(tree, tip=drop)
    
varanus.dist <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Varanus.csv", header=T)
  varanus.dist <- varanus.dist[,c("Name_in_Tree", "Latitude", "Longitude")]
      varanus.dist <- filter(varanus.dist, Name_in_Tree %in% tree$tip.label)
      
      
against <- datum$Body_Mass; names(against) <- rownames(datum) # phyl.resid needs names to determine proper order
for (b in 1:length(trees)) {
  resid.data <- NULL # temporary object
  #resids <- phyl.resid(trees[[b]], against, datum[,c("Brain.size", "M.avgWT")]) # regress brain size and weight against body length
  resids <- phyl.resid(trees[[b]], against, datum[,2:5]) # regress brain size and weight against body length
  resid.data <- cbind(resid.data, resids$resid); resid.data <- cbind(resid.data, against); 
  #colnames(resid.data) <- c("BrainSize", "Weight", "BodyLength")
  colnames(resid.data) <- c(colnames(datum)[2:5], colnames(datum)[1])
  total.data[[b]] <- resid.data
}

traits <- log(traits)
against <- traits$SVL; names(against) <- rownames(traits) # phyl.resid needs names to determine proper order
trait.data <- traits[,c(2:length(traits))]
residual.data <- phyl.resid(tree, against, trait.data)$resid
    residual.data <- residual.data[order(rownames(residual.data)),]
        SVL <- against; SVL <- as.data.frame(SVL); residual.data <- cbind(SVL, residual.data)

compare.traits <- function(data, trait1, trait2){
  xaxis <- colnames(data[trait1])
  yaxis <- colnames(data[trait2])
  fit <- lm(data[,trait1]~data[,trait2])
  (ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1]))
    + geom_point(alpha=0.5, color="red")
    + geom_smooth(method="lm", color="black")
    + theme_classic()
    + labs(x=xaxis, y=yaxis,
           title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         "Intercept =",signif(fit$coef[[1]],5 ),
                         " Slope =",signif(fit$coef[[2]], 5),
                         " P =",signif(summary(fit)$coef[2,4], 5))))
}
compare.traits(residual.data, 1, 2)

test.cor <- cor(residual.data, method="pearson", use = "complete.obs")
res1 <- cor.mtest(residual.data, conf.level = .95)
corrplot(test.cor, method="circle", type="lower", order="alphabet",
         addgrid.col=NA, tl.col="black", title="unburnt continuous",
         tl.cex=0.5, p.mat = res1$p, insig = "label_sig", pch.col = "white")
# more corrplot info at: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html   
# or a PCA value, which is reducing the dimensionality of multiple traits
## if you want to do a PCA value (which you likely do, follow the steps below)
pca.test <- residual.data # make sure to remove incomplete taxa
pca.test <- traits
#test.data <- pca.test[,2:length(trait.data)] # designate the variables of interest (columns)
#test.data <- log(test.data) # log transform data
#ln.test[ln.test=="-Inf"] <- 0
#species <- pca.test[,1] # make note of the species names
#genus <- pca.test[,7]
#parity <- pca.test[,8]
trait.pca <- prcomp(pca.test) # perform your PCA
plot(trait.pca, type="l") # visualize how many axes to keep
trait.pca # determine how much influence each axis holds
summary(trait.pca) # and the cumulative value of each
#(ggbiplot(protea.pca, obs.scale=1, var.scale=1, ellipse=T, circle=T))
#loadings <- as.data.frame(protea.pca$rotation)
axes <- predict(trait.pca, newdata = pca.test) # pull out each taxon's value for each PC axis
varanus.pc1 <- axes[,1] # save the PC loadings as a new variable for the DTT
varanus.pc2 <- axes[,2]
(ggbiplot(trait.pca, obs.scale=1, var.scale=1, ellipse=T, circle=T))

      
geography <- CreateGeoObject_SP(tree, varanus.dist)
      
modelBM <- createModel(tree, "BM")
modelOU <- createModel(tree, "OU")
modelACDC <- createModel(tree, "ACDC")
modelPM <- createModel(tree, "PM")
modelPM_geo <- createGeoModel(tree, geography, "PM+geo")
modelPMOUless <- createModel(tree, "PM_OUless")
modelMC <- createModel_MC(tree)
modelMC_geo <- createGeoModel(tree, geography, "MC+geo")
modelMC_other <- createModel(tree, "MC")

fitBM <- fitTipData(modelBM, varanus.pc1, GLSstyle=T)
fitOU <- fitTipData(modelOU, varanus.pc1, GLSstyle=T)
#fitACDC <- fitTipData(modelACDC, varanus.pc2, GLSstyle=T)
fitPM0 <- fitTipData(modelPM, varanus.pc1, GLSstyle=T)
  fitPM <- fitTipData(modelPM, varanus.pc1, GLSstyle=T, params0=fitPM0$inferredParams)
  fitTipData(modelPM, varanus.pc1, GLSstyle=T, params0=c(0.1, 0.1, mean(varanus.pc1), 0.1, -0.1, 0.1))
fitPM_geo0 <- fitTipData(modelPM_geo, varanus.pc1, GLSstyle=T)
  fitPM_geo <- fitTipData(modelPM_geo, varanus.pc1, GLSstyle=T, params0=fitPM_geo0$inferredParams)
    simulateTipData(modelPM_geo, params=fitPM_geo$inferredParams, method=2)
fitPMOUless <- fitTipData(modelPMOUless, varanus.pc1, GLSstyle=1)
fitMC <- fitTipData(modelMC, varanus.pc1, GLSstyle=T)
fitMC_geo <- fitTipData(modelMC_geo, varanus.pc1, GLSstyle=T)
  fitTipData(modelMC_geo, varanus.pc1, GLSstyle=T, params0=c(0.1,0.1,-0.1))
  simulateTipData(modelMC_geo, params=fitMC_geo$inferredParams, method = 2)
fitMC_other <- fitTipData(modelMC_other, varanus.pc1, GLSstyle=T)

multiphy.AIC("fit", tree, models=c("BM", "OU", "PM", "PM_geo", "PMOUless", "MC", "MC_geo"))
#





