library(phytools)
library(OUwie)
library(mvMORPH)
library(parallel)
library(RPANDA)
library(geiger)
library(Rmisc)
library(dplyr)
library(ggplot2); library(scales)
library(wesanderson)
library(randomcoloR)
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/Source_pulsR.R") # source pulsR code


trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")
datas <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Meliphagoids.MASS.csv", header=T, row.names=1)
    data.mv <- datas[,1]; names(data.mv) <- rownames(datas); data.mv <- log(data.mv); name.check(trees[[1]], datas)

# Files:
    # Pygopodoidea: PB.Pygopodoidea.100.trees; BT.Pygopodoidea.SVL.csv
    # Marsupials: PB.Australian.Marsupials.100.trees; BT.Australian.Marsupials.BL.csv
    # Agamids: PB.Agamids.100.trees; BT.Agamids.SVL.csv
    # Meliphagid Birds: PB.Meliphagides.100.trees; BT.Meliphagoids.MASS.csv
    # Skinks: PB.Skinks.100.trees; BT.Skinks.SVL.csv
    
save.path <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/"

#############################################################
multiphylo.era.simmap <- function(multiphy, split.time) {
  split.trees <- NULL
  for (q in 1:length(multiphy)){
    slices <- 0
    for (y in 1:length(split.time)) {
      temp.slice <- max(nodeHeights(multiphy[[q]])) - split.time[y]
      slices <- append(slices, temp.slice)
    }
    split.trees[[q]] <- make.era.map(multiphy[[q]], slices)
  }
  return(split.trees)
}
#############################################################
# e.g: era.trees <- multiphylo.era.simmap(trees, 5)

#############################################################
two.era.multiphylo <- function(multiphy, split.time) {
  split.trees <- mclapply(1:length(multiphy), function(x) 
    {make.era.map(multiphy[[x]], c(0, (max(nodeHeights(multiphy[[x]]))-split.time)))})
}
#############################################################
# e.g: era.trees <- two.era.multiphylo(trees, 5)

all.models.lik <- NULL
model.fit <- NULL

group.name <- "Meliphagoids"
slice.i <- 4; slice.j <- 20 # designate the beginning and end of your shift times
era.trees.list <- mclapply(slice.i:slice.j, function(x){
  two.era.multiphylo(trees, x)}, mc.cores=8)

# if you're using empirical estimates of measurement error or intraspecific variation:
######################################################################################

      # Fit and Process the OU model
      OU1_list <- mclapply(1:length(trees), function(x){
        mvOU(trees[[x]], data.mv, model="OU1", method="sparse", error=NA, # error=data.me^2
             diagnostic=F, echo=F)}, mc.cores = 8)
      OU_res <- NULL; for(y in 1:length(OU1_list)) {
        OU_temp <- as.data.frame(t(c(OU1_list[[y]]$LogLik,OU1_list[[y]]$AIC,OU1_list[[y]]$AICc, y, "OU1", NA, group.name)))
        OU_res <- rbind(OU_res, OU_temp)
      }
      all.models.lik <- rbind(all.models.lik, OU_res)
      model.fit[["OU1"]] <- OU1_list
      
      # Fit and Process the  BM model
      BM1_list <- mclapply(1:length(trees), function(x){
        mvBM(trees[[x]], data.mv, model="BM1", method="sparse", error=data.me^2,
             diagnostic=F, echo=F)}, mc.cores = 8)
      BM_res <- NULL; for(y in 1:length(BM1_list)) {
        BM_temp <- as.data.frame(t(c(BM1_list[[y]]$LogLik,BM1_list[[y]]$AIC,BM1_list[[y]]$AICc, y, "BM1", NA, group.name)))
        BM_res <- rbind(BM_res, BM_temp)
      }
      all.models.lik <- rbind(all.models.lik, BM_res)
      model.fit[["BM1"]] <- BM1_list
      
      # Fit and Process the EB model
      EB_list <- mclapply(1:length(trees), function(x){
        mvEB(trees[[x]], data.mv, method="sparse", error=data.me^2,
             diagnostic=F, echo=F)}, mc.cores = 8)
      EB_res <- NULL; for(y in 1:length(EB_list)) {
        EB_temp <- as.data.frame(t(c(EB_list[[y]]$LogLik,EB_list[[y]]$AIC,EB_list[[y]]$AICc, y, "EB", NA, group.name)))
        EB_res <- rbind(EB_res, EB_temp)
      }
      all.models.lik <- rbind(all.models.lik, EB_res)
      model.fit[["EB"]] <- EB_list
      
      # Fit and Process the BMOU model
      BMOU_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(era.trees.list[[x]]), function(y){
          mvSHIFT(era.trees.list[[x]][[y]], data.mv, model="BMOU", method="sparse",
                  error=data.me^2, diagnostic=F, echo=F)}, mc.cores = 8)})
      BMOU_list <- NULL
      for (k in 1:length(BMOU_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(BMOU_int_list)){
          int.best[[i]] <- BMOU_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        BMOU_list[[k]] <- int.best[[1]]
      }
      BMOU_res <- NULL; for(y in 1:length(BMOU_list)) {
        BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", BMOU_list[[y]]$shift.time, group.name)))
        BMOU_res <- rbind(BMOU_res, BMOU_temp)
      }
      all.models.lik <- rbind(all.models.lik, BMOU_res)
      model.fit[["BMOU"]] <- BMOU_list
      
      # Fit and Process the BMOUi model
      BMOUi_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(era.trees.list[[x]]), function(y){
          mvSHIFT(era.trees.list[[x]][[y]], data.mv, model="BMOUi", method="sparse",
                  error=data.me^2, diagnostic=F, echo=F)}, mc.cores = 8)})
      BMOUi_list <- NULL
      for (k in 1:length(BMOUi_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(BMOUi_int_list)){
          int.best[[i]] <- BMOUi_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        BMOUi_list[[k]] <- int.best[[1]]
      }
      BMOUi_res <- NULL; for(y in 1:length(BMOUi_list)) {
        BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", BMOUi_list[[y]]$shift.time, group.name)))
        BMOUi_res <- rbind(BMOUi_res, BMOUi_temp)
      }
      all.models.lik <- rbind(all.models.lik, BMOUi_res)
      model.fit[["BMOUi"]] <- BMOUi_list
      
      # Fit and Process the EBOU model
      EBOU_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(era.trees.list[[x]]), function(y){
          mvSHIFT(era.trees.list[[x]][[y]], data.mv, model="EBOU", method="sparse",
                  error=data.me^2, diagnostic=F, echo=F)}, mc.cores = 8)})
      EBOU_list <- NULL;
      for (k in 1:length(EBOU_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(EBOU_int_list)){
          int.best[[i]] <- EBOU_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        EBOU_list[[k]] <- int.best[[1]]
      }
      EBOU_res <- NULL; for(y in 1:length(EBOU_list)) {
        EBOU_temp <- as.data.frame(t(c(EBOU_list[[y]]$LogLik,EBOU_list[[y]]$AIC,EBOU_list[[y]]$AICc, y, "EBOU", EBOU_list[[y]]$shift.time, group.name)))
        EBOU_res <- rbind(EBOU_res, EBOU_temp)
      }
      all.models.lik <- rbind(all.models.lik, EBOU_res)
      model.fit[["EBOU"]] <- EBOU_list
      
      # Fit and Process the EBOUi model
      EBOUi_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(era.trees.list[[x]]), function(y){
          mvSHIFT(era.trees.list[[x]][[y]], data.mv, model="EBOUi", method="sparse",
                  error=data.me^2, diagnostic=F, echo=F)}, mc.cores = 8)})
      EBOUi_list <- NULL
      for (k in 1:length(EBOUi_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(EBOUi_int_list)){
          int.best[[i]] <- EBOUi_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        EBOUi_list[[k]] <- int.best[[1]]
      }
      EBOUi_res <- NULL; for(y in 1:length(EBOUi_list)) {
        EBOUi_temp <- as.data.frame(t(c(EBOUi_list[[y]]$LogLik,EBOUi_list[[y]]$AIC,EBOUi_list[[y]]$AICc, y, "EBOUi", EBOUi_list[[y]]$shift.time, group.name)))
        EBOUi_res <- rbind(EBOUi_res, EBOUi_temp)
      }
      all.models.lik <- rbind(all.models.lik, EBOUi_res)
      model.fit[["EBOUi"]] <- EBOUi_list
      
      # Fit and Process the JN (jump normal) model
      JN_list <- mclapply(1:length(trees), function(x){
        fit_reml_levy(trees[[x]], data.mv, model="JN", maxiter=10, sigma_tip=T, silent=T)}, mc.cores = 8)
      
      JN_res <- NULL; for(y in 1:length(JN_list)) {
        JN_temp <- as.data.frame(t(c(JN_list[[y]]$lnL, JN_list[[y]]$AIC, (-2*(JN_list[[y]]$lnL) + 2*JN_list[[y]]$n_params*(length(JN_list[[y]]$dat)/(length(JN_list[[y]]$dat)-(JN_list[[y]]$n_params)-1))),
                                    y, "JN", NA, group.name))) # we have to manually calculate the AICc for the 'pulsR' models
        JN_res <- rbind(JN_res, JN_temp)
      }
      all.models.lik <- rbind(all.models.lik, JN_res)
      model.fit[["JN"]] <- JN_list
      
      # Fit and Process the NIG (non-Gaussian Jump) model
      NIG_list <- mclapply(1:length(trees), function(x){
        fit_reml_levy(trees[[x]], data.mv, model="NIG", maxiter=10, sigma_tip=T, silent=T)}, mc.cores = 8)
      
      NIG_res <- NULL; for(y in 1:length(NIG_list)) {
        NIG_temp <- as.data.frame(t(c(NIG_list[[y]]$lnL, NIG_list[[y]]$AIC, (-2*(NIG_list[[y]]$lnL) + 2*NIG_list[[y]]$n_params*(length(NIG_list[[y]]$dat)/(length(NIG_list[[y]]$dat)-(NIG_list[[y]]$n_params)-1))),
                                     y, "NIG", NA, group.name))) # we have to manually calculate the AICc for the 'pulsR' models
        NIG_res <- rbind(NIG_res, NIG_temp)
      }
      all.models.lik <- rbind(all.models.lik, NIG_res)
      model.fit[["NIG"]] <- NIG_list
      
      # Fit and Process the Environmental model
      data(InfTemp)
      ENV_list <- mclapply(1:length(trees), function(x){
        fit_t_env(trees[[x]], data.mv, error=data.me^2, env_data=InfTemp, df=50, scale=F, plot=F)}, mc.cores = 8)
      ENV_res <- NULL; for(y in 1:length(ENV_list)) {
        ENV_temp <- as.data.frame(t(c(ENV_list[[y]]$LH,ENV_list[[y]]$aic,ENV_list[[y]]$aicc, y, "ENV", NA, group.name)))
        ENV_res <- rbind(ENV_res, ENV_temp)
      }
      all.models.lik <- rbind(all.models.lik, ENV_res)
      model.fit[["ENV"]] <- ENV_list
      
      # Fit and Process the BioGeoBEARS model
      Allo.mars <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Allopatry.Trend.Australian.Marsupials.csv", header=T)
      BGB_list <- mclapply(1:length(trees), function(x){
        fit_t_env(trees[[x]], data.mv, error=data.me^2, env_data=Allo.mars, df=50, scale=F, plot=F)}, mc.cores = 8)
      BGB_res <- NULL; for(y in 1:length(BGB_list)) {
        BGB_temp <- as.data.frame(t(c(BGB_list[[y]]$LH,BGB_list[[y]]$aic,BGB_list[[y]]$aicc, y, "BGB", NA, group.name)))
        BGB_res <- rbind(BGB_res, BGB_temp)
      }
      all.models.lik <- rbind(all.models.lik, BGB_res)
      model.fit[["BGB"]] <- BGB_list
      
      # Fit and Process the BMS model
      data.ou <- as.data.frame(names(data.mv)); data.ou[,2] <- data.mv; data.ou[,3] <- data.me
      BMS_int_list <- lapply(slice.i:slice.j, function(x){
        mclapply(1:length(trees), function(y){
          OUwie.slice(trees[[y]], data.ou[,c(1:3)],model="BMS", mserr="known",
                      root.station=T, timeslices=c(x))}, mc.cores=8)})
      BMS_list <- NULL;
      for (k in 1:length(BMS_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(BMS_int_list)){
          int.best[[i]] <- BMS_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        BMS_list[[k]] <- int.best[[1]]
      }
      BMS_res <- NULL; for(y in 1:length(BMS_list)) {
        BMS_temp <- as.data.frame(t(c(BMS_list[[y]]$loglik, BMS_list[[y]]$AIC, BMS_list[[y]]$AICc, y, "BMS", BMS_list[[y]]$shift.time, group.name)))
        BMS_res <- rbind(BMS_res, BMS_temp)
      }
      all.models.lik <- rbind(all.models.lik, BMS_res)
      model.fit[["BMS"]] <- BMS_list
      
      # Fit and Process the OUM model
      data.ou <- as.data.frame(names(data.mv)); data.ou[,2] <- data.mv; data.ou[,3] <- data.me
      OUM_int_list <- lapply(slice.i:slice.j, function(x){
        mclapply(1:length(trees), function(y){
          OUwie.slice(trees[[y]], data.ou[,c(1:3)],model="OUM", mserr="known",
                      root.station=T, timeslices=c(x))}, mc.cores=8)})
      OUM_list <- NULL;
      for (k in 1:length(OUM_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(OUM_int_list)){
          int.best[[i]] <- OUM_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$AICc))]
        OUM_list[[k]] <- int.best[[1]]
      }
      OUM_res <- NULL; for(y in 1:length(OUM_list)) {
        OUM_temp <- as.data.frame(t(c(OUM_list[[y]]$loglik, OUM_list[[y]]$AIC, OUM_list[[y]]$AICc, y, "OUM", OUM_list[[y]]$shift.time, group.name)))
        OUM_res <- rbind(OUM_res, OUM_temp)
      }
      all.models.lik <- rbind(all.models.lik, OUM_res)
      model.fit[["OUM"]] <- OUM_list
      
######################################################################################
      

# if you're estimating measurement error or intraspecific variation:
######################################################################################
      # Fit and Process the OU model
      OU_list <- mclapply(1:length(trees), function(x){
        fitContinuous(trees[[x]], as.data.frame(data.mv), model="OU", SE=NA)}, mc.cores = 8)
      OU_res <- NULL; for(y in 1:length(OU_list)) {
        OU_temp <- as.data.frame(t(c(OU_list[[y]]$opt$lnL,OU_list[[y]]$opt$aic,OU_list[[y]]$opt$aicc, y, "OU", NA, group.name)))
        OU_res <- rbind(OU_res, OU_temp)
      }
      all.models.lik <- rbind(all.models.lik, OU_res); unique(all.models.lik$V5)
      model.fit[["OU"]] <- OU_list; names(model.fit)

      # Fit and Process the BM model
      BM_list <- mclapply(1:length(trees), function(x){
        fitContinuous(trees[[x]], as.data.frame(data.mv), model="BM", SE=NA)}, mc.cores = 8)
      BM_res <- NULL; for(y in 1:length(BM_list)) {
        BM_temp <- as.data.frame(t(c(BM_list[[y]]$opt$lnL,BM_list[[y]]$opt$aic,BM_list[[y]]$opt$aicc, y, "BM", NA, group.name)))
        BM_res <- rbind(BM_res, BM_temp)
      }
      all.models.lik <- rbind(all.models.lik, BM_res); unique(all.models.lik$V5)
      model.fit[["BM"]] <- BM_list; names(model.fit)
      
      # Fit and Process the EB model
      EB_list <- mclapply(1:length(trees), function(x){
        fitContinuous(trees[[x]], as.data.frame(data.mv), model="EB", SE=NA)}, mc.cores = 8)
      EB_res <- NULL; for(y in 1:length(BM_list)) {
        EB_temp <- as.data.frame(t(c(EB_list[[y]]$opt$lnL,EB_list[[y]]$opt$aic,EB_list[[y]]$opt$aicc, y, "EB", NA, group.name)))
        EB_res <- rbind(EB_res, EB_temp)
      }
      all.models.lik <- rbind(all.models.lik, EB_res); unique(all.models.lik$V5)
      model.fit[["EB"]] <- EB_list; names(model.fit)
      
      # Fit and Process the SRC model (use in place of BMOU if estimating ME)
      SRC_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(trees), function(y){
          fitContinuous_paleo(trees[[y]], as.data.frame(data.mv), model="SRC_ME", shift.time=x)}, mc.cores = 8)})
      SRC_list <- NULL
      for (k in 1:length(SRC_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(SRC_int_list)){
          int.best[[i]] <- SRC_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$Trait1$aicc))]
        SRC_list[[k]] <- int.best[[1]]
      }
      SRC_res <- NULL; for(y in 1:length(SRC_list)) {
        SRC_temp <- as.data.frame(t(c(SRC_list[[y]]$Trait1$lnl,SRC_list[[y]]$Trait1$aic,SRC_list[[y]]$Trait1$aicc, y, "SRC", SRC_list[[y]]$shift.time, group.name)))
        SRC_res <- rbind(SRC_res, SRC_temp)
      }
      all.models.lik <- rbind(all.models.lik, SRC_res); unique(all.models.lik$V5)
      model.fit[["SRC"]] <- SRC_list; names(model.fit)
      
      # Fit and Process the TRC model (use in place of BMOUi if estimating ME)
      TRC_int_list <- lapply(1:length(era.trees.list), function(x){
        mclapply(1:length(trees), function(y){
          fitContinuous_paleo(trees[[y]], as.data.frame(data.mv), model="TRC_ME", shift.time=x)}, mc.cores = 8)})
      TRC_list <- NULL
      for (k in 1:length(TRC_int_list[[1]])) {
        int.best <- NULL
        for (i in 1:length(TRC_int_list)){
          int.best[[i]] <- TRC_int_list[[i]][[k]] 
          int.best[[i]]["shift.time"] <- slice.i + (i-1)
        }
        int.best <- int.best[order(sapply(int.best, function(x) x$Trait1$aicc))]
        TRC_list[[k]] <- int.best[[1]]
      }
      TRC_res <- NULL; for(y in 1:length(TRC_list)) {
        TRC_temp <- as.data.frame(t(c(TRC_list[[y]]$Trait1$lnl,TRC_list[[y]]$Trait1$aic,TRC_list[[y]]$Trait1$aicc, y, "TRC", TRC_list[[y]]$shift.time, group.name)))
        TRC_res <- rbind(TRC_res, TRC_temp)
      }
      all.models.lik <- rbind(all.models.lik, TRC_res); unique(all.models.lik$V5)
      model.fit[["TRC"]] <- TRC_list; names(model.fit)
      
      # Fit and Process the JN (jump normal) model
      JN_list <- mclapply(1:length(trees), function(x){
        fit_reml_levy(trees[[x]], data.mv, model="JN", maxiter=10, sigma_tip=T, silent=T)}, mc.cores = 8)
      
      JN_res <- NULL; for(y in 1:length(JN_list)) {
        JN_temp <- as.data.frame(t(c(JN_list[[y]]$lnL, JN_list[[y]]$AIC, (-2*(JN_list[[y]]$lnL) + 2*JN_list[[y]]$n_params*(length(JN_list[[y]]$dat)/(length(JN_list[[y]]$dat)-(JN_list[[y]]$n_params)-1))),
                                     y, "JN", NA, group.name))) # we have to manually calculate the AICc for the 'pulsR' models
        JN_res <- rbind(JN_res, JN_temp)
      }
      all.models.lik <- rbind(all.models.lik, JN_res); unique(all.models.lik$V5)
      model.fit[["JN"]] <- JN_list; names(model.fit)
      
      # Fit and Process the NIG (non-Gaussian Jump) model
      NIG_list <- mclapply(1:length(trees), function(x){
        fit_reml_levy(trees[[x]], data.mv, model="NIG", maxiter=10, sigma_tip=T, silent=T)}, mc.cores = 8)
      
      NIG_res <- NULL; for(y in 1:length(NIG_list)) {
        NIG_temp <- as.data.frame(t(c(NIG_list[[y]]$lnL, NIG_list[[y]]$AIC, (-2*(NIG_list[[y]]$lnL) + 2*NIG_list[[y]]$n_params*(length(NIG_list[[y]]$dat)/(length(NIG_list[[y]]$dat)-(NIG_list[[y]]$n_params)-1))),
                                      y, "NIG", NA, group.name))) # we have to manually calculate the AICc for the 'pulsR' models
        NIG_res <- rbind(NIG_res, NIG_temp)
      }
      all.models.lik <- rbind(all.models.lik, NIG_res); unique(all.models.lik$V5)
      model.fit[["NIG"]] <- NIG_list; names(model.fit)
      
      # Fit and Process the Environmental model
      data(InfTemp)
      ENV_list <- mclapply(1:length(trees), function(x){
        fit_t_env(trees[[x]], data.mv, error=NA, env_data=InfTemp, df=50, scale=F, plot=F)}, mc.cores = 8)
      ENV_res <- NULL; for(y in 1:length(ENV_list)) {
        ENV_temp <- as.data.frame(t(c(ENV_list[[y]]$LH,ENV_list[[y]]$aic,ENV_list[[y]]$aicc, y, "ENV", NA, group.name)))
        ENV_res <- rbind(ENV_res, ENV_temp)
      }
      all.models.lik <- rbind(all.models.lik, ENV_res); unique(all.models.lik$V5)
      model.fit[["ENV"]] <- ENV_list; names(model.fit)
      
      # Fit and Process the BioGeoBEARS model
      Allo <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Allopatry.Trend.Skinks.csv", header=T)
      BGB_list <- mclapply(1:length(trees), function(x){
        fit_t_env(trees[[x]], data.mv, error=NA, env_data=Allo, df=50, scale=F, plot=F)}, mc.cores = 8)
      BGB_res <- NULL; for(y in 1:length(BGB_list)) {
        BGB_temp <- as.data.frame(t(c(BGB_list[[y]]$LH,BGB_list[[y]]$aic,BGB_list[[y]]$aicc, y, "BGB", NA, group.name)))
        BGB_res <- rbind(BGB_res, BGB_temp)
      }
      all.models.lik <- rbind(all.models.lik, BGB_res); unique(all.models.lik$V5)
      model.fit[["BGB"]] <- BGB_list; names(model.fit)
      
      
######################################################################################
      
save.all.liks <- paste0(save.path, group.name, "_ModelFit.csv")
    write.csv(all.models.lik, file=save.all.liks, row.names=F)
        # unique(all.models.lik$V5); length(unique(all.models.lik$V5))
        # filter(all.models.lik, V4==49)

  # names(model.fit); length(names(model.fit))
save.all.models <- paste0(save.path, group.name, "_ModelObjects.RDS")
    saveRDS(model.fit, file = save.all.models)
        # model.fit <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Pygopodoidea_ModelObjects.RDS")   
  

total.results <- all.models.lik
#total.results <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Pygopodoidea_ModelFit.csv", header=T)
colnames(total.results) <- c("lnL", "AIC", "AICc", "tree.no", "model", "shift.time", "group")
total.results$tree.no <- as.numeric(total.results$tree.no) # set the tree numbers as actual numbers
sorted.first <- subset(total.results, total.results$tree.no < 10) # subset 1-9
sorted.first <- arrange(sorted.first, tree.no) # sort this subset
sorted.second <- subset(total.results, total.results$tree.no > 9) # subset 10-100
sorted.second <- arrange(sorted.second, tree.no) # sort this subset, this is done to avoid counting errors
sorted.results <- rbind.data.frame(sorted.first, sorted.second) # now bind them together
weight.delta <- NULL
for (j in 1:100){
  targetdata <- dplyr::filter(sorted.results, tree.no==j)
  res <- aicw(as.numeric(as.character(targetdata$AICc)))
  weight.delta <- rbind.data.frame(weight.delta, res)
}
weight.delta <- within(weight.delta, rm(fit))
sorted.results <- cbind.data.frame(sorted.results, weight.delta)
#sub <- subset(sorted.results, sorted.results$tree.no < 101)
outz <- summarySE(sorted.results, measurevar="w", groupvars="model")

outz$group <- "agamid lizards"
FINAL.TABLE <- rbind.data.frame(FINAL.TABLE, outz)
  save.all.liks <- paste0(save.path, "Agamidae", "_ModelFit.csv")
    write.csv(sorted.results, file=save.all.liks, row.names=F)
# unique(all.models.lik$V5); length(unique(all.models.lik$V5))
# filter(all.models.lik, V4==49)

## Otherwise Just Plot the model results
########################################
outz$model <- factor(outz$model, 
                     levels=c("BM", "OU", "EB",
                              "JN", "NIG",
                              "ENV", "BGB",
                              "SRC", "TRC"))# re-order the models in the plot
myplot <- (ggplot(outz, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + geom_text(aes(label=percent(w), vjust=-4)))
  #+ theme(axis.text.x=element_text(angle=45, hjust=1))
  #+ theme_classic())
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.text.x=element_text(angle=45, hjust=1))

# if you just want to do a single plot
myplot <- (ggplot(outz)
           + geom_bar(aes(1, y=w, fill=outz$model), stat="identity")
           + scale_fill_manual( values=wes_palette("Zissou", 9, "continuous")))
myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))


skinks <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Sphenomorphines_ModelFit.csv", header=T)
mars <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupial_ModelFit.csv", header=T)
pygo <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea_ModelFit.csv", header=T)
bird <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Meliphagoids_ModelFit.csv", header=T)
agam <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Agamidae_ModelFit.csv", header=T)

output.total <- rbind.data.frame(skinks, mars, pygo, bird, agam)
colnames(output.total) <- c("lnL", "AIC", "AICc", "tree.no", "model", "shift.time", "group")
# now we need to resolve the AICc weights again, and delta
weight.delta <- NULL
for (k in 1:length(unique(output.total$group))) {
  group.name <- unique(output.total$group)
  target.group <- subset(output.total, output.total$group == group.name[[k]])
  weights <- NULL
  for (j in 1:length(target.group$tree.no)) {
    target.data <- subset(target.group, target.group$tree.no ==j)
    res <- aicw(as.numeric(target.data$AICc))
    weights <- rbind.data.frame(weights, res)
  }
  weight.delta <- rbind.data.frame(weight.delta, weights)
}
weight.delta <- within(weight.delta, rm(fit))
output.total <- cbind.data.frame(output.total, weight.delta)    

## If you want to plot a composite bar graph, subset each radiation, summarize model weights
outz.mars <- summarySE(subset(output.total, group=="Marsupials"), measurevar="w", groupvars="model")
outz.agam <- summarySE(subset(output.total, group=="Agamidae"), measurevar="w", groupvars="model"); #outz.agam[6,c(3:6)]<-.001
outz.pygo <- summarySE(subset(output.total, group=="Pygopodoidea"), measurevar="w", groupvars="model")
outz.skink <- summarySE(subset(output.total, group=="Sphenomorphines"), measurevar="w", groupvars="model")
outz.bird <- summarySE(subset(output.total, group=="Meliphagoids"), measurevar="w", groupvars="model")
## Then combine them back together
outz.mars["group"] <- "marsupial mammals"
outz.agam["group"] <- "agamidae lizards"
outz.pygo["group"] <- "pygopodoid geckos"
outz.skink["group"] <- "sphenomorphine skinks"
outz.bird["group"] <- "meliphagoid birds"
group.outz <- rbind.data.frame(outz.mars, outz.agam, outz.pygo, outz.skink, outz.bird)

## Plot the composite bar graphs
#group.outz$model <- factor(group.outz$model, levels=c("BM", "delta", "kappa", "lambda", "gamma", "EB", "DensDep", "OU", "BMS", "TS", "OUS", "OUMA", "OUMV", "OUMVA", "SRC", "TRC")) # this re-orders the models in the legend
group.outz$model <- factor(group.outz$model, levels=c("BM", "OU", "EB",
                                                      "JN", "NIG",
                                                      "ENV", "BGB",
                                                      "SRC", "TRC")) # this re-orders the models in the legend
myplot <- (ggplot(group.outz)
           + geom_bar(aes(y=w, x=group, fill=group.outz$model), stat="identity")
           + theme(axis.text.x=element_text(angle=25, hjust=1), panel.background=element_blank(), legend.position="bottom")
           + scale_fill_manual( values=wes_palette("Zissou", 9, "continuous")))


output.total$tree.no <- as.numeric(output.total$tree.no) # set the tree numbers as actual numbers
sorted.first <- subset(output.total, output.total$tree.no < 10) # subset 1-9
sorted.first <- arrange(sorted.first, tree.no) # sort this subset
sorted.second <- subset(output.total, output.total$tree.no > 9) # subset 10-100
sorted.second <- arrange(sorted.second, tree.no) # sort this subset, this is done to avoid counting errors
sorted.results <- rbind.data.frame(sorted.first, sorted.second) # now bind them together


bird.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Meliphagoids_ModelObjects.RDS")
pygo.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea_ModelObjects.RDS")
agam.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Agamidae_ModelObjects.RDS")
skink.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Sphenomorphines_ModelObjects.RDS")
mars.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupial_ModelObjects.RDS")


plot.all <- function(res.object, model=NULL, ylim=NULL, xlim=NULL, group=NULL){
  beta.est <- NULL
  if(is.null(model)){
    stop ("which kind of model?")
  }
  else if (model=="ENV"){
    plot(res.object$ENV[[1]], ylim=ylim, xlim=xlim)
    for(y in 2:length(res.object$ENV)){
      lines(res.object$ENV[[y]], col=randomColor(1))
      beta.est <- append(beta.est, res.object$ENV[[y]]$param[[2]])
    }
  }
  else if (model=="BGB"){
    plot(res.object$BGB[[1]], ylim=ylim, xlim=xlim)
    for(y in 2:length(res.object$BGB)){
      lines(res.object$BGB[[y]], col=randomColor(1))
      beta.est <- append(beta.est, res.object$BGB[[y]]$param[[2]])
    }
  }
  beta.est <- as.data.frame(beta.est)
  beta.est[,2] <- group
  return(beta.est)
}
# par(mfrow=c(3,2)) # if you want to multiplot

skink.beta <- plot.all(skink.res, model="ENV", group="skink")
mars.beta <- plot.all(mars.res, model="ENV", group="marsupials")
agam.beta <- plot.all(agam.res, model="ENV", group="agamids")
pygo.beta <- plot.all(pygo.res, model="ENV", group="pygopodoids")
bird.beta <- plot.all(bird.res, model="ENV", group="meliphagoids")
all.beta <- rbind.data.frame(skink.beta, mars.beta, agam.beta, pygo.beta, bird.beta)
      all.beta <- all.beta[which(all.beta[,1] > -2),] # remove the crazy low (or high) estimates
colnames(all.beta) <- c("beta", "group")

(ggplot(all.beta, aes(x=group, y=beta, fill=group)) 
  + geom_violin()
  + theme(axis.text.x=element_text(angle=25, hjust=1), panel.background=element_blank(), legend.position="bottom"))

# Plot all the models with their appropriate weights
aic.test <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/All.AICw.RDS")
(ggplot(aic.test, aes(x=model, y=w, fill=model))
  + geom_bar(stat="identity")
  + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
  + scale_fill_manual( values=wes_palette("Zissou1", 9, "continuous"))
  + theme(panel.background=element_blank(), legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1))
  + facet_wrap(~group, nrow=3, ncol=2)
  + geom_text(aes(label=percent(w), vjust=-1)))


#
