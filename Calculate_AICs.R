require(phytools)

grabAIC <- function (LnL, numparams) {
  AICval = 2 * numparams - 2 * LnL
  return(AICval)
}


phy.AICc <- function(phylo, LnL, numparams) {
  # Calculate AIC
  AICval = 2*numparams - 2*LnL
  
  # Correction for finite sample size
  if(class(phylo)=="multiPhylo") { samplesize = Ntip(phylo[1])+Ntip(phylo[2])}
  else {samplesize = Ntip(phylo)}
  correction_val = (2*numparams*(numparams+1)) / (samplesize - numparams - 1)
  
  AICc_val = AICval + correction_val
  
  return(AICc_val)
}


multiphy.AIC <- function(prefix, phylo, models) {
  results <- NULL
  model.names <- NULL
  for (k in 1:length(models)) {
    if (models[k] == "GMM" || models[k] == "GMM0" || models[k] == "GMM$best.result"
        || models[k] == "GMM_geo" || models[k] == "GMM_geo0" || models[k] == "GMM_all"
        || models[k] == "PM" || models[k] == "PM_geo" || models[k] == "CoPM_geo" || models[k] == "CoPM"
        || models[k] == "JointPM_geo" || models[k] == "MC" || models[k] == "MC_geo"
        || models[k] == "BM" || models[k] == "OU" 
        || models[k] == "CoEvo" || models[k] == "CoEvo_all" || models[k] == "CoEvo_Split"
        || models[k] == "CoEvo_Split$best.result" || models[k] == "JointPM" || models[k] == "JointPM_geo"
        || models[k] == "PMOU" || models[k] == "ACDC") {
      call.model <- paste0(prefix, models[k])
      resLIK <- -(get(call.model)$value)
      resPAR <- length(get(call.model)$inferredParams)
      resAIC <- grabAIC(resLIK, resPAR)
      resAICc <- phy.AICc(phylo, resLIK, resPAR)
      results <- rbind(results, c(resLIK, resPAR, resAIC, resAICc))
      model.names <- append(model.names, models[k])
    } else if (models[k] == "rbt.BM" || models[k] == "rbt.OU") {
      call.model <- paste0(prefix, models[k])
      
      resLIK1 <- get(call.model)$multi.rate.model$logL
      resPAR1 <- get(call.model)$multi.rate.model$k
      resAIC1 <- grabAIC(resLIK1, resPAR1)
      resAICc1 <- phy.AICc(phylo, resLIK1, resPAR1)
      results <- rbind(results, c(resLIK1, resPAR1, resAIC1, resAICc1))
      model.names <- append(model.names, paste0(models[k],"_ind"))
      
      resLIK2 <- get(call.model)$common.rate.model$logL
      resPAR2 <- get(call.model)$common.rate.model$k
      resAIC2 <- grabAIC(resLIK2, resPAR2)
      resAICc2 <- phy.AICc(phylo, resLIK2, resPAR2)
      results <- rbind(results, c(resLIK2, resPAR2, resAIC2, resAICc2))
      model.names <- append(model.names, paste0(models[k], "_com"))
    }
  }
  rownames(results) <- model.names
  colnames(results) <- c("logLik", "no.Param", "AIC", "AICc")
  results <- as.data.frame(results)
  results <- results[order(results$AICc),]
  results$delta.AICc <- results$AICc - results$AICc[1]
  results$AICcWt <- (as.matrix(aic.w(results$AICc)))
  best.model <- subset(results, AICc == min(results$AICc))
  best <- rownames(best.model); 
  if (best == "GMM" || best == "GMM0" || best == "GMM_geo" || best == "GMM_geo0" || best == "GMM_all"
      || best == "PM_geo" || best == "CoPM_geo" || best == "CoPM" || best == "JointPM_geo"
      || best == "CoEvo" || best == "CoEvo_all" || best == "CoEvo_Split"
      || best == "CoEvo_Split" || best == "JointPM" || best == "JointPM_geo" || best == "PMOU" || best == "ACDC"){
    estimates <- get(paste0(prefix, best))$inferredParams
    comment <- get(paste0("model", best))@comment
  } else if (best == "PM" || best == "PM_geo" || best == "MC" || best == "MC_geo"
             || best == "OU" || best == "BM") {
    estimates <- get(paste0(prefix, best))$inferredParams
    
  } else if (best == "rbt.OU_ind" || best == "rbt.OU_com" || best == "rbt.BM_ind" || best == "rbt.BM_com") {
    if (best == "rbt.OU_ind" || best == "rbt.BM_ind") {
      estimates <- get(strsplit(paste0(prefix, best),"_")[[1]][1])$multi.rate.model[1:6]
    } else if (best == "rbt.OU_com" || best == "rbt.BM_com") {
      estimates <- get(strsplit(paste0(prefix, best),"_")[[1]][1])$common.rate.model[1:6]
    }
  }
  output <- list(results, best.model, estimates, comment)
  names(output) <- c("results", "best.model", "parameter.estimates", "model.description")
  return(output)
}


joint.fit <- function(fit1, fit2){
  combo.fit <- NULL
  combo.fit$value <- sum(fit1$value + fit2$value)
  combo.fit$inferredParams <- c(fit1$inferredParams, fit2$inferredParams)
  return(combo.fit)
}