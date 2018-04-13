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
  for (k in 1:length(models)) {
    call.model <- paste0(prefix, models[k])
    resLIK <- -(get(call.model)$value)
    resPAR <- length(get(call.model)$inferredParams)
    resAIC <- grabAIC(resLIK, resPAR)
    resAICc <- phy.AICc(phylo, resLIK, resPAR)
    results <- rbind(results, c(resLIK, resPAR, resAIC, resAICc))
  }
  rownames(results) <- models
  colnames(results) <- c("logLik", "no.Param", "AIC", "AICc")
  results <- as.data.frame(results)
  best.model <- subset(results, AICc == min(results$AICc))
  best <- rownames(best.model); estimates <- get(paste0(prefix, best))$inferredParams
  output <- list(results, best.model, estimates)
  names(output) <- c("results", "best.model", "parameter.estimates")
  return(output)
}