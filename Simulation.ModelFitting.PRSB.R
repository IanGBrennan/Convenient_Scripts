######################################################################################
# Loop 4:
# Simulate data under the BMtrend process onto the extinct trees, then fit models!
######################################################################################
# or read in fossilized trees you've already made
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Meliphagides.100.trees")

trees <- miocene
# trees <- multiphylo.era.simmap(trees, shift.time) # if simulating with mvMORPH!

bm.pars <- 0.1 # set the diffusion parameter of the BM process
trend.pars <- sample(seq(from=0.1, to=0.5, by=0.01), 100, replace=T)
shift.time <- 10







# if simulating using Phytools for BMtrend
simulated.traits <- mclapply(1:length(trees), function(x){
   (fastBM(trees[[x]], mu=trend.pars[[x]], sig2=bm.pars, a=0, nsim=1))})
# if simulating using mvMORPH for BMOUi or BMOU
simulated.traits <- mclapply(1:length(trees), function(x){
  mvSIM(trees[[x]], nsim=1, model="BMOUi", param=test)}, mc.cores=8)
sim.traits.env <- NULL; for(h in 1:length(simulated.traits)){
  sim.traits.env[[h]] <- simulated.traits[[h]][,1]; names(sim.traits.env[[h]]) <- rownames(simulated.traits[[h]])
} 
# if simulating under the ENV model of RPANDA
simulated.traits <- mclapply(1:length(trees), function(x){
  sim_t_env(trees[[1]], param=c(bird.res$ENV[[x]]$param[[1]], bird.res$ENV[[x]]$param[[2]]), env_data=InfTemp, model="EnvExp")})
# make sure the trees are SIMMAP if using BMOU or BMOUi models!


#run.simulations <- function(multiphy, group, nsim) {}

all.models.lik <- NULL
model.fit <- NULL

group.name <- "BMOU"
save.path <- "/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/"

# Fit and Process the OU model
OU1_list <- mclapply(1:length(trees), function(x){
  mvOU(trees[[x]], simulated.traits[[x]], model="OU1", method="sparse", # error=data.me^2
       diagnostic=F, echo=F)}, mc.cores = 8)
OU_res <- NULL; for(y in 1:length(OU1_list)) {
  OU_temp <- as.data.frame(t(c(OU1_list[[y]]$LogLik,OU1_list[[y]]$AIC,OU1_list[[y]]$AICc, y, "OU1", NA, group.name)))
  OU_res <- rbind(OU_res, OU_temp)
}
all.models.lik <- rbind(all.models.lik, OU_res)
model.fit[["OU"]] <- OU1_list

# Fit and Process the  BM model
BM1_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], simulated.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
BM_res <- NULL; for(y in 1:length(BM1_list)) {
  BM_temp <- as.data.frame(t(c(BM1_list[[y]]$LogLik,BM1_list[[y]]$AIC,BM1_list[[y]]$AICc, y, "BM1", NA, group.name)))
  BM_res <- rbind(BM_res, BM_temp)
}
all.models.lik <- rbind(all.models.lik, BM_res)
model.fit[["BM"]] <- BM1_list

# Fit and Process the EB model
EB_list <- mclapply(1:length(trees), function(x){
  mvEB(trees[[x]], simulated.traits[[x]], method="sparse",
       diagnostic=F, echo=F)}, mc.cores = 8)
EB_res <- NULL; for(y in 1:length(EB_list)) {
  EB_temp <- as.data.frame(t(c(EB_list[[y]]$LogLik,EB_list[[y]]$AIC,EB_list[[y]]$AICc, y, "EB", NA, group.name)))
  EB_res <- rbind(EB_res, EB_temp)
}
all.models.lik <- rbind(all.models.lik, EB_res)
model.fit[["EB"]] <- EB_list

# Fit and Process the ENV model, don't use mclapply, as it gets gummed up somehow
data(InfTemp)
ENV_list<-NULL; for(j in 1:length(trees)){
  ENV_list[[j]] <- fit_t_env(trees[[j]], sim.traits.env[[j]], env_data=InfTemp, df=50, scale=F, plot=F)
}
ENV_res <- NULL; for(y in 1:length(ENV_list)) {
  ENV_temp <- as.data.frame(t(c(ENV_list[[y]]$LH,ENV_list[[y]]$aic,ENV_list[[y]]$aicc, y, "ENV", NA, group.name)))
  ENV_res <- rbind(ENV_res, ENV_temp)
}
all.models.lik <- rbind(all.models.lik, ENV_res)
model.fit[["ENV"]] <- ENV_list

# Fit and Process the BMOU model
BMOU_list <- mclapply(1:length(trees), function(x){
  mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-shift.time))), simulated.traits[[x]], model="BMOU", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
  #mvSHIFT(trees[[x]], simulated.traits[[x]], model="BMOU", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
BMOU_res <- NULL; for(y in 1:length(BMOU_list)) {
  #BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", cheese.time[[y]], BMOU_list[[y]]$shift.time, group.name)))
  BMOU_temp <- as.data.frame(t(c(BMOU_list[[y]]$LogLik,BMOU_list[[y]]$AIC,BMOU_list[[y]]$AICc, y, "BMOU", shift.time, BMOU_list[[y]]$shift.time, group.name)))
  BMOU_res <- rbind(BMOU_res, BMOU_temp)
}
all.models.lik <- rbind(all.models.lik, BMOU_res)
model.fit[["BMOU"]] <- BMOU_list

# Fit and Process the BMOUi model
BMOUi_list <- mclapply(1:length(trees), function(x){
  mvSHIFT(make.era.map(trees[[x]], c(0, (max(nodeHeights(trees[[x]]))-shift.time))), simulated.traits[[x]], model="BMOUi", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
  #mvSHIFT(trees[[x]], simulated.traits[[x]], model="BMOUi", method="sparse",diagnostic=F, echo=F)}, mc.cores = 8)
BMOUi_res <- NULL; for(y in 1:length(BMOUi_list)) {
  #BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", cheese.time[[y]], BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_temp <- as.data.frame(t(c(BMOUi_list[[y]]$LogLik,BMOUi_list[[y]]$AIC,BMOUi_list[[y]]$AICc, y, "BMOUi", shift.time, BMOUi_list[[y]]$shift.time, group.name)))
  BMOUi_res <- rbind(BMOUi_res, BMOUi_temp)
}
all.models.lik <- rbind(all.models.lik, BMOUi_res)
model.fit[["BMOUi"]] <- BMOUi_list

# Fit and Process the BMtrend model
Trend_list <- mclapply(1:length(trees), function(x){
  mvBM(trees[[x]], simulated.traits[[x]], model="BM1", method="sparse",
       diagnostic=F, echo=F, param=list(trend=TRUE))}, mc.cores = 8)
Trend_res <- NULL; for(y in 1:length(Trend_list)) {
  Trend_temp <- as.data.frame(t(c(Trend_list[[y]]$LogLik,Trend_list[[y]]$AIC,Trend_list[[y]]$AICc, y, "BMtrend", NA, group.name)))
  Trend_res <- rbind(Trend_res, Trend_temp)
}
all.models.lik <- rbind(all.models.lik, Trend_res)
model.fit[["Trend"]] <- Trend_list


save.all.liks <- paste0(save.path, group.name, "_Simulated_ModelFit.csv")
    write.csv(all.models.lik, file=save.all.liks, row.names=F)
# unique(all.models.lik$V5); length(unique(all.models.lik$V5))
# filter(all.models.lik, V4==49)

# names(model.fit); length(names(model.fit))
save.all.models <- paste0(save.path, group.name, "_Simulated_ModelObjects.RDS")
    saveRDS(model.fit, file = save.all.models)
# model.fit <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Pygopodoidea_ModelObjects.RDS")   


total.results <- all.models.lik
#total.results <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/BMOU_Simulated_ModelFit.csv", header=T)
colnames(total.results) <- c("lnL", "AIC", "AICc", "tree.no", "model", "shift.time", "group")
total.results$tree.no <- as.numeric(total.results$tree.no) # set the tree numbers as actual numbers
sorted.first <- subset(total.results, total.results$tree.no < 10) # subset 1-9
sorted.first <- arrange(sorted.first, tree.no) # sort this subset
sorted.second <- subset(total.results, total.results$tree.no > 9) # subset 10-100
sorted.second <- arrange(sorted.second, tree.no) # sort this subset, this is done to avoid counting errors
sorted.results <- rbind.data.frame(sorted.first, sorted.second) # now bind them together
weight.delta <- NULL
for (j in 1:10){
  targetdata <- dplyr::filter(sorted.results, tree.no==j)
  res <- geiger::aicw(as.numeric(as.character(targetdata$AICc)))
  weight.delta <- rbind.data.frame(weight.delta, res)
}
weight.delta <- within(weight.delta, rm(fit))
sorted.results <- cbind.data.frame(sorted.results, weight.delta)
#sub <- subset(sorted.results, sorted.results$tree.no < 101)
outz <- summarySE(sorted.results, measurevar="w", groupvars="model")


## Otherwise Just Plot the model results
########################################
outz$model <- factor(outz$model, levels=c("BM1", "EB", "OU1", "ENV", "BMOU", "BMOUi", "BMtrend"))

myplot <- (ggplot(outz, aes(x=model, y=w, fill=model))
           + geom_bar(stat="identity")
           + geom_errorbar(aes(ymin=w-se, ymax=w+se), size=0.3, width=0.2)
           + geom_text(aes(label=percent(w), vjust=-2))
           + scale_fill_manual(values=wes_palette("Zissou1", 7, "continuous")))
#+ theme(axis.text.x=element_text(angle=45, hjust=1))
#+ theme_classic())
BMOU.sim <-  myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               axis.text.x=element_text(angle=45, hjust=1))
multiplot(BM.sim, OU.sim,
          EB.sim, BMOU.sim,
          BMOUi.sim, ENV.sim,
          trend.sim, cols=2)

# done: 
