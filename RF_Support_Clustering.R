library(phytools)
library(RPANDA)

enviro.data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Data/Andrae_S1.csv", header=T)
flux.data <-   read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/FossilUncertainty/Data/Aeolian_Flux.csv", header=T)
data(InfTemp)



library(deeptime); library(gridExtra)

pp <- ggplot(enviro.data, aes(Age)) +
  geom_ribbon(aes(ymin = C4_recon_lower, ymax = C4_recon_upper), fill = "pink") +
  geom_line(aes(y = C4_recon_mean), color="red") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0,80), expand = FALSE) 

qq <- gggeo_scale(pp, dat="epochs")


rr <- ggplot(flux.data, aes(Age)) +
  geom_ribbon(aes(ymin = A_Flux-35, ymax = A_Flux+35), fill = "light blue") +
  geom_line(aes(y = A_Flux), color="blue") + scale_x_reverse() + theme_classic() +
  coord_cartesian(xlim = c(0, 13), ylim = c(0,150), expand = FALSE) 

ss <- gggeo_scale(rr, dat="epochs")

grid.arrange(qq, ss, nrow=4, ncol=2)

#sampled.trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/MODEL110_Sampled_Run2_Macropodinae.trees")
sampled.trees <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/REAL_Run5_AllSchemes_508trees_Macropodinae.trees")
fossil.trees <-  read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/REAL_Run4_Fossil_519trees_Macropodinae.trees")
hypsodonty.index <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/Macropodinae_Hypsodonty_Data.csv", header=T)


# Trim tree and data down to overlapping taxa (for SAMPLED TREES)
overlaps <-   intersect(sampled.trees[[1]]$tip.label, hypsodonty.index$Taxon)
tip.drops <-    setdiff(sampled.trees[[1]]$tip.label, overlaps)
sampled.trees <- lapply(sampled.trees, drop.tip, tip=tip.drops)
sampled.data <- filter(hypsodonty.index, Taxon %in% overlaps)
sampled.HI <- sampled.data[,2]; names(sampled.HI) <- sampled.data[,1]; geiger::name.check(sampled.trees[[1]], sampled.HI) 

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/MDS_Clustering_Source.R")
class(sampled.trees) <- "multiPhylo"
inputs <- unroot(sampled.trees)
initial <- topclustMDS(inputs, mdsdim=2, max.k=10, makeplot = T, criterion = "max", trdist = "score")

##### We want to be able to compare differences among trees to differences in model fit
##### so we'll use a few different metrics to try and get at what causes discrepancies

# Calculate the Pairwise Robinson-Foulds distances among all trees
source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/MDS_Clustering_Source.R")
testo <- pairwise.RF(sampled.trees, measure="score")

# Or alternatively use the Quartet Dissimilarity method which is more sensitive
library(Quartet)
pairwise.TQ <- function(input.trees){
  pTQ <- as.data.frame(TQDist(input.trees)); colnames(pTQ) <- rownames(pTQ)
  TQtable <- NULL
  for(i in 1:nrow(pTQ)){
    for(j in i:ncol(pTQ)){
      TQtable <- rbind(TQtable, data.frame(tree1=i, tree2=j, TQdist=pTQ[i,j],
                                           AGEdist=abs(max(nodeHeights(input.trees[[i]])) - max(nodeHeights(input.trees[[j]])))))
    }
  }
  return(list(TQmatrix=pTQ, TQtable=TQtable))
}
chub <- pairwise.TQ(sampled.trees)

# Calculate the 

# From the model fitting data, extract differences in model preference among trees
all.aics <- readRDS("/Users/Ian/Desktop/FINAL_SAMPLED_Model_Fitting_AICCs.RDS")
all.aics <- readRDS("/Users/Ian/Desktop/FINAL_FOSSIL_Model_Fitting_AICCs.RDS")

pairwise.AIC <- function(input.AIC, target.model, comparison.table){
  AICtable <- NULL
  #timer <- progress_estimated(length(tree.span))
  
  #for(i in 1:length(unique(input.AIC$tree))){
  #  for(j in 1:length(unique(input.AIC$tree))){
  #    rowi <- filter(input.AIC, tree==i & model==target.model)
  #    rowj <- filter(input.AIC, tree==j & model==target.model)
  #    AICtable <- rbind(AICtable, data.frame(tree1=i, tree2=j, 
  #                                           AICCWdiff=abs(rowi$aiccw - rowj$aiccw)))
  #  }
  #  print(paste("finished column", i))
  #}
  
  for(i in 1:nrow(comparison.table)){
    rowi <- filter(input.AIC, tree==comparison.table[i,"tree1"] & model==target.model)
    rowj <- filter(input.AIC, tree==comparison.table[i,"tree2"] & model==target.model)
    AICtable <- rbind(AICtable, data.frame(tree1=rowi$tree, tree2=rowj$tree, 
                                           AICCWdiff=abs(rowi$aiccw - rowj$aiccw)))
    print(paste("finished row", i))
    #print(timer$tick())
  }
  
  return(AICtable)
}
testa <- pairwise.AIC(all.aics, target.model="ENVexp", chub$TQtable)

# Now combine the different information together into a single dataframe (remove comparisons of a tree with itself)
besto <- left_join(chub$TQtable, testa); besto <- filter(besto, !tree1==tree2)
#besta <- besto

# Plot AICcWeight difference as a result of Age difference
agedist <- ggplot(besto, aes(x=AGEdist, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#F21A00", shape=16) +
  geom_smooth(col="black") + theme_classic() + geom_smooth(col="black", method="lm", se=F)

# Plot AICcWeight difference as a result of Quartet Dissimilarity or Robinson-Foulds distances
tqdist <- ggplot(besto, aes(x=TQdist, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#EBCC2A", shape=16) +
  geom_smooth(col="black") + theme_classic() + geom_smooth(col="black", method="lm", se=F)

ee.tq <- grid.arrange(agedist, tqdist, nrow=1)
grid.arrange(gl.tq, ge.tq, fe.tq, ee.tq, nrow=4)

# Plot AICcWeight difference as a result of Quartet Dissimilarity or Robinson-Foulds distances
ggplot(besto, aes(x=TQAGE, y=AICCWdiff)) +
  geom_point(alpha=0.5, col="#F21A00") +
  geom_smooth() + theme_classic()

filter(besto, TQdist==0)

# Plot the relationship between AICcWeight and Model Support
sub.res <- filter(all.aics, model %in% c("GRASSexp", "GRASSlin", "FLUXexp", "ENVexp"))
cols <- brewer.pal(5, "Paired"); names(cols) <- c("1", "2", "3", "4", "0")
col.pal <- brewer.pal(8, "Paired")[1:4]; col.pal <- rep(col.pal, max(all.aics$tree))
g.exp <- ggplot(sub.res, aes(x=age, y=aiccw)) +
  geom_point(aes(colour=col.pal)) +
  geom_smooth(aes(colour=col.pal), method="lm", se=F) + theme_classic()

grid.arrange(g.lin, g.exp, f.exp, e.exp)

# And plot the beta values for preferred models
col.pal <- brewer.pal(8, "Paired")
#beta.values <- filter(all.aics, model %in% c("GRASSexp", "FLUXexp", "GRASSlin", "ENVexp") & aiccw >= 0.3); #beta.values <- filter(beta.values, par<70 & par>-25)
ge_res <- filter(all.aics, model=="GRASSexp"); ge <- ggplot(ge_res, aes(x=age, y=aiccw)) + geom_point(colour=col.pal[[3]]) + theme_classic() + geom_smooth(colour=col.pal[[3]]) + geom_smooth(method="lm", se=F, colour=col.pal[[3]]) #+ geom_dotplot(binaxis="y", dotsize=0.1, binwidth=5)
gl_res <- filter(all.aics, model=="GRASSlin"); gl <- ggplot(gl_res, aes(x=age, y=aiccw)) + geom_point(colour=col.pal[[4]]) + theme_classic() + geom_smooth(colour=col.pal[[4]]) + geom_smooth(method="lm", se=F, colour=col.pal[[4]]) #+ geom_dotplot(binaxis="y", dotsize=0.1, binwidth=5)
fe_res <- filter(all.aics, model=="FLUXexp");  fe <- ggplot(fe_res, aes(x=age, y=aiccw)) + geom_point(colour=col.pal[[1]]) + theme_classic() + geom_smooth(colour=col.pal[[1]]) + geom_smooth(method="lm", se=F, colour=col.pal[[1]]) #+ geom_dotplot(binaxis="y", dotsize=0.1, binwidth=5)
ee_res <- filter(all.aics, model=="ENVexp");   ee <- ggplot(ee_res, aes(x=age, y=aiccw)) + geom_point(colour=col.pal[[2]]) + theme_classic() + geom_smooth(colour=col.pal[[2]]) + geom_smooth(method="lm", se=F, colour=col.pal[[2]]) #+ geom_dotplot(binaxis="y", dotsize=0.1, binwidth=5)

fossil.res <- gridExtra::grid.arrange(gl, ge, fe, ee, nrow=1)
grid.arrange(sampled.res, fossil.res, nrow=2)


#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
ggplotRegression <- function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.25, color="red") + # change to 0.25 and "red" for time plots
    stat_smooth(method = "lm", col = "black") + # change to "black" for time plots
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
#######################################################################################
fit <- lm(AICCWdiff ~ TQdist, data=besto) # change this according to the parameter you simulated
plot.fit <- (ggplotRegression(fit))
plot.fit + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

# extract trees which had high support for the FLUXexp model
hi.AICS <- filter(all.aics, model=="GRASSexp" & aiccw >= 0.5)
hi.trees <- sampled.trees[hi.AICS$tree]; class(hi.trees) <- "multiPhylo"
hi.chub <- filter(chub$TQtable, tree1 %in% hi.AICS$tree & tree2 %in% hi.AICS$tree)
hi.testa <- pairwise.AIC(hi.AICS, target.model="FLUXexp", hi.chub)
hi.besto <- left_join(hi.chub, hi.testa); hi.besto <- filter(hi.besto, !tree1==tree2)


gvf <- filter(all.aics, model=="GRASSexp" | model=="FLUXexp")
fiftyAIC <- filter(gvf, aiccw>=0.50)

(ggplot(fiftyAIC, aes(x=age, colour=model))
                    #+ geom_density(alpha=0.75, adjust=0.5)
                    #+ geom_histogram(alpha=0.75)
                    + geom_freqpoly(bins=50)
                    + theme_classic()
                    + scale_fill_manual(values=wes_palette("Zissou1", type="continuous", 3))
                    + scale_x_reverse(lim=c(15,5)))

##################
# Test the relationship between Gamma Stat and AICcWt for Grass Model

gstat <- NULL
for(k in 1:max(all.aics$tree)){
  gstat <- append(gstat, rep(gammaStat(sampled.trees[[k]]),9))
}
all.aics$gamma <- gstat

sub.res <- filter(all.aics, model %in% c("GRASSexp", "GRASSlin", "FLUXexp", "ENVexp"))
sub.res <- filter(testo, model == "GRASSlin")
#cols <- brewer.pal(5, "Paired"); names(cols) <- c("1", "2", "3", "4", "0")
#col.pal <- brewer.pal(8, "Paired")[1:4]; col.pal <- rep(col.pal, max(all.aics$tree))
f.g.lin <- ggplot(sub.res, aes(x=gamma, y=aiccw)) +
  geom_point(colour= col.pal[[4]]) +
  geom_smooth(colour=col.pal[[4]], method="lm", se=F) + 
  geom_smooth(colour=col.pal[[4]], method="loess") +
  theme_classic()
grid.arrange(g.lin, g.exp, f.exp, e.exp, 
             f.g.lin, f.g.exp, f.f.exp, f.e.exp, nrow=2)

gstat <- NULL


fit <- lm(aiccw ~ gamma, data=sub.res) # change this according to the parameter you simulated
plot.fit <- (ggplotRegression(fit))
plot.fit + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))


