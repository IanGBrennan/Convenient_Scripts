#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
ggplotRegression <- function (fit) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=0.75, color="dark Green", size=3) + # change to 0.25 and "red" for time plots
    stat_smooth(method = "lm", col = "black") + # change to "black" for time plots
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + 
    theme_classic() + geom_abline(intercept = 0, slope = 1, color="red", linetype=2)
    
}
#######################################################################################
# the order: agam, bird, mars, pygo, skink


CoEvo.res <- readRDS("/Users/Ian/Desktop/SRecovery_CoEvo.RDS")
CoEvo.res$model <- "CoEvo"

CoEvoall.res <- readRDS("/Users/Ian/Desktop/SRecovery_CoEvoall.RDS")
CoEvoall.res$model <- "CoEvo_all"

CoPMgeo.res <- readRDS("/Users/Ian/Desktop/SRecovery_CoPMgeo.RDS")
CoPMgeo.res$model <- "CoPM_geo"

JointPMgeo.res <- readRDS("/Users/Ian/Desktop/SRecovery_JointPMgeo.RDS")
JointPMgeo.res$model <- "JointPM_geo"

CoEvosplit.res <- readRDS("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/SRecovery_CoEvosplit.RDS")
CoEvosplit.res$model <- "CoEvo_Split"

all.res <- NULL; all.res[[1]] <- CoEvo.res; all.res[[2]] <- CoEvoall.res; all.res[[3]] <- CoPMgeo.res
    all.res[[4]] <- JointPMgeo.res[,1:2]; all.res[[5]] <- JointPMgeo.res[,3:4]; 
    all.res[[1]] <- CoEvosplit.res[,1:2]; all.res[[2]] <- CoEvosplit.res[,3:4];
    colnames(all.res[[4]])<-c("simulated","estimated"); colnames(all.res[[5]])<-c("simulated","estimated")
    colnames(all.res[[1]])<-c("simulated","estimated"); colnames(all.res[[2]])<-c("simulated","estimated")
    

# write loop to automate this process, otherwise you're dumb
plot.res <- NULL
for (i in 1:length(all.res)) { # change this according to the parameter you simulated
  fit <- lm(estimated ~ simulated, data=all.res[[i]]) # change this according to the parameter you simulated
  plot.fit <- (ggplotRegression(fit))
  plot.res[[i]] <- plot.fit
}

jointPMgeo.plot <- grid.arrange(plot.res[[4]], plot.res[[5]], nrow=1)
coevoSplit.plot <- grid.arrange(plot.res[[1]], plot.res[[2]], nrow=1)

grid.arrange(plot.res[[1]], plot.res[[2]],
             plot.res[[3]], jointPMgeo.plot,
             plot.res[[1]],
             nrow=3)

grid.arrange(coevoSplit.plot,coevoSplit.plot,
             coevoSplit.plot,coevoSplit.plot,
             coevoSplit.plot, nrow=3)











