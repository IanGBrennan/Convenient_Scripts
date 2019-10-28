#######################################################################################
# Interlude: the base 'plot' and 'abline' functions are alright, but we want to 
## make it (1) prettier, and (2) include the information from our linear regression
### into the plot, so that we know what our results were. Use custom 'ggplotRegression'
### if you want to change the saturation use 'alpha'
# function inputs:
## 'fit': is the linear model object
## 'alpha': is the opacity of the points
## 'size': the size of the points
## 'color': is a vector of the same length as the point data, providing plot colors
ggplotRegression <- function (fit, alpha=0.5, size=1, color=NULL) {
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(alpha=alpha, color=color, size=size) +
    stat_smooth(method = "lm", col = "black") + 
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))) + theme_classic()
}
#######################################################################################