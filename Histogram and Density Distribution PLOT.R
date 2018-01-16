library(ggplot2)

###########################################
## Doing these plots for a single group
###########################################

## Read in your data
model.res <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Pygopodoidea.Model.Comparison.csv")

## Plot the timing of shifts as a frequency histogram
(ggplot(model.res, aes(x=timing))
       + geom_histogram(binwidth=1))

## Plot the timing of shifts as a density distribution (with mean estimate)
(ggplot(model.res, aes(x=timing)) 
  + geom_density(fill="green", alpha=0.2)
  #+ geom_vline(aes(xintercept=mean(timing, na.rm=T)),
  #             color="red", linetype="dashed", size=1)
  + scale_x_reverse(limits=c(20,0))) # use this to reverse then define the limits of the x axis

## Plot the timing of shifts as both a frequency histogram AND a density distribution
(ggplot(model.res, aes(x=timing))
  + geom_histogram(aes(y=..density..), binwidth=1, colour="black", fill="white")
  + geom_density(alpha=0.2, fill="green")
  + scale_x_reverse(limits=c(20,0))) # use this to reverse then define the limits of the x axis



###########################################
## Doing these plots for multiple groups
###########################################

## Read in your data
model.mult.res <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/All.Radiations.Model.Comparison.csv")

## Plot the timing of shifts as a density distribution, across groups defined in the "radiation" column
(ggplot(model.mult.res, aes(timing, fill=radiation, colour=radiation))
  + geom_density(alpha=0.2) # use alpha to change the opacity
  + scale_x_reverse(limits=c(20,0))) # use this to reverse then define the limits of the x axis






