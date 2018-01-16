setwd("C:/Workspace/king0393/Placoderms/rates_R_final")
source("get.epoch.rates.R")
source("geoplot.epoch.rates.R")
library(OutbreakTools)
library(picante)
library(geoscale)
read.annotated.nexus("placoderms.trees") -> trees
# epoch boundaries should be in ages, not height (converted to height internally)
# The oldest age can make an effectively unbounded epoch by setting this to being older than the oldest sampled root age
# epochs do not have to be all the way to the root, can cut off the uncertain early period if desired
# However note that the plot function sets the age of each epoch at the midpoint of the boundaries. To get around this can rewrite the epochs vector before plotting, as here
# The youngest age should be exactly the age of the youngest taxon
# Can set any number of epochs
c(600, 444, 419, 408, 393, 383, 359, 290) -> epochs
get.epoch.rates(trees, epochs, burnin=0.1) -> run1
# The plot function has four options. the first two are the output from get.epoch.rates and the epoch bounds
# line.col is the colour of the lines. Default black
# plot.mean is whether or not to plot the mean across the whole posterior sample in each epoch. Default true
# mean.col is the colour of the mean plot. Default red
# plot.bounds is whether or not to add dotted vertical lines at the boundaries of each epoch. Default true
# all other arguments are passed to geoscalePlot. See ?geoscalePlot for options
# can rewrite the epochs vector with a new oldest boundary for plotting
c(460, 444, 419, 408, 393, 383, 359, 290) -> epochs
geoplot.epoch.rates(run1, epochs, line.col="grey", plot.mean=T, mean.col="black", plot.bounds=T, age.lim=c(470,280), cex.age=0.8, cex.ts=1, units=c("Age", "Period"), tick.scale=20, ts.col=F)


c(444, 419, 408, 393, 383, 359, 290) -> epochs
get.epoch.rates(trees, epochs, burnin=0.1) -> run1
geoplot.epoch.rates(run1, epochs, line.col="grey", plot.mean=T, mean.col="black", plot.bounds=T, age.lim=c(440,290), cex.age=0.8, cex.ts=1, units=c("Age", "Period"), tick.scale=20, ts.col=F)

