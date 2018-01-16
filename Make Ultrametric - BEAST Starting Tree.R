library(ape)
# http://www.justinbagley.org/1226/update-new-functions-for-generating-starting-trees-for-beast-or-starbeast-in-r

# create a starting tree for BEAST from an ML/BI tree, this will get the analysis past the logLik=-INF stage

tree<-read.tree('Elapidae.concatenated.tre')
tree<-root(tree, c("Epictia_tenella_h0", "Epictia_tenella_h1"))
plot(tree)
#################################################################
#### Make an Ultrametric Tree from ML/BI using PL in APE
#################################################################

# Here are examples of each of the above options, assuming we have the treefile input above and we are interested in 
# using one fossil calibration ranging between 17-14 Ma (million years ago), which we want to place on node 156 of our 
# phylogeny. Because our tree file is large and the data are non-clocklike, we will use a relaxed clock model, except in the 
# 3rd and 4th examples we will use a special tweak to run a strict clock model. We will also specify the bounds of our 
# calibration to be ‘hard’ rather than ‘soft’ bounds (though this appears it may be unused in the current version of ape, so 
# not sure about that). The output trees will be named timetreex, where x = {1,2,3,4}, corresponding to the above set 
# {Option #1, Option #2, Option #3, Option #4}.
plot(tree, cex=0.2); nodelabels(cex=0.4, frame="circle", bg="yellow");
plot(tree, cex=0.2); tiplabels(cex=0.4, frame="circle", bg="yellow");
oz.hyd<-drop.tip(tree, c(1:12)) #remove non-Australian beetles
write.tree(oz.hyd, file="Oz.Hydroporini.tre")

tree<-read.tree("Oz.Hydroporini.tre")

##Option #1. 
timetree1 <- chronos(tree, lambda = 0, model = "relaxed", quiet = FALSE, calibration = makeChronosCalib(tree, node = "132", age.min = 28, age.max = 32, interactive = FALSE, soft.bounds = TRUE)) 
plot(timetree1)
write.tree(timetree1, file="Oz.Beet.tre")

##Option #2
cali<-read.csv("hydroporini.calibrations.csv", header=TRUE) #read in a set of calibrations if you've got em
abs_calib <- makeChronosCalib(tree, interactive = TRUE, soft.bounds = FALSE) #or designate them interactively
timetree3 <- chronos(tree, lambda = 0, model = "relaxed", calibration = cali)
write.tree(timetree3, file="Oz.Beet3.tre")

##Option #3 (making it a strict clock run by specifying a single rate category, and a discrete model)
strict_ctrl <- chronos.control(nb.rate.cat = 1)
timetree1 <- chronos(fishtree, model = "discrete", calibration = makeChronosCalib(fishtree, node = "156", age.min = 14, age.max = 17, interactive = FALSE, soft.bounds = FALSE), control = strict_ctrl) 
##Option #4
timetree1 <- chronos(fishtree, model = "discrete", calibration = makeChronosCalib(fishtree, node = "156", age.min = 14, age.max = 17, interactive = FALSE, soft.bounds = FALSE), control = chronos.control(nb.rate.cat = 1)) 

##Making my 3-calibration time tree
tree146_3Calib_chronogram <- chronos(tree, lambda = 0, model = "relaxed", calibration = tree146_3Calib_cal_absolute)
##Plotting the new time tree
plot.phylo(tree146_3Calib_chronogram, type = 'phylogram', use.edge.length=TRUE, font=3, cex=0.2, no.margin=TRUE)

plot(tree, cex=0.2); nodelabels(cex=0.3, frame="none")
cali<-read.csv("Python.dates.csv")
timepython<-chronos(tree, lambda=0, model="relaxed", calibration=cali)
write.tree(timepython, "Python.starter.tre")

#####################################
#### Stupid Version for TreeSnatcher
####################################
##Option #2
tree<-read.tree("Hydroporini.tre")
branching.times(cstree) #get the branching times, put into text wrangler, and build your csv below (no decimal points)
read.csv("Crayfish.dates.csv")
calib<-read.csv("Hydroporini.dates.csv", header=TRUE) #read in a set of calibrations if you've got em
cali<- makeChronosCalib(cstree, interactive = TRUE, soft.bounds = FALSE) #or designate them interactively
timetree3 <- chronos(cstree, lambda = 0, model = "correlated", calibration = calib)
branching.times(timetree3)
write.tree(timetree3, file="oz.hydroporini.tre")

fix(calib)
?collapse.singles
cstree<-collapse.singles(tree)
plot(cstree)
branching.times(cstree)

