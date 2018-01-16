library(l1ou)
library(phytools)
library(parallel)
library(geiger)

#### Read in the data and the tree
#morph.data <- read.csv("BT.Australian.Marsupials.male.BL.MASS.csv", header=T, row.names=1)
morph.dataz <- read.csv("/Users/Ian/Desktop/Pygo.test.csv", header=T, row.names=1)
tree <- read.nexus("BT.Australian.Marsupials.tre")

#### Check to make sure the tips match the data labels
name.check(new.tree, morph.data);

#### Adjust the data and tree to fit (order matters in l1ou!)
data <- adjust_data(tree, morph.data)

#### Estimate the number and position of shifts a priori 
shift.fit <- estimate_shift_configuration(data$tree, data$Y, nCores=8, quietly=F, criterion="pBIC") # if you want to do only a single trait, 'data$Y[,x]'
plot(shift.fit, cex=0.5)

#### Investigate convergence among the supported shift set
fit.conv <- estimate_convergent_regimes(shift.fit, nCores=8, criterion="pBIC")
plot(fit.conv, cex=0.2)

#### If you need to estimate fit under different scoring criterion
pBIC.score <- configuration_ic(data$tree, data$Y, shift.fit$shift.configuration,
                               criterion="pBIC")

#### Translate the shift values into actual scores that relate to your data
length.regions <- convert_shifts2regions(data$tree, shift.fit$shift.configuration, shift.fit$shift.values[,1]) + shift.fit$intercept[1]
mass.regions <- convert_shifts2regions(data$tree, shift.fit$shift.configuration, shift.fit$shift.values[,2]) + shift.fit$intercept[2]



###### Delete below this line

###### Working through handling the Limbless Squamate Data
morph.data <- read.csv("/Users/Ian/Desktop/Limbless.Data.csv", header=T) # read in the whole CSV
morph.data <- subset(morph.data, All.data=="Yes") # filter it, keeping only taxa that have all the data
morph.data <- subset(morph.data, For.Analysis=="Yes") # filter it for taxa we want to include in the analysis (might have data, but not in tree)
rownames(morph.data) <- morph.data$Name.in.Tree
data <- morph.data[,c("HL.Trunk", "TL.Trunk")] # set it so the data is just the head/tail:trunk ratios
morph.data <- data

keep <- morph.data$Name.in.Tree # collect the names of taxa which we have all the data for (matches data above)
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.PreTrimmed.Limbless.100.trees")
tree <- trees[[1]]
all <- tree$tip.label # collect ALL the names from the tree
drop <- setdiff(all, keep) # filter the tree names to match the data names
new.tree <- drop.tip(tree, tip=drop) # drop those that don't match
tree <- new.tree # set the new tree
# Now up to the top and run the analysis!


test[test == "fossorial"] <- "0"; test[test == "semifossorial"] <- "0"
test[test == "terrestrial"] <- 1; test[test == "arboreal"] <- 1; test[test == "semiarboreal"] <- 1
test[test == "aquatic"] <- 1



