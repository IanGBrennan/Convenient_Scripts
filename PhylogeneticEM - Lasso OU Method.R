library(geiger)
library(PhylogeneticEM)

#### Read in the data and the tree
#morph.data <- read.csv("/Users/Ian/Desktop/Pygo.test.csv", header=T, row.names=1)
morph.data <- read.csv("/Users/Ian/Desktop/Limbless.Data.csv", header=T)
  morph.data <- subset(morph.data, All.data=="Yes" & For.Analysis=="Yes")
    #morph.data <- subset(morph.data, Phylogeny=="Serpentes") # if you want to limit analysis to a specific group
    morph.data <- subset(morph.data, Phylogeny=="Oz.Blindsnakes" | Phylogeny=="Elapidae" | Phylogeny=="Henophidia" | Phylogeny=="Colubridae" | Phylogeny=="Viperids") # if you want to subset more than one group use '|'
      data <- as.data.frame(morph.data[,c("HL.Trunk", "TL.Trunk", "HW.Trunk")])
        #data <- as.data.frame(morph.data[,"HL.Trunk"])
          rownames(data) <- morph.data$Name.in.Tree
      

#tree <- read.nexus("BT.Australian.Marsupials.tre")
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")
tree <- trees[[1]]

# drop tips from the tree to match the data you trimmed above
all <- tree$tip.label
keep <- morph.data$Name.in.Tree
drop <- setdiff(all, keep)
new.tree <- drop.tip(tree, tip=drop)

name.check(new.tree, data)
#####################
### If you want to trim taxa from a larger tree
trees <- read.nexus("PB.Pygopodoidea.100.trees")
tree <- trees[[1]]
all <- tree$tip.label
keep <- rownames(morph.data)
drop <- setdiff(all, keep)
new.tree <- drop.tip(tree, drop)
name.check(new.tree, morph.data)
tree <- new.tree
#####################

morph.data <- t(data) # transpose the data (taxa as columns, traits as rows)
#rownames(morph.data) <- "TL.Trunk" #adjust this to name the trait if you'd like
rownames(morph.data) <- c("HL.Trunk", "TL.Trunk", "HW.Trunk") #adjust this to name the trait if you'd like
res<- PhyloEM(phylo=new.tree, 
              Y_data=morph.data,      # read in the trait data         
              process="scOU",         # identify the process to analyse    
              #random.root=T,         #    
              K_max=30,               # set a maximum limit on the number of shifts to search
              check.tips.names=T,     # check to make sure names/trait data match          
              parallel_alpha=T,       # we want to parallelize the analysis        
              Ncores=8)               # with how many cores?
              independent=F)          # if using multiple traits, are they independent?   
plot(res, show.tip.label=T, label_cex=0.1)

plot(res.head.body.tail.scOU, show.tip.label=T, label_cex=0.1)
plot(res.pygo.head.tail, show.tip.label=T, label_cex=0.5)
plot(res.head, show.tip.label=T, label_cex=0.1)
plot(res.head.tail, show.tip.label=T,label_cex=0.1)
plot(res.head.tail.scOU, show.tip.label=T,label_cex=0.1)
plot(res.tail, show.tip.label=T, label_cex=0.1)
save(res.head,
     res.head.tail,
     res.head.tail.scOU,
     res.pygo.head.tail.scOU,
     res.snakes.head.tail.scOU,
     res.skink.head.tail.scOU,
     res.tail.scOU,
     res.head.scOU,
     file="PhyloEM.all.results.RData")
dataz <- load(file="PhyloEM.all.results.RData")


res.param <- params_process(res.pygo.head.tail, method="LINselect") # pull out parameter estimates
shift.edges <- res.param$shifts$edges # pull out the placement of the shifts
res.head.simm <- shifts_to_simmap(tree, shift.edges) # translate the shifts into a Simmap Tree
plotSimmap(res.head.simm, fsize=0.2) # plot the Simmap Tree

save(res.head.simm, res.head.tail.simm, file="PhyloEM.all.SIMMAP.RData")

phy.bm   <- OUwie(phy.input[[4]], phy.input[[3]], model="BM1",  root.station=T, simmap.tree=T)


plot(equivalent_shifts(tree, res.param)) # plot the alternative shift regimes

resids <- residuals(res.head.tail)
plot(resids[1,])
plot(resids[2,])

###### Working through handling the Limbless Squamate Data
morph.data <- read.csv("/Users/Ian/Desktop/Limbless.Data.csv", header=T) # read in the whole CSV
morph.data <- subset(morph.data, All.data=="Yes") # filter it, keeping only taxa that have all the data
morph.data <- subset(morph.data, For.Analysis=="Yes") # filter it for taxa we want to include in the analysis (might have data, but not in tree)
rownames(morph.data) <- morph.data$Name.in.Tree
data <- morph.data[,c("HL.Trunk", "TL.Trunk")] # set it so the data is just the head/tail:trunk ratios
#data <- morph.data[,"HL.Trunk"]
#names(data) <- morph.data$Name.in.Tree
data <- t(data)
rownames(data) <- "HL.Trunk"

keep <- morph.data$Name.in.Tree # collect the names of taxa which we have all the data for (matches data above)
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")
tree <- trees[[1]]
all <- tree$tip.label # collect ALL the names from the tree
drop <- setdiff(all, keep) # filter the tree names to match the data names
new.tree <- drop.tip(tree, tip=drop) # drop those that don't match
tree <- new.tree # set the new tree
# Now up to the top and run the analysis!

res <- PhyloEM(phylo=tree, 
               Y_data=data,      # read in the trait data         
               process="OU",           # identify the process to analyse    
               random.root=T,          #    
               K_max=10,               # set a maximum limit on the number of shifts to search
               check.tips.names=T,     # check to make sure names/trait data match          
               parallel_alpha=T,       # we want to parallelize the analysis        
               Ncores=8,               # with how many cores?
               independent=T)          # if using multiple traits, are they independent?   
plot(res)
plot.equivalent_shifts

params.res <- params_process(res, method="LINselect")
equivalent_shifts(res, params.res)

