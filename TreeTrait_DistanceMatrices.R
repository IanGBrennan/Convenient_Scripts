library(phytools)
library(DataCombine)
library(sp)
library(adehabitatHR)
library(rgeos)
library(rworldmap); library(ggmap)
library(Rmisc)
library(diversitree)
library(mvtnorm)

tree <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Pygopodoidea.tre")
  trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/PB.Pygopodoidea.100.trees")
trait <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Pygopodoidea.logSVL.csv", header=F)
#range <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Pygopodoidea.Range.Test.csv", header=T, row.names=1)
distribution <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Clean_Pygopodoidea.csv", header=T)
  distribution <- distribution[,c("Name_in_Tree", "Latitude", "Longitude")]
    distribution <- distribution[complete.cases(distribution),] 
      # as.data.frame(table(distribution$Name_in_Tree))

### We don't have distributional data for everything, so we need to drop a few from the tree
  drop <- setdiff(tree$tip.label, distribution$Name_in_Tree)
  tree <- drop.tip(tree, tip=drop)
    trees <- lapply(trees, drop.tip, tip=drop); class(trees) <- "multiPhylo"
### and from the trait data
  trait <- trait[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
    traits <- trait[,2]; names(traits) <- trait[,1] #read in data file in RPANDA format
  sim.trait <- test[which(trait[,1] %in% unique(distribution$Name_in_Tree)),]
    sim.traits <- sim.trait[,2]; names(sim.traits) <- sim.trait[]


all.TTR <- NULL   
for (pp in 1:length(trees)) {
  cat("iteration", pp, "of", length(trees), "\n") #keep track of what tree/loop# we're on
  tree <- trees[[pp]]
  
  ## Build a distance matrix for all taxa in the tree
  dist.mat <- cophenetic.phylo(tree)
  
  ## Build a distance matrix for the trait of interest for all taxa in the tree
  diff.mat <- matrix(NA,length(traits), length(traits))
  colnames(diff.mat) <- names(traits); rownames(diff.mat) <- names(traits)
  diff.mat <- as.data.frame(diff.mat)
  for (j in 1:length(tree$tip.label)) {
    for (k in 1:length(tree$tip.label)) {
      diff.mat[j,k] <- abs(traits[j]-traits[k])
    }
  }
  
  ## Now build a data frame to hold the Tree/Trait/Range distance data
  # First, determine tips in common:
  inboth <- intersect(tree$tip.label, names(traits))
  
  # Get unique combinations of this set:
  unique.combo <-  combn(inboth, m = 2)
  
  # make vectors to hold results
  dist_trait <- rep(NA, ncol(unique.combo)) # the trait distances
  dist_tree <-  rep(NA, ncol(unique.combo)) # the tree distances
  
  TTR <- data.frame(species1 = unique.combo[1,], 
                    species2 = unique.combo[2,], 
                    dist_trait, dist_tree, stringsAsFactors=F)
  
  # now fill it with the pairwise comparisons
  for (ii in 1:nrow(TTR)){
    TTR$dist_trait[ii]  <- diff.mat[TTR$species1[ii], TTR$species2[ii]]
    TTR$dist_tree[ii]   <- dist.mat[TTR$species1[ii], TTR$species2[ii]]    
  }
  TTR$dist_tree <- TTR$dist_tree/2
  ## Save each/all of the Distance/Data Matrices externally
  all.TTR[[pp]] <- TTR
}
saveRDS(all.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.Matrices.RDS")

## Next we need to determine (pairwise) if taxa overlap in their ranges
  ## this step only needs to be done once! (won't change with changes to the tree)
all.taxa <- unique(distribution$Name_in_Tree)
all.sp <- list()
all.hull <- list()
#all.poly <- list()
sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")

for (p in 202:length(unique(distribution$Name_in_Tree))) {
  current.taxon <- all.taxa[[p]]
  current.data <- subset(distribution, distribution$Name_in_Tree == current.taxon)
  if (nrow(current.data) > 1000) {
    current.data <- current.data[sample(nrow(current.data), 1000, replace=F), ]
  } 
  points <- current.data[,c(2,3)]
  all.sp[[p]] <- distribution.sp <- SpatialPoints(points)
  distribution.hull <- LoCoH.k(distribution.sp, k=length(distribution.sp)/2, duplicates="random")
  all.hull[[p]] <- gUnionCascaded(distribution.hull)
  #all.poly[[p]] <- distribution.poly <- mcp(distribution.sp, percent=100)
  
  # only plot the maps below if you need to check/clean the data
    pdf(paste("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Skink_Maps/", current.taxon, ".pdf", sep=""))
    print(ggmap(sq_map) + geom_point(data = current.data, mapping = aes(x = Longitude, y = Latitude), color="red"))
    dev.off()
}
names(all.sp) <- unique(distribution$Name_in_Tree)
names(all.hull) <- unique(distribution$Name_in_Tree)
saveRDS(all.sp,   file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.SpatialPoints.RDS")
saveRDS(all.hull, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.Polygons.RDS")

range_overlap <- rep(NA, ncol(unique.combo)) # make a ma
arrange <- data.frame(species1 = unique.combo[1,], 
                      species2 = unique.combo[2,], 
                      range_overlap, stringsAsFactors=F)
for (t in 1:length(arrange[,1])) {
  overlap <- gOverlaps(all.hull[[arrange[t,1]]], all.hull[[arrange[t,2]]])
  arrange[t,3] <- overlap
}

all.TTR <- lapply(all.TTR, cbind, arrange$range_overlap)
all.TTR <- lapply(all.TTR, setNames, c("species1", "species2", "dist_trait", "dist_tree", "range_overlap"))
saveRDS(all.TTR, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.Matrices.RDS")

curr.t <- all.taxa[[1]]
curr.d <- subset(distribution, distribution$Name_in_Tree == curr.t)
#old.poly <- all.hull$P_Strophurus_ciliaris_ciliaris; plot(old.poly)
pts <- curr.d[,c(2,3)]
test.sp <- SpatialPoints(pts)
test <- gBuffer(test.sp, width=1)
#cas.test <- gUnionCascaded(test)
#all.test <- list()
#all.test[[2]] <- test

test.f <- fortify(test)
#test.f <- merge(test, by.x)
#sbbox <- make_bbox(lon = test.f$long, lat = test.f$lat, f = .1)
#sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
(ggmap(sq_map) 
  + geom_point(data = curr.d, mapping = aes(x=Longitude, y=Latitude), color="red")
  + geom_polygon(data = test.f, aes(x = lat, y = long, group=group, alpha=0.5)))

(ggmap(sq_map) + geom_polygon(aes(x = long, y = lat, alpha=0.5), data = test.f, colour = "black"))

(ggmap(sq_map) + geom_polygon(data = test.o, aes(x = long, y = lat, group=group)))

(ggmap(sq_map) + geom_polygon(data = all.hull[[1]], aes(x = x, y = y), color="red"))


ggplot(sq_map)


## Pull out all sister species for comparison
## in a situation like Hesperoedura, this includes comparison to all:
## Oedura, Amalosia, Nebulifera (make sense?)
all.matches <- NULL # empty object
for(p in 1:length(tree$tip.label)) {
  taxon <- tree$tip.label[p] # choose a current taxon
  sisters <- getSisters(tree, node=taxon, mode="number") # get the sister node/tips to this taxon
  all.desc <- getDescendants(tree, sisters) # get all the descendants of that sister node
  tip.desc <- subset(all.desc, all.desc <= length(tree$tip.label)) # keep only the tips (terminal nodes)
  
  pair.matrix <- matrix(NA, ncol=2, nrow=(length(tip.desc))) # make an empty matrix for the pairwise comparisons
  pair.matrix[,1] <- taxon # who they're being compared against (taxon)
  for (q in 1:length(tip.desc)) {
    pair.matrix[q,2] <- tree$tip.label[tip.desc[q]] # add all the comparisons
  }
  all.matches <- rbind(all.matches, pair.matrix) # keep all the pairwise matches
}
#all.matches <- unique(all.matches[,1:2]) # drop any duplicates, which seem to happen

## Now we need to get the matching data from our Tree/Trait/Range data frame
all.sister.pairs <- NULL
for (tt in 1:length(all.TTR)) {
  cat("iteration", tt, "of", length(all.TTR), "\n") #keep track of what tree/loop# we're on
  sister.pairs <- NULL
  for (i in 1:length(unique(all.TTR[[tt]]$species1))) {
      current <- filter(all.TTR[[tt]], species1==all.matches[i,1])
      matches <- subset(current, current$species2==all.matches[i,2])
      sister.pairs <- rbind(sister.pairs, matches)
  }
  all.sister.pairs <- rbind(all.sister.pairs, sister.pairs)
}
saveRDS(all.sister.pairs, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.AllTree.Comparisons.RDS")

(ggplot(all.sister.pairs, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.1)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse()
  + ggtitle("Pygopodoid Geckos:
Trait Distance between Overlapping and Non-overlapping Taxa"))

short.d <- filter(skink, dist_tree <= 23)
skink.plot <- (ggplot(short.d, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.1)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse()
  + ggtitle("Sphenomorphine Skinks:
Trait Distance between Overlapping and Non-overlapping Taxa"))



agam <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Agamids.AllTree.Comparisons.RDS")
  agam$clade <- "Agamidae"
    agam.1 <- agam[1:48,]
pygo <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.AllTree.Comparisons.RDS")
  pygo$clade <- "Pygopodoidea"
    pygo.1 <- pygo[1:114,]
bird <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Meliphagoids.AllTree.Comparisons.RDS")
  bird$clade <- "Meliphagoidea"
    bird.1 <- bird[1:52,]
mars <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Marsupials.AllTree.Comparisons.RDS")
  mars$clade <- "Marsupials"
    mars.1 <- mars[1:108,]
skink <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Skinks.AllTree.Comparisons.RDS")
  skink$clade <- "Skinks"
      skink.1 <- skink[1:124,]
all.clades <- rbind(agam, pygo)
  all.clades <- rbind(all.clades, bird)
    all.clades <- rbind(all.clades, mars)
      all.clades <- rbind(all.clades, skink)
        saveRDS(all.clades, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/AllClades.AllTree.Comparisons.RDS")
    all.1 <- rbind(agam.1, bird.1)
      all.1 <- rbind(all.1, mars.1)
        all.1 <- rbind(all.1, skink.1)
          all.1 <- rbind(all.1, pygo.1)
multiplot(agam.plot, pygo.plot, bird.plot, mars.plot, skink.plot, ncol=1)


short.d <- filter(all.clades, dist_tree <= 23)
(ggplot(short.d, aes(x=dist_tree, y=dist_trait, color=range_overlap))
              #+ geom_point(alpha=0.1)
              + geom_smooth(aes(fill=clade))
              + scale_x_reverse()
              + ggtitle("All Clades:
Trait Distance between Overlapping and Non-overlapping Taxa")
              + theme_classic())



short.dff <- subset(dff, dff$dist_tree <= 20)

(ggplot(dff, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.25)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse()
  + ggtitle("Pygopodoids Trait Distance between Overlapping and Non-overlapping Taxa"))
(ggplot(fuller, aes(x=dist_trait, fill=range_overlap))
  + geom_density(alpha=0.5))


# if you want to plot a series of species within a given genus:
test <- distribution[which(distribution$Name_in_Tree %in% all.taxa[175:178]),]
sbbox <- make_bbox(lon = distribution$Longitude, lat = distribution$Latitude, f = .1)
sq_map <- get_map(location = sbbox, maptype = "terrain", source = "google")
(ggmap(sq_map) + geom_point(data = test, mapping = aes(x = Longitude, y = Latitude, color=Name_in_Tree)))


pygo <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Pygopodoids.Matrices.RDS")

test <- dff
test <- rbind(test, agam)



## Set your parameters relevant to your empirical parameter estimates
######################################################
# alpha
#alpha = 2 # either a static value
diff.alpha <- seq(from=0.5, to=5, by=0.1) # or set it as a vector of sampled values
diff.alpha <- sample(diff.alpha, size=100, replace=T) # or set it as a vector of sampled values
# pre/post shift sigma
preshift.sigma  = 1
postshift.sigma = 5
# and shift time
sim.shifts <- seq(from=1, to=20, by=1)
sim.shifts <- sample(sim.shifts, size=100, replace=T)
######################################################
bm.pars <- 0.1 # set the diffusion parameter of the BM process
ou.pars <- c(0.5, sample(diff.alpha, 1), 1) # set the diffusion parameter, the alpha, and the optimum

##########################################################################################################
## LOOP 1:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either the SRC ('Single-rate-constraint', BM-OU, 1 sig2) or the 
### TRC ('Two-rate-constraint', BM-OU, 2 sig2) model, then comparatively fit a set of standard models 
### (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a different model use the second loop, to simulate under a Brownian
### motion or Ornstein-Uhlenbeck process, using Diversitree (although you could also use Geiger).

### This loop will simulate data onto the tree with fossil tips, fit a series of models, then
### drop the extinct tips and associated data, refit the same models, and provide a summary for each
##########################################################################################################
for (z in 1:1) {
  traits.geiger <- list(); traits.ouwie <- list() # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- tree # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(input.trees), "\n") #keep track of what tree/loop# we're on
  
  for (i in 1:num.sims) {
    # Option A (comment out top when simulating different shift times, comment out bottom when simulating different alphas)
    m <- split.matrices <- split.vcv(phy, sim.shifts[z]) # divide the vcv matrix at a static time
    #m <- split.matrices <- split.vcv(phy, sim.shifts[i]) # or at differing times
    
    # Option B (adjust to change simulating model, comment out both to get the SRC model)
    #m[[1]] <- m[[1]] * preshift.sigma # adjust BM (old era) vcv according to a rate scalar (usually = 1)
    #m[[2]] <- m[[2]] * postshift.sigma # adjust OU (new era) vcv according to a rate scalar (faster or slower)
    
    # Option C (comment out top when simulating different times, comment out bottom when simulating different alphas)
    m[[2]] <- ouMatrix(m[[2]], alpha=diff.alpha[z]) # transform the second half of the vcv matrix according to your alpha
    #m[[2]] <- ouMatrix(m[[2]], alpha=alpha) # transform the second half of the vcv matrix according to your alpha
    
    m.rev.rel.rad <- m[[1]] + m[[2]] # combine the two matrices back together
    
    # OR, do it like the 'ecological release' model
    # m <- lapply(m, function(x) x*0.1)
    # m.rev.rel <- m[[1]] + m[[2]]
    
    # draw simulated data from a multivariate normal distribution, with appropriate root state (mean)
    traitz[[i]] <- setNames(rmvnorm(n=1, mean=rep(1, nrow(m.rev.rel.rad)), 
                                    sigma=m.rev.rel.rad), rownames(m.rev.rel.rad))
    t.ouwie <- NULL
    t.ouwie <- as.data.frame(names(traitz[[i]]))
    t.ouwie[,2] <- as.data.frame(t(traitz[[i]]))
    #traitz.ouwie[[i]] <- t.ouwie
    
    t.geiger <- t.ouwie; names <- t.ouwie[,1]
    rownames(t.geiger) <- names; t.geiger[,1] <- NULL
    #traitz.geiger[[i]] <- t.geiger
    
    traits.geiger[[i]] <- t.geiger
    traits.ouwie[[i]]  <- t.ouwie
  }
  sim.traits.geiger[[z]] <- traits.geiger
  sim.traits.ouwie[[z]] <- traits.ouwie
  
  
  save.sim.traits[[z]] <- as.data.frame(sim.traits.geiger[[z]]); 
  save.sim.traits[[z]][,"tree.num"] <- z
  save.sim.traits[[z]][,"gen.model"] <- "SRC"
  #save(save.sim.traits, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/Simulated.Traits.PlioPleistocene.SRC.RData")
}
datas <- as.data.frame(traits.ouwie)
rpanda.traits <- datas[,2]; names(rpanda.traits) <- datas[,1]
traits <- rpanda.traits


##########################################################################################################
## LOOP 2:
##########################################################################################################
### this loop will simulate data onto trees you've created with extinct tips (sim.fossil.trees)
### under either Brownian Motion (BM) or Ornstein-Uhlenbeck (OU) models then comparatively fit a set 
### of standard models (BM,EB, OU, SRC, TRC) to the simulated data.
### If you want to simulate data under a mode-variable model, use the first loop.
##########################################################################################################
bm.pars <- 0.1 # set the diffusion parameter of the BM process
ou.pars <- c(0.5, sample(diff.alpha, 1), 1) # set the diffusion parameter, the alpha, and the optimum
######################################################
## Set empty object to hold your results
sim.traits.geiger <- list(); sim.traits.ouwie <- list() # make trait lists for all trees in geiger and ouwie data format
save.sim.traits <- NULL
sim.traits <- NULL
extant.data <- list()
# make sure to the outputs depending on if you're simulating different times, or alphas
# this means changing the 'm' and the 'm[[2]]' objects below!
num.sims <- 1 # designate the number of simulated data sets you like to create per tree
### if trees are sensitive to shift dates: (this is for the Plio-Pleistocene trees)
######################################################

for (z in 1:1) {
  traits.geiger <- list(); traits.ouwie <- NULL # make intermediate trait lists for each tree in geiger and ouwie data format
  phy <- input.trees[[z]] # designating the target tree
  traitz <- list(); #traitz.ouwie <- list(); traitz.geiger <- list() # make intermediary data lists
  #traitz.geiger <- NULL; traitz.ouwie <- NULL
  cat("iteration", z, "of", length(input.trees), "\n") #keep track of what tree/loop# we're on
  
  for (i in 1:num.sims) {
    simulated.traits <- NULL
    simulated.traits <- as.data.frame(sim.character(phy, model="bm", bm.pars))
    #simulated.traits <- as.data.frame(sim.character(phy, model="ou", ou.pars))
    traits.geiger[[i]] <- simulated.traits
  }
  
  sim.traits.geiger[[z]] <- traits.geiger
  sim.traits.ouwie[[z]][,1] <- rownames(simulated.traits)
  #sim.traits.ouwie[[z]] <- traits.ouwie
  
  save.sim.traits[[z]] <- as.data.frame(sim.traits.geiger[[z]]); 
  save.sim.traits[[z]][,"tree.num"] <- z
  save.sim.traits[[z]][,"gen.model"] <- "BM" # OU, lowBM, hiBM, SRC
  saveRDS(save.sim.traits, file="/Users/Ian/Google.Drive/R.Analyses/MioceneAustralia_TreeTraitRange/Simulated.Pygopodoidea.BM.Traits.RDS")
  
  datas <- as.data.frame(save.sim.traits[[z]])
  sim.traits <- datas[,1]; 
  names(sim.traits) <- rownames(datas)
  
}
extinct <- subset(stree.res, stree.res$tree.type == "fossil tree")
extant  <- subset(stree.res, stree.res$tree.type == "extant tree")
outz <- summarySE(extant, measurevar="weight$w", groupvars="model")
colnames(outz) <- c("model", "N", "w", "sd", "se", "ci")

write.csv(stree.res, file="/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Trait Simulations/PlioPleistocene.OU.results.csv")












########## The below section seems to be more work than it's worth, stick to Loess

#### Make a summary of the Disparity estimates within a sliding window (mean, lowerCI, upperCI)
sim.nums <- seq(0, max(short.d$dist_tree), .1) #set your sliding window by (start, end, distance of window movement)
overlap.estimates <- NULL; noverlap.estimates <- NULL
for (t in sim.nums) {
  #for (i in 0:23) { #use this if you want just a subset of the age
  time.min <- t
  time.max <- t+1 # adjust the window width if you're getting NaN values (no observations within a time period)
  
  over <- subset(short.d, short.d$range_overlap == "TRUE")
  nover<- subset(short.d, short.d$range_overlap == "FALSE")
  
  overlap.chunk <- subset(over, over$dist_tree>=time.min & over$dist_tree<=time.max)
  noverlap.chunk <- subset(nover, nover$dist_tree>=time.min & nover$dist_tree<=time.max)
  
  over.est <- as.data.frame(t(CI(overlap.chunk$dist_trait)))
  nover.est<- as.data.frame(t(CI(noverlap.chunk$dist_trait)))
  
  
  timing <- paste(time.min, "to", time.max); colnames(timing)
  over.output.frame <- NULL; nover.output.frame <- NULL
  over.output.frame <- as.data.frame(c(timing, over.est))
  nover.output.frame <- as.data.frame(c(timing, nover.est))
  colnames(over.output.frame) <- c("time.window", "overlap.upper95", "overlap.mean", "overlap.lower95")
  colnames(nover.output.frame) <- c("time.window", "noverlap.upper95", "noverlap.mean", "noverlap.lower95")
  #colnames(output.frame) <- NULL
  overlap.estimates <- rbind(overlap.estimates, over.output.frame)
  noverlap.estimates <- rbind(noverlap.estimates, nover.output.frame)
}
overlap.estimates[,"timing"] <- sim.nums # add a column with the times
#overlap.estimates <- emp.mean.estimates[-c(101),] # drop the last row which is superfluous
noverlap.estimates[,"timing"] <- sim.nums # add a column with the times
#sim.mean.estimates <- sim.mean.estimates[-c(101),] # drop the last row which is superfluous
total.estimates <- cbind(overlap.estimates, c(noverlap.estimates[,2:4]))

#### Now let's plot the trend in disparity from the posterior of trees as a confidence ribbon
rib <- (ggplot(data=total.estimates)
        + geom_ribbon(aes(x=timing, ymin=overlap.lower95, ymax=overlap.upper95, fill="overlap"))
        + geom_ribbon(aes(x=timing, ymin=noverlap.lower95, ymax=noverlap.upper95, fill="noverlap"))
        + theme_classic())

(ggplot(total.estimates) 
  + geom_line(aes(x=timing, y=overlap.mean))
  + geom_line(aes(x=timing, y=noverlap.mean)))








##### Below this line is the previous attempt to determine overlap using biomes, DISREGARD
##########################################################################################
range_overlap <- rep(NA, ncol(ucomb))
arrange <- data.frame(species1 = ucomb[1,], species2 = ucomb[2,], range_overlap, stringsAsFactors=F)
for (t in 1:length(arrange[,1])) {
  taxon1 <- subset(range, rownames(range)==arrange[t,1])
  t1 <- c(as.character(taxon1[[1]]), 
          as.character(taxon1[[2]]), 
          as.character(taxon1[[3]]), 
          as.character(taxon1[[4]]),
          as.character(taxon1[[5]]))
  taxon2 <- subset(range, rownames(range)==arrange[t,2])
  t2 <- c(as.character(taxon2[[1]]), 
          as.character(taxon2[[2]]), 
          as.character(taxon2[[3]]), 
          as.character(taxon2[[4]]),
          as.character(taxon2[[5]]))
  
  overlap <- intersect(t1, t2)
  overlap <- overlap[!is.na(overlap)]
  if (length(overlap > 0)) {
    arrange[t,3] <- "yes"
  } else {
    arrange[t,3] <- "no"
  }
}

dff <- cbind(dff, arrange$range_overlap)
colnames(dff) <- c("species1", "species2", "dist_trait", "dist_tree", "range_overlap")
short.dff <- subset(dff, dff$dist_tree <= 20)
saveRDS(dff, file="/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Pygopodoidea.Tree.Trait.Range.Matrices.RDS")


(ggplot(dff, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point()
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse())

(ggplot(short.dff, aes(x=dist_trait, fill=range_overlap))
  + geom_density(alpha=0.5))


saveRDS(all.taxa)
all.taxa <- readRDS(file="/Users/Ian/Google.Drive/R.Analyses/BayesTraits/Agamids.Tree.Trait.Range.Matrices.RDS")
all.taxa <- rbind.data.frame(all.taxa, dff)
shorty <- subset(all.taxa, all.taxa$dist_tree <= 20)
(ggplot(all.taxa, aes(x=dist_tree, y=dist_trait, color=range_overlap))
  + geom_point(alpha=0.5)
  + geom_smooth(aes(fill=range_overlap))
  + scale_x_reverse())


distributions <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Distribution_Data/Agamid_PointData/Agamids_PointData.csv", header=T)
distribution <- read.csv("/Users/Ian/Downloads/Clean_Agamids/Clean_Agamids.csv", header=T)
test <- unique(distribution$Species...matched)
tests <- unique(distributions$Genus_species)
setdiff(tests, test)
