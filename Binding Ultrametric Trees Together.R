library(ape)
library(phytools)
library(TreePar)
library(phangorn)

basis <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Zheng.Wiens.Squamates.tre")
tips <- basis$tip.label
write.csv(tips, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/ZW.Squamate.Tips.csv")

drops <- c("Acrochordus_granulatus",
           "Anilius_scytale",
           "Anolis_carolinensis",
           "Bothrops_asper",
           "Casarea_dussumieri",
           "Iguana_iguana",
           "Lampropeltis_getula",
           "Leptotyphlops_humilis",
           "Liotyphlops_albirostris",
           "Micrurus_fulvius",
           "Thamnophis_marcianus",
           "Trachyboa_boulengeri",
           "Trachyboa_gularis",
           "Tropidophis_curtus_curtus",
           "Tropidophis_feicki",
           "Tropidophis_greenwayi",
           "Tropidophis_haetianus_haetianus",
           "Tropidophis_melanurus_melanurus",
           "Tropidophis_pardalis",
           "Tropidophis_taczanowskyi",
           "Tropidophis_wrighti",
           "Typhlops_jamaicensis",
           "Xenophidion_schaeferi")
           
pythons <- lapply(henos, drop.tip, tip=drops)
pythons[[1]]
class(pythons) <- "multiPhylo"
write.nexus(pythons, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Henophidians.100.trees")


#######################################################################
# Binding trees together to make a single ultrametric tree
#######################################################################

# we need to start from the shallowest trees, and work our way backward.
# trying to imprint the age afterwards (with 'chronos') will adjust our branching times!

# what I have:
# 100 PB Gecko Trees
# 100 PB Skink Trees
# 100 PB Blindsnake Trees
# 100 PB Henophidian Trees (Boidae, Uro/Cylindro, Xeno/Loxo/Pythonidae)
# 100 PB Caenophidian Trees (Elapidae + Colubridae)
# 100 PB Aniliidae+Tropidophiidae Trees

# these trees are finalized
elapids <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Elapids.100.trees")
geckos <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Pygopodoidea.100.trees")
anilios <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Anilios.Blindies.100.trees")
viperids <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Viperids.100.trees")
colubrids <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Colubrids.100.trees")
henos <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Henophidians.100.trees")
anil.trop <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Aniliidae.Tropidophiidae.100.trees")
gerro.cord <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Gerro.Cordylidae.tre")

# these trees could be expanded, by extracting larger clades from the Z/W tree
## requires collecting more data from the literature
skinks <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Skinks.100.trees")
anguids <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Anguidae.tre")
gymno.amphis <- read.tree("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Gymno.Amphis.tre")
acontines <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Acontines.tre")
lygosomines <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/Lygosomines.tre")



########
# Elapid and Colubrid trees are non-Ultrametric due to rounding errors, here's the fix:
# instructions: http://blog.phytools.org/2016/08/fixing-ultrametric-tree-whose-edges-are.html
########
fixed.list <- list()
for (i in 1:length(viperids)) {
  fixed <- nnls.tree(cophenetic(viperids[[i]]), viperids[[i]], rooted=T)
  fixed.list[[i]] <- fixed
}
class(fixed.list) <- "multiPhylo"
is.ultrametric(fixed.list) # check to make sure it worked
write.nexus(fixed.list, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/PB.Viperids.100.trees")

########
# Step 1a: Elapids and Colubrids
########
# add roots to the elapid and colubrid trees (diverged ~50 Million years ago)
ages <- rnorm(100, mean=55, sd=1) # pull the age of their MRCA from a normal dist. (provide uncertainty)

elapid.list <- list()
for (i in 1:length(elapids)) {
  tree.max <- max(nodeHeights(elapids[[i]]))
  root.length <- (ages[i]-tree.max) 
  elapid.list[[i]] <- addroot(elapids[[i]], root.length)
}
class(elapid.list) <- "multiPhylo" # this is important

colubrid.list <- list()
for (i in 1:length(colubrids)) {
  tree.max <- max(nodeHeights(colubrids[[i]]))
  root.length <- (ages[i]-tree.max) 
  colubrid.list[[i]] <- addroot(colubrids[[i]], root.length)
}
class(colubrid.list) <- "multiPhylo" # this is important

# loop through and add all the colubrids and elapids together
ce.list <- list()
for (i in 1:length(elapid.list)) {
  tree1 <- elapid.list[[i]]
  tree2 <- colubrid.list[[i]]
  ce.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(ce.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(ce.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(ce.list) # check as we go along to make sure it's staying ultrametric


########
# Step 1b: Colubroids (Elapids/Colubrids) and Viperids
########
# add roots to the colubroid and viperid trees (diverged ~60 Million years ago)
ages <- rnorm(100, mean=60, sd=2)

# first deal with the colubroid tree we made above
ce.root <- list()
for (i in 1:length(ce.list)) {
  tree.max <- max(nodeHeights(ce.list[[i]]))
  root.length <- (ages[i]-tree.max)
  ce.root[[i]] <- addroot(ce.list[[i]], root.length)
}
class(ce.root) <- "multiPhylo"

# now deal with the henophidians
vip.list <- list()
for (i in 1:length(viperids)) {
  tree.max <- max(nodeHeights(viperids[[i]]))
  root.length <- (ages[i]-tree.max) 
  vip.list[[i]] <- addroot(viperids[[i]], root.length)
}
class(vip.list) <- "multiPhylo" # this is important

# loop through and add all the colubroids and viperids together
caeno.list <- list()
for (i in 1:length(vip.list)) {
  tree1 <- vip.list[[i]]
  tree2 <- ce.root[[i]]
  caeno.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(caeno.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(caeno.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(caeno.list) # check as we go along to make sure it's staying ultrametric


########
# Step 2: Caenophidians (Elapids/Colubrids/Viperids) and Henophidians (Boas, Pythons, etc)
########
# add roots to the elapid and python trees (diverged ~85 Million years ago)
ages <- rnorm(100, mean=92, sd=1) # pull the age of their MRCA from a normal dist. (provide uncertainty)

# first deal with the caenophidian tree we made above
caeno.root <- list()
for (i in 1:length(caeno.list)) {
  tree.max <- max(nodeHeights(caeno.list[[i]]))
  root.length <- (ages[i]-tree.max)
  caeno.root[[i]] <- addroot(caeno.list[[i]], root.length)
}
class(caeno.root) <- "multiPhylo"

# now deal with the henophidians
heno.list <- list()
for (i in 1:length(henos)) {
  tree.max <- max(nodeHeights(henos[[i]]))
  root.length <- (ages[i]-tree.max) 
  heno.list[[i]] <- addroot(henos[[i]], root.length)
}
class(heno.list) <- "multiPhylo" # this is important

# loop through and add all the caenophidians and henophidians together (caeno+heno)
ch.list <- list()
for (i in 1:length(heno.list)) {
  tree1 <- caeno.root[[i]]
  tree2 <- heno.list[[i]]
  ch.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(ch.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(ch.list[[1]], show.tip.label=F) # have a look to make sure it worked


########
# Step 3: Add Anilius/Tropidophiids to the Caeno/Henophidian Tree
########
# add roots to the anil.trop and caeno/heno trees (diverged ~95 Million years ago)
ages <- rnorm(100, mean=95, sd=1)

# first deal with the caeno/heno tree we made above
ch.root <- list()
for (i in 1:length(ch.list)) {
  tree.max <- max(nodeHeights(ch.list[[i]]))
  root.length <- (ages[i]-tree.max)
  ch.root[[i]] <- addroot(ch.list[[i]], root.length)
}
class(ch.root) <- "multiPhylo"

# now add roots on the aniliids/tropidophiids
anil.trop.list <- list()
for (i in 1:length(anil.trop)) {
  tree.max <- max(nodeHeights(anil.trop[[i]]))
  root.length <- (ages[i]-tree.max) 
  anil.trop.list[[i]] <- addroot(anil.trop[[i]], root.length)
}
class(anil.trop.list) <- "multiPhylo" # this is important

# loop through and add all the caeno/henophidians and aniliids/tropidophiids together (macrostomata)
macro.list <- list()
for (i in 1:length(anil.trop.list)) {
  tree1 <- ch.root[[i]]
  tree2 <- anil.trop.list[[i]]
  macro.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(macro.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(macro.list[[1]], show.tip.label=F) # have a look to make sure it worked


########
# Step 4: Anilios
########
# add roots to the macrostomatan and anilios trees (diverged ~125 Million years ago)
ages <- rnorm(100, mean=125, sd=4)

# first deal with the macrostomatan tree we made above
macro.root <- list()
for (i in 1:length(macro.list)) {
  tree.max <- max(nodeHeights(macro.list[[i]]))
  root.length <- (ages[i]-tree.max)
  macro.root[[i]] <- addroot(macro.list[[i]], root.length)
}
class(macro.root) <- "multiPhylo"

# next add roots to the blindsnakes
anilios.list <- list()
for (i in 1:length(anilios)) {
  tree.max <- max(nodeHeights(anilios[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  anilios.list[[i]] <- addroot(anilios[[i]], root.length)
}
class(anilios.list) <- "multiPhylo" # this is important


# loop through and add all the blindsnakes to the macrostomatans (ophidians)
snake.list <- list()
for (i in 1:length(anilios.list)) {
  tree1 <- macro.root[[i]]
  tree2 <- anilios.list[[i]]
  snake.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(snake.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(snake.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(snake.list) # check to make sure we're doing this right

########
# Step 5: Anguids
########
# add roots to the ophidian tree (diverged ~160 Million years ago)
ages <- rnorm(100, mean=162, sd=1)

# first deal with the macrostomatan tree we made above
ophid <- list()
for (i in 1:length(macro.list)) {
  tree.max <- max(nodeHeights(snake.list[[i]]))
  root.length <- (ages[i]-tree.max)
  ophid[[i]] <- addroot(snake.list[[i]], root.length)
}
class(ophid) <- "multiPhylo"

# next add roots to the anguids
anguid.list <- list()
for (i in 1:length(anguids)) {
  tree.max <- max(nodeHeights(anguids[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  anguid.list[[i]] <- addroot(anguids[[i]], root.length)
}
class(anguid.list) <- "multiPhylo" # this is important


# loop through and add all the anguids to the ophidians
ang.ophid.list <- list()
for (i in 1:length(anguid.list)) {
  tree1 <- ophid[[i]]
  tree2 <- anguid.list[[i]]
  ang.ophid.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(ang.ophid.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(ang.ophid.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(ang.ophid.list) # check to make sure we're doing this right


########
# Step 6: Sphenomorphines + Lygosomines
########
# add roots to the lygosomine and sphenomorphine trees (diverged ~80 Million years ago)
ages <- rnorm(100, mean=80, sd=1)

spheno.list <- list()
for (i in 1:length(skinks)) {
  tree.max <- max(nodeHeights(skinks[[i]]))
  root.length <- (ages[i]-tree.max) 
  spheno.list[[i]] <- addroot(skinks[[i]], root.length)
}
class(spheno.list) <- "multiPhylo" # this is important

lygo.list <- list()
for (i in 1:length(lygosomines)) {
  tree.max <- max(nodeHeights(lygosomines[[i]]))
  root.length <- (ages[i]-tree.max) 
  lygo.list[[i]] <- addroot(lygosomines[[i]], root.length)
}
class(lygo.list) <- "multiPhylo" # this is important

# loop through and add all the colubrids and elapids together
sl.list <- list()
for (i in 1:length(lygo.list)) {
  tree1 <- lygo.list[[i]]
  tree2 <- spheno.list[[i]]
  sl.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(sl.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(sl.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(sl.list) # check as we go along to make sure it's staying ultrametric


########
# Step 7: Acontine skinks + Other Skinks
########
# add roots to the skink tree (diverged ~90 Million years ago)
ages <- rnorm(100, mean=90, sd=1.5)

# first deal with the macrostomatan tree we made above
oid <- list()
for (i in 1:length(sl.list)) {
  tree.max <- max(nodeHeights(sl.list[[i]]))
  root.length <- (ages[i]-tree.max)
  oid[[i]] <- addroot(sl.list[[i]], root.length)
}
class(oid) <- "multiPhylo"

# next add roots to the acontines
acont.list <- list()
for (i in 1:length(acontines)) {
  tree.max <- max(nodeHeights(acontines[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  acont.list[[i]] <- addroot(acontines[[i]], root.length)
}
class(acont.list) <- "multiPhylo" # this is important


# loop through and add all the anguids to the ophidians
skink.list <- list()
for (i in 1:length(acont.list)) {
  tree1 <- oid[[i]]
  tree2 <- acont.list[[i]]
  skink.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(skink.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(skink.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(skink.list) # check to make sure we're doing this right


########
# Step 8: Skinks+Cordyloids
########
# add roots to the skink+cordyloid tree (diverged ~150 Million years ago)
ages <- rnorm(100, mean=150, sd=1)

# first deal with the skinks tree we made above
scincoid <- list()
for (i in 1:length(skink.list)) {
  tree.max <- max(nodeHeights(skink.list[[i]]))
  root.length <- (ages[i]-tree.max)
  scincoid[[i]] <- addroot(skink.list[[i]], root.length)
}
class(scincoid) <- "multiPhylo"

# next add roots to the cordyloids
cord.list <- list()
for (i in 1:length(gerro.cord)) {
  tree.max <- max(nodeHeights(gerro.cord[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  cord.list[[i]] <- addroot(gerro.cord[[i]], root.length)
}
class(cord.list) <- "multiPhylo" # this is important


# loop through and add all the anguids to the ophidians
cordski.list <- list()
for (i in 1:length(cord.list)) {
  tree1 <- cord.list[[i]]
  tree2 <- scincoid[[i]]
  cordski.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(cordski.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(cordski.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(cordski.list) # check to make sure we're doing this right

########
# Step 9: Gymnopthalmids/Amphisbaenians + Anguid/Snakes
########
# add roots to this tree (diverged ~150 Million years ago)
ages <- rnorm(100, mean=170, sd=2)

# first deal with the anguid/snakes tree we made above
ao <- list()
for (i in 1:length(ang.ophid.list)) {
  tree.max <- max(nodeHeights(ang.ophid.list[[i]]))
  root.length <- (ages[i]-tree.max)
  ao[[i]] <- addroot(ang.ophid.list[[i]], root.length)
}
class(ao) <- "multiPhylo"

# next add roots to the gymno/amphis
ga.list <- list()
for (i in 1:length(gymno.amphis)) {
  tree.max <- max(nodeHeights(gymno.amphis[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  ga.list[[i]] <- addroot(gymno.amphis[[i]], root.length)
}
class(ga.list) <- "multiPhylo" # this is important


# loop through and add all the gymno/amphis to the anguid/snakes
toid.list <- list()
for (i in 1:length(ga.list)) {
  tree1 <- ga.list[[i]]
  tree2 <- ao[[i]]
  toid.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(toid.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(toid.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(toid.list) # check to make sure we're doing this right


########
# Step 10: Skinks to Higher Squamates
########
# add roots to this tree (diverged ~150 Million years ago)
ages <- rnorm(100, mean=180, sd=2)

# first deal with the gymno/amphis/anguid/snakes tree we made above
gaas <- list()
for (i in 1:length(toid.list)) {
  tree.max <- max(nodeHeights(toid.list[[i]]))
  root.length <- (ages[i]-tree.max)
  gaas[[i]] <- addroot(toid.list[[i]], root.length)
}
class(gaas) <- "multiPhylo"

# next add roots to the cordylid/skinks
skicord.list <- list()
for (i in 1:length(cordski.list)) {
  tree.max <- max(nodeHeights(cordski.list[[i]]))
  root.length <- (ages[i]-tree.max) # this is the important step to change!
  skicord.list[[i]] <- addroot(cordski.list[[i]], root.length)
}
class(skicord.list) <- "multiPhylo" # this is important


# loop through and add all the gymno/amphis to the anguid/snakes
nogek.list <- list()
for (i in 1:length(skicord.list)) {
  tree1 <- gaas[[i]]
  tree2 <- skicord.list[[i]]
  nogek.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(nogek.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(nogek.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.ultrametric(nogek.list) # check to make sure we're doing this right

########
# Step 11: Add on the geckos (pygopodoids)
########
ages <- rnorm(100, mean=195, sd=3)

# first add roots to the lizard (skank) trees (all that exclude the geckos)
skank.root <- list()
for (i in 1:length(nogek.list)) {
  tree.max <- max(nodeHeights(nogek.list[[i]]))
  root.length <- (ages[i]-tree.max)
  skank.root[[i]] <- addroot(nogek.list[[i]], root.length)
}
class(skank.root) <- "multiPhylo"


# and then add roots to the gecko tree
gecko.list <- list()
for (i in 1:length(geckos)) {
  tree.max <- max(nodeHeights(geckos[[i]]))
  root.length <- (ages[i]-tree.max)
  gecko.list[[i]] <- addroot(geckos[[i]], root.length)
}
class(gecko.list) <- "multiPhylo"


# loop through and add all the lizards to the geckos
all.list <- list()
for (i in 1:length(gecko.list)) {
  tree1 <- gecko.list[[i]]
  tree2 <- skank.root[[i]]
  all.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(all.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(all.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.binary(all.list); is.ultrametric(all.list) # make sure they're all binary and ultrametric

########
# Step 12: Add on the dibamus
########
ages <- rnorm(100, mean=205, sd=2)

# first add roots to the lizard (skank) trees (all that exclude the geckos)
all.root <- list()
for (i in 1:length(all.list)) {
  tree.max <- max(nodeHeights(all.list[[i]]))
  root.length <- (ages[i]-tree.max)
  all.root[[i]] <- addroot(all.list[[i]], root.length)
}
class(all.root) <- "multiPhylo"


# make a simple two taxon tree as surrogate for Dibamus
twotree <- pbtree(n=2, scale=15)
twotree$tip.label[2] <- "Dibamus_novaeguineae"
twotree$tip.label[1] <- "Dibamus_taylori"

# and then add roots to the "dibamus tree"
dibamus <- list()
for (i in 1:100) {
  tree.max <- max(nodeHeights(twotree))
  root.length <- (ages[i]-tree.max)
  dibamus[[i]] <- addroot(twotree, root.length)
}
class(dibamus) <- "multiPhylo"


# loop through and add all the lizards to the geckos
final.all.list <- list()
for (i in 1:length(dibamus)) {
  tree1 <- all.root[[i]]
  tree2 <- dibamus[[i]]
  final.all.list[[i]] <- bind.tree(tree1, tree2, where="root")
}
class(final.all.list) <- "multiPhylo" # designate it as a multiPhylo object
plot(final.all.list[[1]], show.tip.label=F) # have a look to make sure it worked
is.binary(final.all.list); is.ultrametric(final.all.list) # make sure they're all binary and ultrametric

final.all <- lapply(final.all.list, drop.tip, tip="t1")

write.nexus(final.all, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")


########
# Step 13: Add on another Dibamus, because I'm an idiot
########
trees <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL.UnTrimmed.Limbless.100.trees")
tip <- "Dibamus_taylori"
sister <- "Dibamus_novaeguineae"
tree.test <- bind.tip(trees[[1]],tip,where=which(trees[[1]]$tip.label==sister),
               position=0.1*trees[[1]]$edge.length[which(trees[[1]]$edge[,2]==
                                                     which(trees[[1]]$tip.label==sister))])
new.trees <- NULL
for (i in 1:length(trees)) {
  current.tree <- trees[[i]]
  added.tree <- bind.tip(current.tree,tip,where=which(current.tree$tip.label==sister),
                         position=0.1*current.tree$edge.length[which(current.tree$edge[,2]==
                                                                     which(current.tree$tip.label==sister))])
  new.trees[[i]] <- added.tree
}
class(new.trees) <- "multiPhylo" # designate it as a multiPhylo object
write.nexus(new.trees, file="/Users/Ian/Google.Drive/R.Analyses/Limb.Reduced.Data/FINAL+1.UnTrimmed.Limbless.100.trees")


