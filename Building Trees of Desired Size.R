library(phytools)
library(TreeSim)

test <- pbtree(n=100, b=1.5, d=1)
length(test$tip.label)
plot(test)


get.desired.trees <- function(ntrees, extant.tips, b.rate, d.rate) {
  out.trees <<- NULL
  num.trees = 1
  while (num.trees <= ntrees) {
    xtree <- pbtree(n=extant.tips, b.rate, d.rate)
    tip.num <- length(xtree$tip.label)
    if (tip.num < extant.tips) {
      num.trees = num.trees
    } else {
      out.trees[[num.trees]] <<- xtree
      num.trees <- num.trees+1
    }
  }
}

test <- get.desired.trees(ntrees=3, extant.tips=50, b.rate=1.5, d.rate=1)
