clade <- read.nexus("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.tre")
datum <-   read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.MlogBL.csv", header=F, row.names = 1)

source("/Users/Ian/Google.Drive/R.Analyses/Convenient Scripts/New.Models.adapted.from.Slater.2013.R"); ## now source the function from your WD


fitContinuous(clade, datum, model="EB")
fitContinuous_paleo(clade, datum, model="TRC", shift.time=4)
fitContinuous_paleo(clade, dat, model="BM")
fitContinuous_paleo(clade, datum, model="Burst.Bind", shift.time=7)
fitContinuous_paleo(clade, datum, model="Burst.Bind.Time")
fitContinuous_paleo(clade, datum, model="EB")
fitContinuous_paleo(clade, datum, model="Burst.Burst", shift.time=5)



require(TreeSim);
require(geiger);
clade$edge.length

ebclade <- rescale(clade, "EB", a=-0.4)
ebclade$edge.length
ebclade.re <- rescale(ebclade, "depth", depth=max(nodeHeights(clade)))
plot(ebclade)

fitContinuous(ebclade, datum, model="EB")
fitContinuous(ebclade.re, datum, model="EB")


ebphy <- exponentialchangeTree(phy, a = -0.4); ## now transform this tree under a strong EB process
ebphy$edge.length; # take a look at these edge lengths
plot(ebphy); # and plot it

### Notes from Graham Slater
#### http://grokbase.com/t/r/r-sig-phylo/12b4s8bczv/fitcontinuous-early-burst-model :
# what you'll probably find is that the exponentialchangeTree has a lot of very very short edges 
# and from plotting it you'll see that the only long edges are those towards the root. This is 
# such a strong burst that it basically requires all the evolutionary change in our trait to occur
# in the two edges leading from the root. Such a burst is unlikely in real data. I find it helpful 
# to think about rate half lives. An EB parameter of -0.4 gives a rate half life of 
# log(2) / 0.4 = 1.732868 time units - that is it takes ~1.8 million years for the rate to halve 
# from its initial value. To give an idea of what this translates to, the average mammalian order 
# is ~ 50 my old, which would result in approximately 30 half lives elapsing. That, at least to me, 
# has an exceedingly low prior probability!

# So while I agree with Dave that some paleo trees will have short internal edges or terminals 
# that might mess up the EB model, I think what is more important is to recognize the effect 
# that the EB parameter has on your edge lengths and to chose an appropriate lower bound when 
# fitting the model to avoid this. Probably, for most paleo trees, an exponential change parameter 
# of -0.4 is never going to be reasonable.


