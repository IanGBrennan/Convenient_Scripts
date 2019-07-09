# Combine information on missing data and sample names/phylogeny

# this is the sample info pulled from the tip names
t1 <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/Egernia_SampleInfo.csv", header=T)
# this is a data frame of missing data pulled from IQtree and plopped into a csv file
t2 <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/Egernia_TaxonSampleInfo.csv", header=T)

library(dplyr)

t3 <- left_join(t1, t2)
head(t3)
write.csv(t3, "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/Egernia_SampleInfo_ALL.csv")
t4 <- sort(t1, "Gap.Ambiguity", decreasing=T)

t5 <- arrange(t1, )
t5 <- arrange(t1, desc(percent_gap_ambiguity))
head(t5)


missing <- data.frame(percent_missing = t5$percent_gap_ambiguity*100)
rownames(missing) <- t5$tip_label
missing$other <- 100 - missing$percent_missing
missed <- missing[match(egtree$tip.label, rownames(missing)),] # probably not necessary, but maybe useful

egtree <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/RAxML/RAxML_bipartitions_T545.tre")

# plot the data against the tree, but make sure you haven't done some thing like ladderize(egtree)
plotTree.barplot(egtree, missed, add=TRUE, 
                 args.barplot=list(xlab="Percent Missing Data", col=c("darkblue", "lightBlue"), border=F))



