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

##################################################################################
## Here's another way which may be simpler and more valuable
##################################################################################
# read in your alignment
align.file <- read.dna("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/mtGenomes/Egernia_RefAligned_Assemblies_1_2_3_RENAMED.fasta", format="fasta")

# make a function to check how much missing data there is
checkMissing <- function(alignment, taxon.file, count.gaps.as.missing=TRUE, print.df=TRUE, missing.threshold=NULL){
  library(ape); library(seqinr)

  # I'm going to use the quick and dirty approach of getting a percentage of missing data
  # to identify taxa we might want to remove
  # But it might be worth it to think about using a base composition or content method
  # like 'baseContent' in the 'spiderDev' package to do this as well
  # the only worry is that it'll consistently flag outgroups
  
  #full.align <- read.dna(paste0(sub.dir, paste0(alignment, ".fasta")), format="fasta")
  full.align <- alignment
  if(count.gaps.as.missing==TRUE) {missing.bases <- c(2, 240, 4)} else missing.bases <- c(2,240)
  missing.sum <- apply(full.align, MARGIN = 1, FUN = function(x) length(which(as.numeric(x) %in% missing.bases)))
  missing.percent <- missing.sum/ncol(full.align)
  missing.df <- data.frame(new_tip_label=names(missing.percent), mitoGenome_percent_incomplete=missing.percent); rownames(missing.df)<-NULL
  if(print.df==TRUE){print(missing.df)}
  
  new.data <- dplyr::left_join(taxon.file, missing.df)
  return(new.data)
  
  if(!is.null(missing.threshold)){
    bad.names <- names(which(missing.percent >= missing.threshold))
    print(paste("dropping", bad.names))
    good.names <- setdiff(rownames(full.align), bad.names)
    new.align <- full.align[which(rownames(full.align) %in% good.names), ]
    write.FASTA(new.align, paste0(sub.dir, paste0(alignment,"_Reduced"),".fasta"))
    cat(c("your reduced alignment is now called:\n", 
          paste0(sub.dir, paste0(alignment,"_Reduced"),".fasta")))
  } else cat(c("your alignment has not been changed"))
}

# spit out the sample info file with the 
new.info <- checkMissing(alignment=align.file, taxon.file=t1, count.gaps.as.missing=TRUE, missing.threshold=NULL, print.df=FALSE)
#write.csv(new.info, "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/Egernia_SampleInfo_2.csv")

incomplete <- new.info[complete.cases(new.info$mitoGenome_percent_incomplete),]
missing <- data.frame(percent_missing = incomplete$mitoGenome_percent_incomplete*100)
rownames(missing) <- incomplete$new_tip_label
missing$other <- 100 - missing$percent_missing
egtree <- read.tree("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T545_Egernia/mtGenomes/Egernia_RefAligned_Assemblies_1_2_3_RENAMED.tre")

missed <- missing[match(egtree$tip.label, rownames(missing)),] # probably not necessary, but maybe useful
#missed["Scincella_reference",] <- c(0,100)

# plot the data against the tree, but make sure you haven't done some thing like ladderize(egtree)
plotTree.barplot(egtree, missed, 
                 args.barplot=list(xlab="Percent Missing Data", col=c("darkblue", "lightBlue"), border=F))
