library(Biostrings)
library(ips)


introns <- readDNAStringSet("/Users/Ian/Google.Drive/R.Analyses/Chen.Data/Intron_with_pangolin/Laurasiatheria_intron_with_pangolin_3638genes_Nuc.fas")
introns <- readDNAMultipleAlignment("/Users/Ian/Google.Drive/R.Analyses/Chen.Data/Intron_with_pangolin/Laurasiatheria_intron_with_pangolin_3638genes_Nuc.fas")
names(introns) = c("cat.Felis", "shrew.Sorex", "pangolin.Manis", "cow.Bos", "dog.Canis", "microbat.Myotis",
                               "panda.Ailuropoda", "rhino.Ceratotherium", "megabat.Pteropus", "human.Homo", "pig.Sus", 
                               "weasel.Mustela", "dolphin.Tursiops", "vicuna.Vicugna", "elephant.Loxodonta", "horse.Equus",
                               "sheep.Ovis", "armadillo.Dasypus", "sloth.Choloepus", "tenrec.Echinops", "hyrax.Procavia",
                               "hedgehog.Erinaceus", "mouse.Mus")

test <- subseq(introns[1], start=1, end=10)
test <- subseq(introns, start=1, end=1)
base.freq <- alphabetFrequency(introns, as.prob = T) # using 'baseOnly' won't distinguish between 'N' and '-'
rownames(base.freq) = seq_names
plot(base.freq[,"N"])
axis(side=1,at=c(1:length(seq_name)),labels=seq_names, tick=NULL)

masktest <- introns
rowmask(masktest) <- IRanges(start=1, end=3)




masked.gaps.introns <- maskGaps(introns, min.fraction=0.13)
m.sdist <- stringDist(as(masked.gaps.introns, "DNAStringSet"), method="hamming")
m.clust <- hclust(m.sdist, method="single")
plot(m.clust)

origMAlign <- readDNAMultipleAlignment(filepath=system.file("extdata", "msx2_mRNA.aln", package="Biostrings", format="clustal"))


sdist <- stringDist(as(introns, "DNAStringSet"), method="hamming")
clust <- hclust(sdist, method="single")
plot(clust)

intron.df <- as.data.frame(introns)
cat.shrew <- subset(intron.df, rownames(intron.df)=="cat.Felis" | intron.df$names=="shrew.Sorex")







