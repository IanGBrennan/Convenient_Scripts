go.tree <- read.tree("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Treefiles/ASTRAL_Komodo_NEWICK.tre.cf.tree")
    #go.tree <- reroot(go.tree, node=which(go.tree$tip.label == "Sphenodon_punctata"))
bval <- read.csv("/Users/Ian/Desktop/GenomeStripper/Komodo_assembly/Existing_Alignments/combined_alignments/Trimmed_Alignments/Treefiles/ASTRAL_Komodo_gCFs.csv", header=T)
  bval[1,"gCF"] <- 100; bval <- rbind(bval, c(min(bval$child)-1, 1, NA,NA,NA,NA,NA))
    #bval$child <- bval$child + 1; rbind(bval, c(min(bval$child)-1, 1, NA,NA,NA,NA,NA)); rbind(bval, c(min(bval$child)-2, 1, NA,NA,NA,NA,NA))
        

plotBranchbyTrait(go.tree, bval$gCF, mode="edges")

testo <- rtree(7); plot(testo); edgelabels(fr="c", cex=1); nodelabels(fr="c", cex=1)
edge.cols <- c(rep("green", 6), rep("blue",6))
plot(testo, edge.color=edge.cols)

testa <- "(A,(B,((C,D),(E,(F,G)))));"
cat(testa, file = "testa.tre", sep = "\n")
testa <- read.tree("testa.tre")
plot(testa); edgelabels(fr="c", cex=1); nodelabels(fr="c", cex=1); tiplabels(fr="c", cex=1)
testa$edge

edge.info <- data.frame("ID"=c(1,2,3,4,5,6,7,8,9,10,11,12), 
                        "GCF"=c(100,90,50,60,90,80,60,60,80,20,40,10),
                        "GCF_col"=c("black","black","red","red","black",
                                    "black","red","red","black","red","red","red"))
plot(testa, edge.color=edge.info$GCF_col, edge.width=2); edgelabels(fr="c", cex=1)


#clz <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- clz(max(edge.info$GCF) - min(edge.info$GCF))
clz <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- clz(100)

#edge.info <- edge.info[order(edge.info$GCF),] # don't need to reorder them, it's annoying
edge.cols <- point.colors[edge.info$GCF]
edge.info$new_col <- edge.cols
plot(testa, edge.color=edge.info$new_col, edge.width=5); edgelabels(fr="c", cex=1.5, col="black", bg="white")


edge.info <- data.frame("child"=c(9,10,11,12,13), 
                        "gCF"=c(10,90,40,20,60))
                        #"gCF_col"=c("black","black","red","red","black"))


plot(go.tree, cex=0.3); nodelabels(go.tree$node.label,node=2:go.tree$Nnode+Ntip(go.tree),
                                adj=c(1,-0.2),frame="none")

all.edges <- data.frame(parent=testa$edge[,1], child=testa$edge[,2])
tip.edges <- subset(all.edges, all.edges$child <= Ntip(testa))
    tip.edges$gCF <- 100
int.edges <- subset(all.edges, all.edges$child > Ntip(testa))
    int.edges <- left_join(int.edges, edge.info)
combo.edges <- rbind(tip.edges, int.edges)
combo.edges <- combo.edges[order(match(combo.edges$child, all.edges$child)),]

clz <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- clz(100)
edge.cols <- point.colors[combo.edges$gCF]
combo.edges$gCF_col <- edge.cols

plot(testa, edge.color=combo.edges$gCF_col, edge.width=5); edgelabels(fr="c", cex=1.5, col="black", bg="white"); nodelabels(fr="c", cex=1.5, col="white", bg="black")


all.edges <- data.frame(parent=go.tree$edge[,1], child=go.tree$edge[,2])
tip.edges <- subset(all.edges, all.edges$child <= Ntip(go.tree))
    tip.edges$gCF <- 100
int.edges <- subset(all.edges, all.edges$child > Ntip(go.tree))
    int.edges <- left_join(int.edges, bval)
        int.edges <- int.edges[,c("parent", "child", "gCF")]
combo.edges <- rbind(tip.edges, int.edges)
combo.edges <- combo.edges[order(match(combo.edges$child, all.edges$child)),]     

clz <- colorRampPalette(brewer.pal(6,"RdYlBu")); point.colors <- clz(100)
edge.cols <- point.colors[combo.edges$gCF]
combo.edges$gCF_col <- edge.cols

plot(go.tree, edge.color=combo.edges$gCF_col, edge.width=5); 
edgelabels(fr="c", cex=1.5, col="black", bg="white"); nodelabels(fr="c", cex=1.5, col="white", bg="black")



#################################################
# read the data in the better way!
# I should build this into the function directly
#################################################
input.tree <- read.tree("~/Desktop/GenomeStripper/Elapids/Existing_Alignments/combined_alignments/gCF.cf.tree")
d = read.delim("~/Desktop/GenomeStripper/Elapids/Existing_Alignments/combined_alignments/gCF.cf.stat", header=T, comment.char="#")
names(d)[1] = "child"
#names(d)[11] = "branchlength"

d[1,"gCF"] <- 100; 
d <- rbind(d, c(min(d$child)-1, 100, NA,NA,NA,NA,NA))
color.by.CF(curr.tree, d)   

source("~/Google.Drive/R.Analyses/Convenient Scripts/color.by.CF.R")
color.by.CF(input.tree, "~/Desktop/GenomeStripper/Elapids/Existing_Alignments/combined_alignments/gCF.cf.stat", 
            terminal.color=NULL, legend=T)

plot(rep(1,100),col=colorRampPalette(brewer.pal(6, "RdYlBu"))(100),pch=15,cex=10, axes=F);


#plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))
layout(mat = matrix(c(1, 2, 1, 0), 
                    nrow = 2, 
                    ncol = 2),
       heights = c(4, 1),    # Heights of the two rows
       widths = c(5, 0))     # Widths of the two columns




color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

color.bar(colorRampPalette(brewer.pal(6, "RdYlBu"))(100))

# now let's try it out
testo <- color.by.CF(input.tree, d)
# now combine the original and new data frames, and drop incomplete fields
besto <- left_join(testo, d); besto <- besto[complete.cases(besto),]
# and plot gCFs as a function of branch length
library(ggplot2)
ggplot(besto, aes(x=Length, y=gCF, label=child)) +
  geom_point(shape=21, fill=besto$gCF_col, size=10) + theme_classic() + geom_smooth(color="black",linetype="dashed", level=0)

# we could fit a linear regression to the gCF/branch length data
source("~/Google.Drive/R.Analyses/Convenient Scripts/ggplotRegression.R")
ggplotRegression(lm(besto$gCF ~ besto$Length), alpha=0.75, size=5, color=besto$gCF_col)

# or we could just plot it with a best fit Loess line (which probably makes more sense)
ggplot(besto, aes(x=Length, y=gCF, label=child)) +
  geom_point(shape=21, fill=besto$gCF_col, size=10) + theme_classic() + geom_smooth(color="black",linetype="dashed", level=0)

