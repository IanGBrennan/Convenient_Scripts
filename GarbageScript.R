go.tree <- read.tree("~/Google.Drive/R.Analyses/Varanus_Project/Varanidae_STRICT_HKY_270_con.newick")
go.data <- read.csv("~/Google.Drive/R.Analyses/Varanus_Project/Varanus_AllSVL.csv", header=T)

go.data <- go.data[complete.cases(go.data),]

keeps <- intersect(go.tree$tip.label, go.data$Name_in_Tree)
go.tree <- drop.tip(go.tree, setdiff(go.tree$tip.label, keeps))
go.data <- filter(go.data, Name_in_Tree %in% go.tree$tip.label)

go.svl <- go.data$SVL; names(go.svl) <- go.data$Name_in_Tree
go.logsvl <- log(go.svl)

# Plot SVL by sample:
all.data <- read.csv("~/Google.Drive/R.Analyses/Varanus_Project/Goanna.Marsupial.DATA.csv", header=T)
all.data <- all.data[order(all.data$Body_Length),]
all.data$Name_in_Tree <- factor(all.data$Name_in_Tree, levels = all.data$Name_in_Tree)
#all.data$Body_Length <- log(all.data$Body_Length)
ggplot(all.data, aes(x=Name_in_Tree, y=Body_Length, fill=Group), colour=brewer.pal(3,"Paired")) +
  geom_col() + theme_classic() + coord_flip() + scale_x_discrete(limits = rev(all.data$Name_in_Tree)) + scale_fill_brewer(palette="Paired", "Group")

# Continuous trait map
obj1 <- contMap(go.tree, go.logsvl, plot=FALSE, outline=F); 
n<-length(obj1$cols); 
obj1$cols[1:n] <- rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(n));
plot(obj1,type="fan", legend=0.7*max(nodeHeights(obj1$tree)),
     fsize=c(0.7,0.9), lwd=5, border=F); axisPhylo(1, backward=T)

# Barplot of trait value
plotTree.barplot(go.tree, go.svl, args.barplot=list(beside=TRUE, border=F))


go.body <- go.data[,c("SVL", "Tail")]; rownames(go.body) <- go.data$Name_in_Tree
plotTree.barplot(go.tree, go.body, args.barplot=list(beside=T, border=F))
# plot the data against the tree, but make sure you haven't done some thing like ladderize(egtree)
plotTree.barplot(go.tree, go.body, args.barplot=list(xlab="Varanus Total Length", col=c("darkblue", "lightBlue"), border=F))
plotTree(go.tree, type="fan", part=0.5)


test.tree <- drop.tip(go.tree, tip=c("Lanthanotus.borneensis", "Necrosaurus.cayluxi"))
plotTree(test.tree, type="fan", part=0.5) 
test.data <- go.svl[2:length(go.svl)]

plotTree.wBars(test.tree, test.data, type="fan", part=0.5)





tl <- go.body$SVL + go.body$Tail; go.body$Total <- tl
go.body <- go.body[,c("SVL", "Total")]
  
go.body <- as.matrix(go.body)
Z <- go.body

cols<-c("darkgrey","white")
scale<-0.5*max(nodeHeights(tree))/max(Z)
plotTree.wBars(tree,Z[,ncol(Z)],type="fan",scale=scale,col=cols[ncol(Z)], part=0.5)
obj<-get("last_plot.phylo", envir = .PlotPhyloEnv)
plotTree.wBars(tree,Z[,1],type="fan",scale=scale,add=TRUE,lims=obj$x.lim,
               col=cols[1], part=0.5)
obj<-legend(x="topleft",legend=colnames(Z),pch=22,
            pt.cex=1.5,pt.bg=cols,cex=0.9,bty="n",horiz=FALSE)
leg.length<-400 ## legend length
polygon(0.97*obj$rect$left+c(0,0,leg.length*scale,
                             leg.length*scale),-obj$text$y[2]+c(0,
                                                                diff(par()$usr[3:4])/100,diff(par()$usr[3:4])/100,0))
polygon(0.97*obj$rect$left+c(0,0,leg.length*scale/2,
                             leg.length*scale/2),-obj$text$y[2]+c(0,
                                                                  diff(par()$usr[3:4])/100,diff(par()$usr[3:4])/100,0),
        col="darkgrey")
for(i in 0:2){
  lines(rep(0.97*obj$rect$left+i/2*leg.length*scale,2),
        -obj$text$y[2]+c(0,-diff(par()$usr[3:4])/100))
  text(0.97*obj$rect$left+i/2*leg.length*scale,-obj$text$y[2],
       i*leg.length/2,cex=0.6,pos=1)
}

pdata <- go.body[,1]; names(pdata) <- rownames(go.body)
plotTree.wBars(tree, pdata, type="fan")


total<-setNames(round(runif(n=Ntip(tree),min=50,max=300)),tree$tip.label)
sampled<-round(runif(n=Ntip(tree),min=0.1,max=0.6)*total)
Z<-cbind(sampled,total)





tiger.size <- sample(800:1000, 50, replace=T)
tiger.lat <- sample(115:126, 50, replace=T)

suta.size <- sample(200:400, 50, replace=T)
suta.lat <- sample(138:147, 50, replace=T)

brachy.size1 <- sample(400:600, 25, replace=T)
brachy.size2 <- sample(500:700, 25, replace=T)
brachy.size3 <- sample(600:800, 25, replace=T)
brachy.lat1 <- sample(116:126, 25, replace=T)
brachy.lat2 <- sample(127:138, 25, replace=T)
brachy.lat3 <- sample(138:145, 25, replace=T)
brachy <- data.frame(SVL=c(brachy.size1, brachy.size2, brachy.size3), Latitude=c(brachy.lat1, brachy.lat2, brachy.lat3), Species=rep("Brachyurophis", 75))

ends <- data.frame(SVL=c(tiger.size, suta.size), Latitude=c(tiger.lat, suta.lat), Species=c(rep("Notechis", 50), rep("Suta", 50)))
allup <- rbind(ends, brachy)

ggplot(allup, aes(Latitude, SVL, height=SVL, fill=Species)) +
  geom_ridgeline() + theme_bw()

################


ggplot(allup, aes(x = Latitude, y = SVL)) +
  geom_density_ridges(aes(fill = Species))

ggplot(iris, aes(x = Sepal.Length, y = Species)) +
  geom_density_ridges(aes(fill = Species)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))



ggplot(allup, aes(Latitude, Species,group=Species)) +
  geom_density_ridges()

ggplot(allup, aes(Latitude, SVL, col=Species)) +
  geom_point(size=3) + theme_bw()





#
