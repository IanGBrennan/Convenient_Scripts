library(geosphere)
library(gmapsdistance)



test <- gmapsdistance(origin="40.7128+74.0060", destination="35.2809+149.1300", mode="driving", key="...")

Lat1 <- 40.7128
Lon1 <- 74.0060
Lat2 <- -35.2809
Lon2 <- 149.1300

testdist1 <- data.frame(longitude=Lon1, latitude=Lat1)
testdist2 <- data.frame(longitude=Lon2, latitude=Lat2)
testdist <- rbind(testdist1, testdist2)

distm(x=testdist, y=testdist)


cdist <- read.csv("/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/Continental_Distance_Matrices.csv", header=T)


ages <- c(0,5,10,15,20,25,30,35,40,50,60,70,80,90,100)
dist.mat <- NULL
for (i in 1:length(ages)){
  curr.mat <- dplyr::filter(cdist, Age==ages[[i]])
  dist.mat[[i]] <- distm(x=curr.mat[,3:4], y=curr.mat[,3:4])
}

red.dist <- NULL
for (i in 1:length(dist.mat)){
  min.val <- min(dist.mat[[i]][which(dist.mat[[i]]>0)])
  red.dist[[i]] <- dist.mat[[i]] / min.val
}

tot.dist <- NULL
for (i in 1:length(red.dist)){
  tot.dist <- rbind(tot.dist, red.dist[[i]])
}

write.csv(tot.dist, file="/Users/Ian/Google.Drive/R.Analyses/Varanus_Project/Continental_Distances.csv" )



iqtrees <- read.tree("/Users/Ian/Desktop/Species_Tree/FINAL_FullSampling/IQTree_GeneTrees/IQTree_Gene.trees")
"Varanus.komodoensis" %in% iqtrees[[1]]$tip.label
lapply(iqtrees, function(x) {"Varanus_komodoensis_75731" %in% iqtrees[[x]]$tip.label})
kin <- sapply(1:length(iqtrees), function(x) {"Varanus_komodoensis_75731" %in% iqtrees[[x]]$tip.label})
kintrees <- iqtrees[which(kin==TRUE)]

write.tree(kintrees, file="/Users/Ian/Desktop/Species_Tree/FINAL_FullSampling/IQTree_GeneTrees/IQTree_KomodoGene.trees")

#



















cdist