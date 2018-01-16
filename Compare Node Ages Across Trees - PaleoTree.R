library(paleotree)
citation("paleotree")


mcc <- read.nexus("Pygopodoidea.nuclear.all.1k.pb.tre")
compare <- read.nexus("Fossil.Calibrations.trees")
test <- compareNodeAges(mcc, compare, dropUnshared=T)
test<-as.data.frame(test)
write.table(test, file="Pygopodoidea.missing1.CompareNodeAge.xls", quote=F, sep="\t")

pb <- read.nexus("Original.Pygopodoidea.all.nuclear.100mil.trees")
test <- compareNodeAges(mcc, pb.files, dropUnshared=F)
test<-as.data.frame(test)
write.table(test, file="Pygopodoidea.500.pb.CompareNodeAge.xls", quote=F, sep="\t")

fish <- NULL
pygopodoidea <- as.data.frame(runif(100, -10,10))
mean(pygopodoidea[,1])
car.pygo <- as.data.frame(runif(100, -9, 10))
carpho <- as.data.frame(runif(100, -7,7))
pygo <- as.data.frame(runif(100, -5,5))
all.diplo <- as.data.frame(runif(100, -7,7))
diplo <- as.data.frame(runif(100, -5,5))
gekkota <- as.data.frame(runif(100, -2, 7))

fish <- cbind(fish, diplo[,1])
write.table(fish, file="NUMBERS.xls", quote=F, sep="\t")

#### Import the data for node ages
table <- read.csv("Pygopodoidea.missing.focal.nodes.csv", row.names = 1, header=T)
daters <- data.frame(x=(1:105), y=table)
treenames <- row.names(daters) #set the rownames for the colors
treenames[1:100] <- "pb" #set the posterior distribution values as the same color

#### Plot the fossil calibration simulations and the empirical posterior node ages
(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.crown.pygopodoidea, colour = treenames, size=4)))

(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.pygo.carpho, colour = treenames, size=4)))

(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.crown.pygo, colour = treenames, size=4)))

(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.crown.carpho, colour = treenames, size=4)))

(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.all.diplo, colour = treenames, size=4)))

(ggplot(dat=daters)
  + geom_point(aes(x=1, y=daters$y.crown.diplo, colour = treenames, size=4)))



var(test[,2])
