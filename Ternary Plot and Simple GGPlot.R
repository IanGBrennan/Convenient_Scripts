#install.packages('ggtern')
library(ggtern)

####### Let's make a Ternary Plot (or try)
#############################################
data<-read.csv("Ternary.test.csv")
ecology<-data$Ecology #read in the ecology data
trunk<-data$trunk.total #read in the second frame (trunk/total length)
head<-data$head.total #read in the next fram (head/total length)
tail<-data$tail.total #read in the next frame (tail/total length)
head.eighteen<-data$head.eighteen
head.ten <- (data$head.total)*10

#Create a vector of the colors
colors<-c("lightgreen", "blue", "red", "pink", "orange", "darkgray")
names(colors) <- levels(ecology)
col<-colors[match(ecology,names(colors))]

#Create a dataframe of the above information
df<-data.frame(trunk=trunk,
               head=head,
               tail=tail)
plot<-ggtern(data=df,
       mapping=aes(x=trunk, y=head, z=tail)) + geom_point(color=col, size=2.5)
plot + tern_limits(L=1, T=0.2, R=1)

df<-data.frame(trunk=trunk,
               head=head.ten,
               tail=tail)
ggtern(data=df,
       mapping=aes(x=trunk, y=head, z=tail)) + geom_point(color=col, 
                             size=2.5)

####### Let's make a Ternary Plot of the Pygopodidae
#############################################
data<-read.csv("Pygo.Ternary.csv")
ecology<-data$Ecology #read in the ecology data
trunk<-data$trunk.total #read in the second frame (trunk/total length)
trunk<-data$log.trunk
head<-data$head.total #read in the next fram (head/total length)
head <- data$log.head
tail<-data$tail.total #read in the next frame (tail/total length)
tail <- data$log.tail
genus<-data$Genus
head.abs<-data$head
trunk.abs<-data$trunk
tail.abs<-data$tail
svl<-data$svl
logsvl<-data$log.svl
logtail<-data$log.tail

#Create a vector of the colors
colors<-c("lightgreen", "blue", "red", "black", "orange", "darkgray")
names(colors) <- levels(genus)
col<-colors[match(genus,names(colors))]

#Create a dataframe of the above information
df<-data.frame(trunk=trunk,
               head=head,
               tail=tail)
plot<-ggtern(data=df,
       mapping=aes(x=trunk, y=head, z=tail)) + geom_point(color=col, size=3)
plot + tern_limits(L=0.65, T=0.65, R=0.75) #define the max limits of the plot

plot +
  scale_T_continuous(limits=c(.1,.3)) +
  scale_R_continuous(limits=c(.4,.6)) +
  scale_L_continuous(limits=c(.3,.5))

ggtern(df,aes(trunk, head, tail)) + 
  geom_density_tern(aes(color=..level..),bins=5) +
  geom_point(color=col, size=3) +
  scale_T_continuous(limits=c(.1,.3)) +
  scale_R_continuous(limits=c(.4,.6)) +
  scale_L_continuous(limits=c(.3,.5))

plot +
  scale_T_continuous(breaks=c(0,.6),breaks=c(0,1,.1)) +
  scale_R_continuous(breaks=c(.4,1),breaks=c(0,1,.1)) +
  scale_L_continuous(breaks=c(0,.6),breaks=c(0,1,.1))


############################################
#### Create Two Dimension plot of X vs. Y
############################################
btdf<-data.frame(x=logsvl, y=logtail) #create dataframe from data above
#plot the object below, control specifics in the geom_point() section
ggplot(data=btdf,mapping=aes(logsvl, logtail))+geom_point(color=col, size=3)


#### Alternative Ternary Plot Functions
library(ade4)
triangle.plot(data$head, data$tail, data$trunk)

library(rgl)
plot3d(data$head, data$tail, data$trunk)