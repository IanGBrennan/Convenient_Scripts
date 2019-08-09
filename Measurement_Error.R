tree <- read.nexus("/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T392_Neobatrachus/Neobatrachus_RAxML_ALl_Loci.tre")
tree <- midpoint.root(tree)
plot(tree)

calibrations <- makeChronosCalib(tree,interactive=T)
test <- chronos(tree, model="relaxed", calibration=calibrations)
plot(test, cex=0.4)
vcv(test)
write.tree(test, "/Users/Ian/Google.Drive/ANU Herp Work/Lemmon Projects/T392_Neobatrachus/Neobatrachus_PL.tre")


test.data <- c(75, 82, 83, 78, 80, 88, 77, 82, 81, 81, 80, 75)
# lognormal mean: log((mean(data)^2)/(sqrt(sd(data)^2+mean(data)^2)))
log((mean(test.data)^2)/(sqrt((sd(test.data)^2)+(mean(test.data)^2))))
# lognormal standard deviation: sqrt(log((1+(sd(data)^2)/(mean(data)^2))))
sqrt(log(1+((sd(test.data)^2)/(mean(test.data)^2))))


test.data <- c(60,68.5,77)
mean(test.data)
sd(test.data)

mars.data <- read.csv("/Users/Ian/Google.Drive/R.Analyses/BayesTraits/BT.Australian.Marsupials.BL.csv", header=T, row.names=1)
    colnames(mars.data) <- c("minimum", "maximum", "average")
macro.data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/CrownHeight_Macropodinae_spMEANS_ME.csv",
                       header=T, row.names=1)
    colnames(macro.data) <- c("minimum", "maximum", "average")
#mars.data <- log(mars.data)
ME.df <- function(mes.df, transformation=c("natural", "base10")){
  test <- NULL
  for(i in 1:nrow(mes.df)){
    measurements <- c(mes.df[i,"minimum"], mes.df[i,"maximum"], mes.df[i,"average"])
    #value <- (log(mean(measurements))) - (0.5 * log(((mean(measurements)/sd(measurements))^2)+1))
    if(transformation=="natural"){
      value <- log((mean(measurements)^2)/(sqrt((sd(measurements)^2)+(mean(measurements)^2))))
    }
    else if(transformation=="base10"){
      value <- log10((mean(measurements)^2)/(sqrt((sd(measurements)^2)+(mean(measurements)^2))))
    }
    test <- rbind(test, value)
  }
  colnames(test) <- "lognormal_mean"
  mes.df <- cbind(mes.df, test)
  
  tost <- NULL
  for(i in 1:nrow(mes.df)){
    measurements <- c(mes.df[i,"minimum"], mes.df[i,"maximum"], mes.df[i,"average"])
    if(transformation=="natural"){
      standiv <- sqrt(log(1+((sd(measurements)^2)/(mean(measurements)^2))))
    }
    else if (transformation=="base10"){
      standiv <- sqrt(log10(1+((sd(measurements)^2)/(mean(measurements)^2))))
    }
    tost <- rbind(tost, standiv)
  }
  colnames(tost) <- "lognormal_sd"
  mes.df <- cbind(mes.df, tost)
  return(mes.df)
}

mars.data <- ME.df(mars.data, transformation="natural")
data.mv <- mars.data$lognormal_mean; names(data.mv) <- rownames(mars.data)
data.me <- mars.data$lognormal_sd; names(data.me) <- rownames(mars.data)

macro.data <- ME.df(macro.data, transformation="natural")
write.csv(macro.data, "/Users/Ian/Google.Drive/ANU Herp Work/Macropod_Dating/CrownHeight_Macropodinae_spMEANS_ME.csv")

