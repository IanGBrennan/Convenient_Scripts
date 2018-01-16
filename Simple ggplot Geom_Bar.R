library(ggplot2)

data <- read.csv("/Users/Ian/Google.Drive/ANU Herp Work/Cyrtodactylus/Cyrtodactylus Barcoding/Cyrtodactylus Publications by YEAR.csv")

colorz <- 
class(data$class)
data$class <- factor(data$class, c("hi", "medhi", "medium", "medmed", "medlow", "low"))
(ggplot(data, aes(x=year, weight=results, fill=class, order=class))
  + geom_bar())

