library(DDD)

tree<-read.tree("oz.marsupials.tre")
mar<-getBtimes(file="oz.marsupials.tre")

bd_ML(mar)
ntdBD<-bd_loglik(mar, pars1=c(0.125, 0.072, 0, 0), pars2=c(0, 0, 0, 1, 2), missnumspec=0)
#loglik = 86.575042
eds<-bd_loglik(mar, pars1=c(0.125, 0.072, 0, 0), pars2=c(1, 0, 0, 1, 2), missnumspec=0)

dd_ML(mar)
ntdDD<-dd_loglik(brts=mar, pars1=c(0.126, 0.073, 10504), pars2=c(170, 1, 0, 0, 1, 2), missnumspec=0)
#logLik = 86.558643

tree<-read.tree("pygopodoidea.tre")
pygo<-getBtimes(file="pygopodoidea.tre")
bd_ML(pygo)
ntdBD<-bd_loglik(mar, pars1=c(0.074, 0.000015, 0, 0), pars2=c(0, 0, 0, 1, 2), missnumspec=0)
#logLik = 82.205358
dd_ML(pygo)
ntdDD<-dd_loglik(brts=mar, pars1=c(0.173605, 0.000099, 177.63), pars2=c(170, 1, 0, 0, 1, 2), missnumspec=0)
#logLik = 53.707000
fitdAICrc(pygo, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate", "yule3rate"), ints=120)
