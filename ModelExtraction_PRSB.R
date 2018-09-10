bird.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Meliphagoids_ModelObjects.RDS")
pygo.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea_ModelObjects.RDS")
agam.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Agamidae_ModelObjects.RDS")
skink.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Sphenomorphines_ModelObjects.RDS")
mars.res <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupial_ModelObjects.RDS")

agam.i <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Agamidae_BMOUi_ModelObjects.RDS")
pygo.i <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Pygopodoidea_BMOUi_ModelObjects.RDS")
skink.i <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Sphenomorphines_BMOUi_ModelObjects.RDS")
mars.i <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Marsupials_BMOUi_ModelObjects.RDS")
bird.i <- readRDS("/Users/Ian/Google.Drive/ANU Herp Work/Adaptive Radiation/Body Size Model LOOP/Model.Fit.ME/Meliphagoids_BMOUi_ModelObjects.RDS")

extract.parameter <- function(results.object.model, method=c("newModels", "geiger", "mvMORPH"),
                              parameter=c("alpha", "sigma", "ME", "shift.time", "scalar", "ps.sigma",
                                          "stationary.var", "DC")) {
  param.est <- NULL
  
  if (method=="newModels"){
    for (k in 1:length(results.object.model)){
      if (parameter=="alpha"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$alpha)
      } else if (parameter=="sigma"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$beta)
      } else if (parameter=="ME"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$MErr)
      } else if (parameter=="shift.time"){
        param.est <- append(param.est, results.object.model[[k]]$shift.time)
      } else if (parameter=="scalar"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$post.shift.scalar)
      } else if (parameter=="ps.sigma"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$beta/results.object.model[[k]]$Trait1$post.shift.scalar)
      } else if (parameter=="stationary.var"){
        param.est <- append(param.est, results.object.model[[k]]$Trait1$beta/(2*results.object.model[[k]]$Trait1$alpha))
      }
    }
  } else if (method=="mvMORPH"){
    for (k in 1:length(results.object.model)){
      if (parameter=="alpha"){
        param.est <- append(param.est, results.object.model[[k]]$alpha)
      } else if (parameter=="sigma"){
        param.est <- append(param.est, results.object.model[[k]]$sigma)
      } else if (parameter=="ps.sigma"){
        param.est <- append(param.est, results.object.model[[k]]$sig)
      }else if (parameter=="stationary.var"){
        param.est <- append(param.est, results.object.model[[k]]$sig/(2*results.object.model[[k]]$alpha))
      } else if (parameter=="shift.time"){
        param.est <- append(param.est, results.object.model[[k]]$shift.time)
      }
    }
  } else if (method=="geiger"){
    for (k in 1:length(results.object.model)){
      if (parameter=="alpha"){
        param.est <- append(param.est, results.object.model[[k]]$opt$alpha)
      } else if (parameter=="DC"){
        param.est <- append(param.est, results.object.model[[k]]$opt$a)
      } else if (parameter=="sigma"){
        param.est <- append(param.est, results.object.model[[k]]$opt$sigsq)
      }
    }
  }
  return(param.est)
}

agam.alpha <- extract.parameter()


test <- extract.parameter(agam.bmoui$BMOUi, parameter="sigma", method="mvMORPH")
rate.sigma <- extract.parameter(pygo.bmoui$BMOUi, parameter="sigma", method="mvMORPH")
    test <- extract.parameter(pygo.bmoui$BMOUi, parameter="ps.sigma", method="mvMORPH")
    alpha <- extract.parameter(pygo.bmoui$BMOUi, parameter="alpha", method="mvMORPH")
    test/(2*alpha)
    s.sigma <- extract.parameter(pygo.bmoui$BMOUi, parameter="stationary.sigma", method="mvMORPH")
bm.sigma <- extract.parameter(skink.bmoui$BMOUi, parameter="sigma", method="mvMORPH")
st.sigma <- extract.parameter(skink.bmoui$BMOUi, parameter="stationary.sigma", method="mvMORPH")
bm.sigma <- extract.parameter(mars.i$BMOUi, parameter="sigma", method="mvMORPH")
st.sigma <- extract.parameter(mars.i$BMOUi, parameter="stationary.sigma", method="mvMORPH")

eb.sigma <- extract.parameter(bird.res$EB, parameter="sigma", method="geiger")
eb.dc <- extract.parameter(bird.res$EB, parameter="DC", method="geiger")


SRC_alpha <- extract.parameter(bird.res$SRC, parameter="alpha")
TRC_alpha <- extract.parameter(bird.res$TRC, parameter="alpha")
SRC_sigma <- extract.parameter(bird.res$SRC, parameter="sigma")
SRC_ME <- extract.parameter(bird.res$SRC, parameter="ME")
TRC_scalar <- extract.parameter(bird.res$TRC, parameter="scalar")

agam.sigma <- extract.parameter(agam.res$TRC, parameter="sigma")
agam.ps.sigma <- extract.parameter(agam.res$TRC, parameter="ps.sigma")

bird.alpha <- extract.parameter(bird.res$SRC, parameter="alpha", method="newModels")
  bird.alpha <- as.data.frame(bird.alpha); colnames(bird.alpha) <- "alpha"; bird.alpha$group <- "Meliphagoid"
agam.alpha <- extract.parameter(agam.res$SRC, parameter="alpha", method="newModels")
  agam.alpha <- as.data.frame(agam.alpha); colnames(agam.alpha) <- "alpha"; agam.alpha$group <- "Agamid"
pygo.alpha <- extract.parameter(pygo.res$SRC, parameter="alpha", method="newModels")
  pygo.alpha <- as.data.frame(pygo.alpha); colnames(pygo.alpha) <- "alpha"; pygo.alpha$group <- "Pygopodoid"
mars.alpha <- extract.parameter(mars.res$SRC, parameter="alpha", method="newModels")
  mars.alpha <- as.data.frame(mars.alpha); colnames(mars.alpha) <- "alpha"; mars.alpha$group <- "Marsupial"
skink.alpha <- extract.parameter(skink.res$SRC, parameter="alpha", method="newModels")
  skink.alpha <- as.data.frame(skink.alpha); colnames(skink.alpha) <- "alpha"; skink.alpha$group <- "Sphenomorphine"
all.alpha <- rbind(bird.alpha, agam.alpha, pygo.alpha, mars.alpha, skink.alpha)
  all.alpha <- all.alpha[which(all.alpha[,1] < 1),]
(ggplot(all.alpha, aes(x = alpha, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))
  
bird.ME <- extract.parameter(bird.res$SRC, parameter="stationary.var", method="newModels")
    bird.ME <- as.data.frame(bird.ME); colnames(bird.ME) <- "ME"; bird.ME$group <- "Meliphagoid"
agam.ME <- extract.parameter(agam.res$SRC, parameter="stationary.var", method="newModels")
    agam.ME <- as.data.frame(agam.ME); colnames(agam.ME) <- "ME"; agam.ME$group <- "Agamid"
pygo.ME <- extract.parameter(pygo.res$SRC, parameter="stationary.var", method="newModels")
    pygo.ME <- as.data.frame(pygo.ME); colnames(pygo.ME) <- "ME"; pygo.ME$group <- "Pygopodoid"
mars.ME <- extract.parameter(mars.res$SRC, parameter="stationary.var", method="newModels")
    mars.ME <- as.data.frame(mars.ME); colnames(mars.ME) <- "ME"; mars.ME$group <- "Marsupial"
skink.ME <- extract.parameter(skink.res$SRC, parameter="stationary.var", method="newModels")
    skink.ME <- as.data.frame(skink.ME); colnames(skink.ME) <- "ME"; skink.ME$group <- "Sphenomorphine"
all.ME <- rbind(bird.ME, agam.ME, pygo.ME, mars.ME, skink.ME)
(ggplot(all.ME, aes(x = ME, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))

bird.sigma <- extract.parameter(bird.res$SRC, parameter="sigma", method="newModels")
    bird.sigma <- as.data.frame(bird.sigma); colnames(bird.sigma) <- "sigma.bm"; bird.sigma$group <- "Meliphagoid"
agam.sigma <- extract.parameter(agam.res$SRC, parameter="sigma", method="newModels")
    agam.sigma <- as.data.frame(agam.sigma); colnames(agam.sigma) <- "sigma.bm"; agam.sigma$group <- "Agamid"
pygo.sigma <- extract.parameter(pygo.res$SRC, parameter="sigma", method="newModels")
    pygo.sigma <- as.data.frame(pygo.sigma); colnames(pygo.sigma) <- "sigma.bm"; pygo.sigma$group <- "Pygopodoid"
mars.sigma <- extract.parameter(mars.res$SRC, parameter="sigma", method="newModels")
    mars.sigma <- as.data.frame(mars.sigma); colnames(mars.sigma) <- "sigma.bm"; mars.sigma$group <- "Marsupial"
skink.sigma <- extract.parameter(skink.res$SRC, parameter="sigma", method="newModels")
    skink.sigma <- as.data.frame(skink.sigma); colnames(skink.sigma) <- "sigma.bm"; skink.sigma$group <- "Sphenomorphine"
all.sigma <- rbind(bird.sigma, agam.sigma, pygo.sigma, mars.sigma, skink.sigma)
    #all.sigma <- all.ME[which(all.ME[,1] < 1),]
(ggplot(all.sigma, aes(x = sigma.bm, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))

bird.sigma <- extract.parameter(bird.res$TRC, parameter="sigma")
bird.sigma <- as.data.frame(bird.sigma); colnames(bird.sigma) <- "sigma"; bird.sigma$group <- "Meliphagoid"
agam.sigma <- extract.parameter(agam.res$TRC, parameter="sigma")
agam.sigma <- as.data.frame(agam.sigma); colnames(agam.sigma) <- "sigma"; agam.sigma$group <- "Agamid"
pygo.sigma <- extract.parameter(pygo.res$TRC, parameter="sigma")
pygo.sigma <- as.data.frame(pygo.sigma); colnames(pygo.sigma) <- "sigma"; pygo.sigma$group <- "Pygopodoid"
mars.sigma <- extract.parameter(mars.res$TRC, parameter="sigma")
mars.sigma <- as.data.frame(mars.sigma); colnames(mars.sigma) <- "sigma"; mars.sigma$group <- "Marsupial"
skink.sigma <- extract.parameter(skink.res$TRC, parameter="sigma")
skink.sigma <- as.data.frame(skink.sigma); colnames(skink.sigma) <- "sigma"; skink.sigma$group <- "Sphenomorphine"
all.sigma <- rbind(bird.sigma, agam.sigma, pygo.sigma, mars.sigma, skink.sigma)
#all.sigma <- all.ME[which(all.ME[,1] < 1),]
(ggplot(all.sigma, aes(x = sigma, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=3, ncol=2, scales="free"))

bird.ps.sigma <- extract.parameter(bird.res$TRC, parameter="ps.sigma")
    bird.ps.sigma <- as.data.frame(bird.ps.sigma); colnames(bird.ps.sigma) <- "ps.sigma"; bird.ps.sigma$group <- "Meliphagoid"
agam.ps.sigma <- extract.parameter(agam.res$TRC, parameter="ps.sigma")
    agam.ps.sigma <- as.data.frame(agam.ps.sigma); colnames(agam.ps.sigma) <- "ps.sigma"; agam.ps.sigma$group <- "Agamid"
pygo.ps.sigma <- extract.parameter(pygo.res$TRC, parameter="ps.sigma")
    pygo.ps.sigma <- as.data.frame(pygo.ps.sigma); colnames(pygo.ps.sigma) <- "ps.sigma"; pygo.ps.sigma$group <- "Pygopodoid"
mars.ps.sigma <- extract.parameter(mars.res$TRC, parameter="ps.sigma")
    mars.ps.sigma <- as.data.frame(mars.ps.sigma); colnames(mars.ps.sigma) <- "ps.sigma"; mars.ps.sigma$group <- "Marsupial"
skink.ps.sigma <- extract.parameter(skink.res$TRC, parameter="ps.sigma")
    skink.ps.sigma <- as.data.frame(skink.ps.sigma); colnames(skink.ps.sigma) <- "ps.sigma"; skink.ps.sigma$group <- "Sphenomorphine"
all.ps.sigma <- rbind(bird.ps.sigma, agam.ps.sigma, pygo.ps.sigma, mars.ps.sigma, skink.ps.sigma)
#all.sigma <- all.ME[which(all.ME[,1] < 1),]
(ggplot(all.ps.sigma, aes(x = ps.sigma, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=3, ncol=2, scales="free"))

all.both <- cbind(all.sigma, all.ps.sigma[1])
(ggplot(all.both, aes(x = c(ps.sigma,sigma), fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=3, ncol=2, scales="free"))


bird.sigma <- extract.parameter(bird.i$BMOUi, parameter="sigma", method="mvMORPH")
    bird.sigma <- as.data.frame(bird.sigma); colnames(bird.sigma) <- "sigma.bm"; bird.sigma$group <- "Meliphagoid"
agam.sigma <- extract.parameter(agam.i$BMOUi, parameter="sigma", method="mvMORPH")
    agam.sigma <- as.data.frame(agam.sigma); colnames(agam.sigma) <- "sigma.bm"; agam.sigma$group <- "Agamid"
pygo.sigma <- extract.parameter(pygo.i$BMOUi, parameter="sigma", method="mvMORPH")
    pygo.sigma <- as.data.frame(pygo.sigma); colnames(pygo.sigma) <- "sigma.bm"; pygo.sigma$group <- "Pygopodoid"
mars.sigma <- extract.parameter(mars.i$BMOUi, parameter="sigma", method="mvMORPH")
    mars.sigma <- as.data.frame(mars.sigma); colnames(mars.sigma) <- "sigma.bm"; mars.sigma$group <- "Marsupial"
skink.sigma <- extract.parameter(skink.i$BMOUi, parameter="sigma", method="mvMORPH")
    skink.sigma <- as.data.frame(skink.sigma); colnames(skink.sigma) <- "sigma.bm"; skink.sigma$group <- "Sphenomorphine"
all.sigma <- rbind(bird.sigma, agam.sigma, pygo.sigma, mars.sigma, skink.sigma)
    #all.sigma <- all.sigma[which(all.sigma[,1] < 0.5),]
(ggplot(all.sigma, aes(x = sigma.bm, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))

bird.alpha <- extract.parameter(bird.i$BMOUi, parameter="alpha", method="mvMORPH")
    bird.alpha <- as.data.frame(bird.alpha); colnames(bird.alpha) <- "alpha"; bird.alpha$group <- "Meliphagoid"
agam.alpha <- extract.parameter(agam.i$BMOUi, parameter="alpha", method="mvMORPH")
    agam.alpha <- as.data.frame(agam.alpha); colnames(agam.alpha) <- "alpha"; agam.alpha$group <- "Agamid"
pygo.alpha <- extract.parameter(pygo.i$BMOUi, parameter="alpha", method="mvMORPH")
    pygo.alpha <- as.data.frame(pygo.alpha); colnames(pygo.alpha) <- "alpha"; pygo.alpha$group <- "Pygopodoid"
mars.alpha <- extract.parameter(mars.i$BMOUi, parameter="alpha", method="mvMORPH")
    mars.alpha <- as.data.frame(mars.alpha); colnames(mars.alpha) <- "alpha"; mars.alpha$group <- "Marsupial"
skink.alpha <- extract.parameter(skink.i$BMOUi, parameter="alpha", method="mvMORPH")
    skink.alpha <- as.data.frame(skink.alpha); colnames(skink.alpha) <- "alpha"; skink.alpha$group <- "Sphenomorphine"
all.alpha <- rbind(bird.alpha, agam.alpha, pygo.alpha, mars.alpha, skink.alpha)
    all.alpha <- all.alpha[which(all.alpha[,1] < 1),]
(ggplot(all.alpha, aes(x = alpha, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))
    
bird.stat <- extract.parameter(bird.i$BMOUi, parameter="stationary.var", method="mvMORPH")
    bird.stat <- as.data.frame(bird.stat); colnames(bird.stat) <- "stationary.var"; bird.stat$group <- "Meliphagoid"
agam.stat <- extract.parameter(agam.i$BMOUi, parameter="stationary.var", method="mvMORPH")
    agam.stat <- as.data.frame(agam.stat); colnames(agam.stat) <- "stationary.var"; agam.stat$group <- "Agamid"
pygo.stat <- extract.parameter(pygo.i$BMOUi, parameter="stationary.var", method="mvMORPH")
    pygo.stat <- as.data.frame(pygo.stat); colnames(pygo.stat) <- "stationary.var"; pygo.stat$group <- "Pygopodoid"
mars.stat <- extract.parameter(mars.i$BMOUi, parameter="stationary.var", method="mvMORPH")
    mars.stat <- as.data.frame(mars.stat); colnames(mars.stat) <- "stationary.var"; mars.stat$group <- "Marsupial"
skink.stat <- extract.parameter(skink.i$BMOUi, parameter="stationary.var", method="mvMORPH")
    skink.stat <- as.data.frame(skink.stat); colnames(skink.stat) <- "stationary.var"; skink.stat$group <- "Sphenomorphine"
all.stat <- rbind(bird.stat, agam.stat, pygo.stat, mars.stat, skink.stat)
    #all.stat <- all.alpha[which(all.alpha[,1] < 1),]
(ggplot(all.stat, aes(x = stationary.var, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))
    

bird.sigma <- extract.parameter(bird.res$EB, parameter="sigma", method="geiger")
    bird.sigma <- as.data.frame(bird.sigma); colnames(bird.sigma) <- "sigma.eb"; bird.sigma$group <- "Meliphagoid"
agam.sigma <- extract.parameter(agam.res$EB, parameter="sigma", method="geiger")
    agam.sigma <- as.data.frame(agam.sigma); colnames(agam.sigma) <- "sigma.eb"; agam.sigma$group <- "Agamid"
pygo.sigma <- extract.parameter(pygo.res$EB, parameter="sigma", method="geiger")
    pygo.sigma <- as.data.frame(pygo.sigma); colnames(pygo.sigma) <- "sigma.eb"; pygo.sigma$group <- "Pygopodoid"
mars.sigma <- extract.parameter(mars.res$EB, parameter="sigma", method="geiger")
    mars.sigma <- as.data.frame(mars.sigma); colnames(mars.sigma) <- "sigma.eb"; mars.sigma$group <- "Marsupial"
skink.sigma <- extract.parameter(skink.res$EB, parameter="sigma", method="geiger")
    skink.sigma <- as.data.frame(skink.sigma); colnames(skink.sigma) <- "sigma.eb"; skink.sigma$group <- "Sphenomorphine"
all.sigma <- rbind(bird.sigma, agam.sigma, pygo.sigma, mars.sigma, skink.sigma)
    #all.sigma <- all.sigma[which(all.sigma[,1] < 0.5),]
(ggplot(all.sigma, aes(x = sigma.eb, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))

bird.DC <- extract.parameter(bird.res$EB, parameter="DC", method="geiger")
    bird.DC <- as.data.frame(bird.DC); colnames(bird.DC) <- "DC"; bird.DC$group <- "Meliphagoid"
agam.DC <- extract.parameter(agam.res$EB, parameter="DC", method="geiger")
    agam.DC <- as.data.frame(agam.DC); colnames(agam.DC) <- "DC"; agam.DC$group <- "Agamid"
pygo.DC <- extract.parameter(pygo.res$EB, parameter="DC", method="geiger")
    pygo.DC <- as.data.frame(pygo.DC); colnames(pygo.DC) <- "DC"; pygo.DC$group <- "Pygopodoid"
mars.DC <- extract.parameter(mars.res$EB, parameter="DC", method="geiger")
    mars.DC <- as.data.frame(mars.DC); colnames(mars.DC) <- "DC"; mars.DC$group <- "Marsupial"
skink.DC <- extract.parameter(skink.res$EB, parameter="DC", method="geiger")
    skink.DC <- as.data.frame(skink.DC); colnames(skink.DC) <- "DC"; skink.DC$group <- "Sphenomorphine"
all.DC <- rbind(bird.DC, agam.DC, pygo.DC, mars.DC, skink.DC)
    #all.sigma <- all.sigma[which(all.sigma[,1] < 0.5),]
(ggplot(all.DC, aes(x = DC, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  #+ xlim(0,0.5)
  + facet_wrap(~group, nrow=5, ncol=1, scales="free"))


bird.s <- extract.parameter(bird.res$TRC, parameter="shift.time", method="newModels")
    bird.s <- as.data.frame(bird.s); colnames(bird.s) <- "time"; bird.s$group <- "Meliphagoid"
agam.s <- extract.parameter(agam.res$TRC, parameter="shift.time", method="newModels")
    agam.s <- as.data.frame(agam.s); colnames(agam.s) <- "time"; agam.s$group <- "Agamid"
pygo.s <- extract.parameter(pygo.res$TRC, parameter="shift.time", method="newModels")
    pygo.s <- as.data.frame(pygo.s); colnames(pygo.s) <- "time"; pygo.s$group <- "Pygopodoid"
mars.s <- extract.parameter(mars.res$TRC, parameter="shift.time", method="newModels")
    mars.s <- as.data.frame(mars.s); colnames(mars.s) <- "time"; mars.s$group <- "Marsupial"
skink.s <- extract.parameter(skink.res$TRC, parameter="shift.time", method="newModels")
    skink.s <- as.data.frame(skink.s); colnames(skink.s) <- "time"; skink.s$group <- "Sphenomorphine"
all.s <- rbind(bird.s, agam.s, pygo.s, mars.s, skink.s)
(ggplot(all.s, aes(x = time, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  + xlim(5,15)
  + facet_wrap(~group, nrow=5, ncol=1))

bird.s <- extract.parameter(bird.i$BMOUi, parameter="shift.time", method="mvMORPH")
    bird.s <- as.data.frame(bird.s); colnames(bird.s) <- "time"; bird.s$group <- "Meliphagoid"
agam.s <- extract.parameter(agam.i$BMOUi, parameter="shift.time", method="mvMORPH")
    agam.s <- as.data.frame(agam.s); colnames(agam.s) <- "time"; agam.s$group <- "Agamid"
pygo.s <- extract.parameter(pygo.i$BMOUi, parameter="shift.time", method="mvMORPH")
    pygo.s <- as.data.frame(pygo.s); colnames(pygo.s) <- "time"; pygo.s$group <- "Pygopodoid"
mars.s <- extract.parameter(mars.i$BMOUi, parameter="shift.time", method="mvMORPH")
    mars.s <- as.data.frame(mars.s); colnames(mars.s) <- "time"; mars.s$group <- "Marsupial"
skink.s <- extract.parameter(skink.i$BMOUi, parameter="shift.time", method="mvMORPH")
    skink.s <- as.data.frame(skink.s); colnames(skink.s) <- "time"; skink.s$group <- "Sphenomorphine"
all.s <- rbind(bird.s, agam.s, pygo.s, mars.s, skink.s)
(ggplot(all.s, aes(x = time, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  + xlim(5,15)
  + facet_wrap(~group, nrow=5, ncol=1))


bird.shift <- extract.parameter(bird.i$BMOUi, parameter="shift.time", method="mvMORPH")
    bird.shift <- as.data.frame(bird.shift); colnames(bird.shift) <- "shift"; bird.shift$group <- "Meliphagoid"
agam.shift <- extract.parameter(agam.i$BMOUi, parameter="shift.time", method="mvMORPH")
    agam.shift <- as.data.frame(agam.shift); colnames(agam.shift) <- "shift"; agam.shift$group <- "Agamid"
pygo.shift <- extract.parameter(pygo.i$BMOUi, parameter="shift.time", method="mvMORPH")
    pygo.shift <- as.data.frame(pygo.shift); colnames(pygo.shift) <- "shift"; pygo.shift$group <- "Pygopodoid"
mars.shift <- extract.parameter(mars.i$BMOUi, parameter="shift.time", method="mvMORPH")
    mars.shift <- as.data.frame(mars.shift); colnames(mars.shift) <- "shift"; mars.shift$group <- "Marsupial"
skink.shift <- extract.parameter(skink.i$BMOUi, parameter="shift.time", method="mvMORPH")
    skink.shift <- as.data.frame(skink.shift); colnames(skink.shift) <- "shift"; skink.shift$group <- "Sphenomorphine"
all.shift <- rbind(bird.shift, agam.shift, pygo.shift, mars.shift, skink.shift)
#all.sigma <- all.sigma[which(all.sigma[,1] < 0.5),]
(ggplot(all.shift, aes(x = shift, fill = group)) 
  + geom_density(alpha = 1)
  + theme(panel.background=element_blank(), legend.position="bottom")
  + scale_fill_manual(values=wes_palette("Darjeeling1", 5, "continuous"))
  + xlim(0,15))
  #+ facet_wrap(~group, nrow=5, ncol=1, scales="free"))
