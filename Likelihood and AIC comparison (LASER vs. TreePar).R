library(BioGeoBEARS)

#### Method 1 #############################
###########################################
#LASER Output:use as is, calculate AIC with:
calc_AIC_vals(109.258, 7) #calc_AIC_vals(Likelihood, #parameters)
#TESS Output: gotta adjust it to the LASER output
comp<-(out[1])+sum(log(2:length(y))) #out[1] is the output from tess.likelihood estimate
calc_AIC_vals(comp, 7) #calc_AIC_vals(Laser-adjusted Likelihood, #parameters)
#TreePar Output: gotta adjust it to the LASER output
tp.out<--(p.threeshift[[4]][1])+sum(log(2:length(y))) #p.threeshift is the most likely treepar scenario
        #-(likelihood optimal treepar model) + sum(log(2:length(y))) ....make sure that the likelihood is negative
calc_AIC_vals(tp.out, 11) #calc_AIC_vals(Laser-adjusted Likelihood, #params)


#### Method 2 #############################
###########################################
#LASER Output: must adjust as below
comp<-(y4r[1])-sum(log(2:length(pygo)))
        #-(Likelihood of best model) - sum(log(2:length(y)))
calc_AIC_vals (comp, x)
#TESS Output: use as is, calculate AIC with:
calc_AIC_vals(out, x)
#TreePar Output: use as is, snd calculate AIC with -likelihood:
calc_AIC_vals(-(p.threeshift[[4]][1]), x)

# x is the number of parameters