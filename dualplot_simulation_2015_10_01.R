rm(list=ls())
library("genlasso")
source("dualplot_source.R")
date = "2015_10_01"

sigma = 1
n.length = 10
n.vec = seq(100,1000,length.out=n.length)
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)


for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.location.org, jump.mean.org*log(n.vec[i]), sigma)
  png(paste("simul_",date,"_lambda-sqrtnlogn_scale-logn_n-",n.vec[i],".png",sep=""), width=8, height=3, units="in", res=300)
  dualplot_suite(y, sqrt(n.vec[i]*(log(n.vec[i])+1)*sigma),jump.mean.org*log(n.vec[i]),jump.location.org) 
  dev.off()
}

