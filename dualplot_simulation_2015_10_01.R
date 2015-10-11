rm(list=ls())
library("genlasso")
source("dualplot_source.R")
date = "2015_10_01"

sigma = 1
n.length = 10
n.vec = seq(100,1000,length.out=n.length)
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

#scale signal at log n, lambda at sqrt(n log n)
for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.mean.org*log(n.vec[i]), jump.location.org,  sigma)
  png(paste("figures/simul_",date,"_lambda-sqrtnlogn_scale-logn_n-",n.vec[i],".png",sep=""), width=8, height=3, units="in", res=300)
  dualplot_suite(y, jump.mean.org*log(n.vec[i]), jump.location.org, lambda=sqrt(n.vec[i]*(log(n.vec[i])+1)*sigma)) 
  dev.off()

  print("*")
}

#scale signal at log n, lambda at sqrt(log n)
for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.mean.org*log(n.vec[i]), jump.location.org,  sigma)
  png(paste("figures/simul_",date,"_lambda-sqrtlogn_scale-logn_n-",n.vec[i],".png",sep=""), width=8, height=3, units="in", res=300)
  dualplot_suite(y,jump.mean.org*log(n.vec[i]), jump.location.org, lambda=sqrt((log(n.vec[i])+1)*sigma))
  dev.off()

  print("*")
}

#don't scale signal, lambda at sqrt(log n)
for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.mean.org,  jump.location.org, sigma)
  png(paste("figures/simul_",date,"_lambda-sqrtlogn_n-",n.vec[i],".png",sep=""), width=8, height=3, units="in", res=300)
  dualplot_suite(y, jump.mean.org, jump.location.org, lambda=sqrt((log(n.vec[i])+1)*sigma))
  dev.off()

  print("*")
}


