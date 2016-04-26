rm(list=ls())
library("genlasso")
source("source.R")
source("source_plotter.R")
date = Sys.Date()
add.info = ""

sigma = 1
n = 200
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

trials = 1000

dual.list = list(trials)
lambda = 5

for(i in 1:trials){
  set.seed(10*i)
  y = generate.problem(n, jump.mean.org, jump.location.org,  sigma)
  
  res = dualplot_suite(y, jump.mean.org, jump.location.org, lambda = 15)
  res = dualplot_suite(y, jump.mean.org, jump.location.org, lambda = 15, plot.res = FALSE)
  dual.list[[i]] = compute.dual(y,res$coef)
  
  if(i%%floor(trials/10)==0)  cat('*')
}

save(dual.list,file=paste0("dual_multipletrials_",date,".RData"))
