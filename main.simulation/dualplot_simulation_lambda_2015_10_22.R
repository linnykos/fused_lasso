rm(list=ls())
library("genlasso")
source("source.R")
source("source_plotter.R")
date = Sys.Date()
add.info = ""

sigma = 1
n.length = 10
n.vec = round(10^seq(2,3,length.out=n.length))
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

num.true.jumps = count.jumps(form.truth(jump.mean.org, jump.location.org, 100))
fixed.const = 3

res.list = list(n.length)
lambda.vec = numeric(n.length)

for(i in 1:n.length){
  set.seed(10*i)
  y = generate.problem(n.vec[i], jump.mean.org, jump.location.org,  sigma)
  
  png(paste("~/DUMP_figures/simul-",date,"_fixed.const-",fixed.const,"_n-",n.vec[i],add.info,".png",sep=""), width=8, height=3, units="in", res=300)
  res.list[[i]] = dualplot_suite(y, jump.mean.org, jump.location.org, num.est.jumps=fixed.const*num.true.jumps)
  dev.off()
  
  lambda.vec[i] = res.list[[i]]$lambda
  cat('*')
}


png(paste("~/DUMP_figures/simul-",date,"_lambda-vs-n",add.info,".png",sep=""), width=4,height=4,units="in",res=300)
plot(x=1:n.length, y=lambda.vec)
lines(x=1:n.length, y=lambda.vec, lwd=2)
dev.off()



