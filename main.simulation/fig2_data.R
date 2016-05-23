rm(list=ls())
setwd("~/ryan/fused.git")
source("source_header.R")

#set up parameters
sigma = 2
n.length = 10
n.vec = round(10^seq(2,4,length.out=n.length))
trials = 50
jump.mean = c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

pdf(file = paste0("plots/data-", Sys.Date(), ".pdf"), width = 9,
  height = 2.5)

par(mfcol = c(1,2))
par(mar = c(2,4,1,1))

i.vec = c(5,10)
trial = 1

for(i in i.vec){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  true.jumps = enumerate.jumps(truth)

  set.seed(i*trial*10)
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

  plot(y, col=rgb(.5,.5,.5), pch=16, cex=1.25, ylab = "Value")

  truth.split = split.signal(truth)
  .plot.helper(truth.split$location, truth.split$mean, n.vec[i], col = 3)
}

dev.off()
