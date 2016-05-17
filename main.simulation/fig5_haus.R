rm(list=ls())
load("~/ryan/fused.git/results/final-2016-05-15.RData")
setwd("~/ryan/fused.git/")
source("source_header.R")

compute.haus <- function(mat1, mat2){
  haus = pmax(mat1, mat2)
  med = apply(haus, 2, median, na.rm = T)
  upper = apply(haus, 2, quantile, prob = 0.75, na.rm = T)
  lower = apply(haus, 2, quantile, prob = 0.25, na.rm = T)

  list(median = med, upper = upper, lower = lower)
}

haus.list = vector("list", 3)
haus.list[[1]] = compute.haus(res.list$left.org, res.list$right.org)
haus.list[[2]] = compute.haus(res.list$left.oracle, res.list$right.oracle)
haus.list[[3]] = compute.haus(res.list$left.adapt, res.list$right.adapt)

png(file = paste0("~/DUMP/haus-dist-", Sys.Date(), ".png"), width = 4,
 height = 4, units = "in", res = 300)

par(mar = c(4,4,1,2))

plot(NA, xlim = c(0, max(n.vec)), ylim = c(0, max(haus.list[[1]]$upper)),
  xlab = "Number of observations", ylab = "Hausdorff Distance")

col.vec = c(1, 2, 4)
lty.vec = c(3, 1, 2)

for(i in 1:3){
  points(x = n.vec, y = haus.list[[i]]$med, pch = 16, col = col.vec[i],
   cex = 2)
  lines(x = n.vec, y = haus.list[[i]]$med, lty = lty.vec[i], col = col.vec[i],
   lwd = 3)
  
  for(j in 1:length(n.vec)){
    lines(x = rep(n.vec[j], 2), y = c(haus.list[[i]]$lower[j], 
     haus.list[[i]]$upper[j]), lty = lty.vec[i], col = col.vec[i], lwd = 2)
  }
}

dev.off()
