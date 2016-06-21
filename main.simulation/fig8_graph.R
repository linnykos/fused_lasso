rm(list=ls())
load("~/ryan/fused.git/results/graph-dim_100_2016-06-16_reduced.RData")
setwd("~/ryan/fused.git")
source("source_header.R")


true.edge.mat = enumerate.jumps.2dgraph(beta0)
est.edge.mat = enumerate.jumps.2dgraph(beta.est)

y.clipped = y
y.clipped[y.clipped < quantile(y, 0.1)] = quantile(y, 0.1)
y.clipped[y.clipped > quantile(y, 0.9)] = quantile(y, 0.9)

pdf(paste0("plots/graph_", Sys.Date(), ".pdf"), height = 3, width = 9)

par(mfrow = c(1,3))
par(mar = c(0.25, 0.25, 2, 0.25))
color.vec = gray.colors(30)
zlim = range(c(beta0, beta.est, y.clipped))

plot.2dgraph(beta0, true.edge.mat = true.edge.mat, color.vec = color.vec,
  zlim = zlim)
title(main = "Mean")

image(y.clipped, col = color.vec, asp = T, xlab = "", ylab = "", yaxt = "n",
 xaxt = "n", bty = "n", zlim = zlim)
title(main = "Data")

plot.2dgraph(beta.est, true.edge.mat = true.edge.mat, est.edge.mat =
 est.edge.mat, zlim = zlim, color.vec = color.vec)
title(main = "Estimate")

graphics.off()
quit()
