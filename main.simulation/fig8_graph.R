rm(list=ls())
load("~/ryan/fused.git/results/graph-dim_100_2016-06-01_reduced.RData")
setwd("~/ryan/fused.git")
source("source_header.R")


true.edge.mat = enumerate.jumps.2dgraph(beta0)
est.edge.mat = enumerate.jumps.2dgraph(beta.est)

pdf(paste0("plots/graph_", Sys.Date(), ".pdf"), height = 3, width = 9)

par(mfrow = c(1,3))
par(mar = c(1,1,4,1))

plot.2dgraph(beta0, true.edge.mat = true.edge.mat)
title(main = "True values")

image(y, col = gray.colors(30), asp = T, xlab = "", ylab = "", yaxt = "n",
 xaxt = "n", bty = "n")
title(main = "Observed values")

plot.2dgraph(beta.est, true.edge.mat = true.edge.mat, est.edge.mat =
 est.edge.mat)
title(main = "Estimated values")

dev.off()
