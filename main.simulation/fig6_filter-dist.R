rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

load("~/ryan/fused.git/results/final-ROC2-2016-05-13.RData")

mat1 = res.list$oracle.left
mat2 = res.list$oracle.right

png(paste0("~/DUMP/Haus_", Sys.Date(), ".png"), height = 5, width = 5,
 units = "in", res = 300)
plot.hausdorff(mat1, mat2)
dev.off()
