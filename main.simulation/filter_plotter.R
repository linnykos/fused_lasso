rm(list=ls())
setwd("..")
source("source_header.R")

#WARNING: MESSY CODE!!
library(assertthat)
load("~/ryan/fused.git/results/filterExperiment_hausdorff_2016-05-11.RData")
mat1 = res$oracle.left.list[[1]]
mat2 = res$oracle.right.list[[1]]

png(paste0("~/DUMP/Haus_", Sys.Date(), ".png"), height = 5, width = 5,
 units = "in", res = 300)
plot.hausdorff(mat1, mat2)
dev.off()


load("~/ryan/fused.git/results/filterExperiment_2016-05-11.RData")
sigma = res$setup$sigma
n.length = 1
n.vec = res$setup$n.vec
jump.mean = res$setup$jump.mean
jump.location = res$setup$jump.location

#plot data
i = 1
truth = form.truth(jump.mean, jump.location, n.vec[i])
true.jumps = enumerate.jumps(truth)

filter.bandwidth = ceiling(0.25*(log(length(truth)))^2)

trial = 1
set.seed(i*trial*10)
y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

res = fusedlasso1d(y)
cv = cv.trendfilter(res, verbose = FALSE)
fit = coef(res, lambda = cv$lambda.1se)$beta


