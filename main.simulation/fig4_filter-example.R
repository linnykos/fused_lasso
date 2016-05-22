rm(list=ls())
load("~/ryan/fused.git/results/final-ROC2-2016-05-13.RData")
setwd("~/ryan/fused.git/")
source("source_header.R")

i = 5
trial = 1

truth = form.truth(jump.mean, jump.location, n.vec[i])
set.seed(i*trial*10)
y = generate.problem(n.vec[i], jump.mean, jump.location, sigma)

res = fusedlasso1d(y)
cv = cv.trendfilter(res, verbose = F)

fit = coef(res, lambda = cv$lambda.1se)$beta

threshold = filter.mat[1,which(colnames(filter.mat) == paste0(haus.quant))]

#png(file = paste0("~/DUMP/filter-example-", Sys.Date(), ".png"), width = 6,
# height = 4, units = "in", res = 300)
pdf(file = paste0("~/DUMP/filter-example-", Sys.Date(), ".pdf"), width = 6,
  height = 4)


plotfused(jump.mean, jump.location, y, fit, truth = truth,
  filter.bandwidth = ceiling(0.25*log(n.vec[i])^2), threshold = 
  threshold, plotDual = F)

dev.off()
