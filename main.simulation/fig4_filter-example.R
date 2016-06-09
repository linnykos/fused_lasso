rm(list=ls())
setwd("~/ryan/fused.git/")

load("results/final-ROC2-2016-05-21.RData")
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

pdf(file = paste0("plots/filter-example-", Sys.Date(), ".pdf"), width = 4.5,
  height = 5)

plotfused(jump.mean, jump.location, y, fit, truth = truth,
  filter.bandwidth = ceiling(0.25*log(n.vec[i])^2), threshold = 
  threshold, plotDual = F)

dev.off()
