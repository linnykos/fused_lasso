#remove the "junk" in the workspace of the results to make life easier

rm(list=ls())
load("~/ryan/fused.git/results/graph-dim_100_2016-06-01.RData")

truth = as.numeric(beta0)
mse.vec = apply(out$beta, 2, function(x){sum((x-truth)^2)})
idx = which.min(mse.vec)
beta.est = matrix(out$beta[,idx], dim1, dim2)
lambda.est = out$lambda[idx]
df.est = out$df[idx]

rm(list = "out")
save.image(file = 
 "~/ryan/fused.git/results/graph-dim_100_2016-06-01_reduced.RData")
