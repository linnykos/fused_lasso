rm(list=ls())
source("source_header.R")

registerDoMC(cores = 12)

sigma = 1
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean, 
             jump.location = jump.location)

dist.mat = matrix(0, ncol = trials, nrow = n.length)
lambda.mat = matrix(0, ncol = trials, nrow = n.length)
mse.mat = matrix(0, ncol = trials, nrow = n.length)
haus.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.mat = matrix(0, ncol = trials, nrow = n.length)
haus.filter.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.filter.mat = matrix(0, ncol = trials, nrow = n.length)

simulation_suite <- function(trial){
  set.seed(i*trial*10)
    
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)
    
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = FALSE)
    
  fit = coef(res, lambda=cv$lambda.1se)$beta
    
  dist = compute.hausdorff(true.jumps, enumerate.jumps(fit), 
                                          one.sided = TRUE)
  haus = compute.hausdorff(true.jumps, enumerate.jumps(fit))
  lambda = cv$lambda.1se
  mse = compute.mse(fit, true.seq = truth)
  jumps = count.jumps(fit)
  
  filter.bandwidth = ceiling(0.25*(log(length(y)))^2)
  jumps.filter.idx = apply.filter(fit, filter.bandwidth, 0.5, return.type = "location")
  haus.filter = compute.hausdorff(true.jumps, jumps.filter.idx)
  jumps.filter = length(jumps.filter.idx)

  #plot
  png(paste0("~/DUMP/fused_lasso_n-", n.vec[i], "_trial-", trial, "_", DATE, ".png"), 
   width = 8, height = 3, units = "in", res = 300)
  plotfused(jump.mean, jump.location, y, fit, cv$lambda.1se, count.jumps(fit), 
            filter.bandwidth = filter.bandwidth, plotDual = F)
  dev.off()
 
  c(dist, haus, lambda, mse, jumps, haus.filter, jumps.filter)
}

for(i in 1:n.length){
  truth = form.truth(jump.mean, jump.location, n.vec[i])
  true.jumps = enumerate.jumps(truth)
 
  res.tmp = foreach(trial = 1:trials) %dopar% simulation_suite(trial)
  res.tmp = do.call(cbind, res.tmp)  

  dist.mat[i,] = res.tmp[1,]
  haus.mat[i,] = res.tmp[2,]
  lambda.mat[i,] = res.tmp[3,]
  mse.mat[i,] = res.tmp[4,]
  jumps.mat[i,] = res.tmp[5,]
  haus.filter.mat[i,] = res.tmp[6,]
  jumps.filter.mat[i,] = res.tmp[7,]

  res = list(dist.mat = dist.mat, lambda.mat = lambda.mat, mse.mat = mse.mat,
   haus.mat = haus.mat, jumps.mat = jumps.mat, setup = setup)
  save(res, file=paste0("~/DUMP/screening_experiment_", DATE, ".RData"))
 
  cat('*')
}


png("~/DUMP/truth_2016-01-27.png", height = 5, width = 10, units = "in", res = 300)
plot(y, col = rgb(.5,.5,.5), pch = 16, cex = 1.25)
plot.helper(c(1,true.jumps,n.vec[i]), jump.mean, n.vec[i], col="red")
dev.off()

lambda.mat = res$lambda.mat
n.vec = res$setup$n.vec
logn.vec = log(n.vec, base = 10)

#take median
lambda.vec = apply(lambda.mat, 1, mean)
lambda.std = apply(lambda.mat, 1, sd)

png("~/DUMP/screening_2016-01-14.png", height = 5, width = 10, units = "in", res = 300)

par(mfrow = c(1,2))

plot(x = n.vec, y = lambda.vec, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylab = "Median Lambda value (1-se rule, 5-fold cv)", ylim = c(0, max(lambda.vec+lambda.std)))
for(i in 1:length(lambda.std)) {
  lines(x = rep(n.vec[i], 2), y = c(lambda.vec[i] - lambda.std[i], lambda.vec[i] + lambda.std[i]),
   cex = 2)
}

#form the best sqrt(n) fit
x = seq(min(n.vec), max(n.vec), length.out = 1000)
const = max(lambda.vec)/max(sqrt(x))
y = const*sqrt(x)
lines(x, y, col="red", lwd = 2)
#dev.off()

mse.vec = apply(mse.mat, 1, mean)
mse.std = apply(mse.mat, 1, sd)

#png("~/DUMP/mse_2016-01-14.png", height = 5, width = 5, units = "in", res = 300)
plot(x = n.vec, y = mse.vec, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylab = "Median MSE value (1-se rule, 5-fold cv)", ylim = c(0, max(mse.vec + mse.std)))
for(i in 1:length(mse.std)) {
  lines(x = rep(n.vec[i], 2), y = c(mse.vec[i] - mse.std[i], mse.vec[i] + mse.std[i]),
   cex = 2)
}

const = min(mse.vec)/(min(log(x)*log(log(x))/x))
y = const*log(x)*log(log(x))/x
lines(x, y, col="red", lwd = 2)
dev.off()
