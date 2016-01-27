rm(list=ls())
source("source_header.R")

sigma = 1
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

setup = list(sigma = sigma, n.vec = n.vec, jump.mean = jump.mean.org, 
             jump.location = jump.location.org)

dist.mat = matrix(0, ncol = trials, nrow = n.length)
lambda.mat = matrix(0, ncol = trials, nrow = n.length)
mse.mat = matrix(0, ncol = trials, nrow = n.length)
haus.mat = matrix(0, ncol = trials, nrow = n.length)
jumps.mat = matrix(0, ncol = trials, nrow = n.length)

for(i in 1:n.length){
  truth = form.truth(jump.mean.org, jump.location.org, n.vec[i])
  true.jumps = enumerate.jumps(truth)
  
  for(trial in 1:trials){
    set.seed(i*trial*10)
    
    y = generate.problem(n.vec[i], jump.mean.org, jump.location.org,  sigma)
    
    res = fusedlasso1d(y)
    cv = cv.trendfilter(res, verbose = FALSE)
    
    fit = coef(res, lambda=cv$lambda.1se)$beta
    
    dist.mat[i,trial] = compute.hausdorff(true.jumps, enumerate.jumps(fit), 
                                          one.sided = TRUE)
    haus.mat[i,trial] = compute.hausdorff(true.jumps, enumerate.jumps(fit))
    lambda.mat[i,trial] = cv$lambda.1se
    mse.mat[i,trial] = compute.mse(fit, true.seq = truth)
    jumps.mat[i,trial] = count.jumps(fit)

    res = list(dist.mat = dist.mat, lambda.mat = lambda.mat, mse.mat = mse.mat,
     haus.mat = haus.mat, jumps.mat = jumps.mat, setup = setup)
    save(res, file=paste0("~/DUMP/screening_experiment_", DATE, ".RData"))
  }
  
  cat('*')
}


lambda.mat = res$lambda.mat
n.vec = res$setup$n.vec
logn.vec = log(n.vec, base = 10)

#remove uncomputed rows
idx = which(apply(lambda.mat, 1, function(x){all(x==0)}))
if(length(idx)) {
  lambda.mat = lambda.mat[-idx,]
  n.vec = n.vec[-idx]
}

#take median
lambda.vec = apply(lambda.mat, 1, function(x){median(x[x!=0])})

png("~/DUMP/lambda_2016-01-14.png", height = 5, width = 5, units = "in", res = 300)
plot(x = n.vec, y = lambda.vec, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylab = "Median Lambda value (1-se rule, 5-fold cv)")

#form the best sqrt(n) fit
x = seq(min(n.vec), max(n.vec), length.out = 1000)
const = max(lambda.vec)/max(sqrt(x))
y = const*sqrt(x)
lines(x, y, col="red", lwd = 2)
dev.off()

mse.mat = res$mse.mat
if(length(idx)) mse.mat = mse.mat[-idx,]
mse.vec = apply(mse.mat, 1, function(x){median(x[x!=0])})

png("~/DUMP/mse_2016-01-14.png", height = 5, width = 5, units = "in", res = 300)
plot(x = n.vec, y = mse.vec, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylab = "Median MSE value (1-se rule, 5-fold cv)")
const = min(mse.vec)/(min(log(x)*log(log(x))/x))
y = const*log(x)*log(log(x))/x
lines(x, y, col="red", lwd = 2)
dev.off()
