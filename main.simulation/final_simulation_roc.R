#load sourcese
rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

#set up parameters
sigma = 2
n.length = 10
n.vec = round(10^seq(2,4,length.out=n.length))
trials = 500
jump.mean = c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)
haus.quant = 0.8
i = 5

registerDoMC(cores = 20)

#testing function
run.test <- function(y, truth){
  res = fusedlasso1d(y)
  cv = cv.trendfilter(res, verbose = FALSE)

  fit = coef(res, lambda=cv$lambda.min)$beta

  n = length(y)
  filter.bandwidth = ceiling(0.25*(log(n))^2)

  filter.quant = staircase.threshold(y, fit,
   filter.bandwidth, cv$lambda.1se, controls = list(type =
    "original", quant = seq(0, 1, length.out = 101), 
    trials = 500))

  filter.quant
}

filter.mat = matrix(0, ncol = 101, nrow = trials)
colnames(filter.mat) = seq(0, 1, length.out = 101)


#run the simulations
truth = form.truth(jump.mean, jump.location, n.vec[i])

for(trial in 1:trials){
  set.seed(i*trial*10)
  y = generate.problem(n.vec[i], jump.mean, jump.location,  sigma)

  res = run.test(y, truth)

  #unpack the results
  filter.mat[trial,] = res

  save.image(file = paste0("results/final-ROC-", Sys.Date(), ".RData"))
 
}


