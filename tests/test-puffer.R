#more of just simple experiments than tests
rm(list=ls())
source("source_header.R")

sigma = 0.2
n.length = 10
n.vec = round(10^seq(2, 4, length.out = n.length))
trials = 50
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

i = 6
trial = 1
truth = form.truth(jump.mean, jump.location, n.vec[i])
true.jumps = enumerate.jumps(truth)

set.seed(i*trial*10)

y = generate.problem(n.vec[i], jump.mean, jump.location, sigma)

res = puffer.estimate(y)
