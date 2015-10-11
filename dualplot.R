rm(list=ls())
library("genlasso")
source("dualplot_source.R")

set.seed(10)
n = 200
sigma = 1
dist = 1
jump.mean =     c(0, 2,  4, 1, 4)*dist
jump.location = c(0, .2, .4, .6, .8)

#denoised
y = generate.problem(n, jump.mean, jump.location,  sigma)
res = dualplot_suite(y,jump.mean, jump.location,num.est.jumps=12)

#standard scenario
set.seed(10)
y.org = generate.problem(n, jump.location, jump.mean, sigma)
dualplot_suite(y.org,jump.mean,jump.location,3)

#rescaled
dualplot_suite(y.org*2,jump.mean*2,jump.location,3)

dualplot_suite(y.org+10,jump.mean+10,jump.location,3)

#for experiment, shift the largest 20 points
ordering = order(y,decreasing=TRUE)
y.mod = y.org
select.num = 10
y.mod[ordering[1:select.num]] = y.org[ordering[1:select.num]]+20
dualplot_suite(y.mod,jump.mean,jump.location,3) #jump.mean is not necessarily meaningful

#re-center y.org and then rescale
y.mean = mean(y.org)
y.centered = y.org - y.mean
dualplot_suite(2*y.centered,2*(jump.mean-y.mean),jump.location,3)

#try rescaling the mean and noise but a sqrt level
set.seed(10)
random.seq = rnorm(n,sd=sigma)
random.seq = random.seq - mean(random.seq)
baseline = generate.problem(n, jump.location, jump.mean, sigma = 0)
baseline.mean = mean(baseline)
baseline = baseline - baseline.mean
dualplot_suite(baseline+random.seq,3,jump.mean-baseline.mean,jump.location)

dualplot_suite(2*baseline+sqrt(2)*random.seq,2*(jump.mean-baseline.mean),jump.location,3)
dualplot_suite(2*baseline+random.seq,2*(jump.mean-baseline.mean),jump.location,3)

sigma = 1
n.length = 5
n.vec = seq(100,500,length.out=n.length)
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org = c(0, .2, .4, .6, .8)

for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.location.org, jump.mean.org*log(n.vec[i]), sigma)
  dualplot_suite(y, jump.mean.org*log(n.vec[i]),jump.location.org, sqrt(n.vec[i]*(log(n.vec[i])+1)*sigma)) 
}


