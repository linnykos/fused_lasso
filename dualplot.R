library("genlasso")
source("dualplot_source.R")

set.seed(10)
n = 100
sigma = 1
dist = 1
jump.mean =     c(0, 2,  4, 1, 4)*dist
jump.location = c(0, .2, .4, .6, .8)

#denoised
y = generate.problem(n, jump.location, jump.mean, 0)
dualplot_suite(y,3,jump.mean,jump.location)

#standard scenario
set.seed(10)
y.org = generate.problem(n, jump.location, jump.mean, sigma)
dualplot_suite(y.org,3,jump.mean,jump.location)

#rescaled
dualplot_suite(y.org*2,3,jump.mean*2,jump.location)

dualplot_suite(y.org+10,3,jump.mean+10,jump.location)

#for experiment, shift the largest 20 points
ordering = order(y,decreasing=TRUE)
y.mod = y.org
select.num = 10
y.mod[ordering[1:select.num]] = y.org[ordering[1:select.num]]+20
dualplot_suite(y.mod,3,jump.mean,jump.location) #jump.mean is not necessarily meaningful

#re-center y.org and then rescale
y.mean = mean(y.org)
y.centered = y.org - y.mean
dualplot_suite(2*y.centered,3,2*(jump.mean-y.mean),jump.location)

#try rescaling the mean and noise but a sqrt level
set.seed(10)
random.seq = rnorm(n,sd=sigma)
random.seq = random.seq - mean(random.seq)
baseline = generate.problem(n, jump.location, jump.mean, sigma = 0)
baseline.mean = mean(baseline)
baseline = baseline - baseline.mean
dualplot_suite(baseline+random.seq,3,jump.mean-baseline.mean,jump.location)

dualplot_suite(2*baseline+sqrt(2)*random.seq,3,2*(jump.mean-baseline.mean),jump.location)
dualplot_suite(2*baseline+random.seq,3,2*(jump.mean-baseline.mean),jump.location)

#use the fact that our theory shows lambda = sqrt(log(en)*sigma) should be good
sigma = 1
n.length = 5
n.vec = seq(100,500,length.out=n.length)
jump.mean =     c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)

for(i in 1:n.length){
  set.seed(10)
  y = generate.problem(n.vec[i], jump.location, jump.mean*log(n.vec[i]), sigma)
  dualplot_suite(y, sqrt(n.vec[i]*(log(n.vec[i])+1)*sigma),jump.mean*log(n.vec[i]),jump.location) 
}


