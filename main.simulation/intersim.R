#load sourcese
rm(list=ls())
setwd("~/ryan/fused.git/")
source("source_header.R")

#set up parameters
sigma = 2
n.length = 10
n.vec = round(10^seq(2,4,length.out=n.length))
trials = 50
jump.mean = c(0, 2,  4, 1, 4)
jump.location = c(0, .2, .4, .6, .8)
haus.quant = 0.8

i = 5
trial = 1

truth = form.truth(jump.mean, jump.location, n.vec[i])
set.seed(i*trial*10)
y = generate.problem(n.vec[i], jump.mean, jump.location, sigma)

res = fusedlasso1d(y)
cv = cv.trendfilter(res, verbose = F)

fit = coef(res, lambda = cv$lambda.1se)$beta

#construct the lowerinterpolant
interpolant = lower.interpolant(fit, truth)

#make an arbitrary example
n = 100
test = rep(0,n)
test[1:5] = .4
test[6:8] = .5
test[9:20] = .1
test[21:25] = .2
test[26:30] = .05
test[31:35] = .15
test[36:50] = -.1
test[51:80] = -.2
test[60:63] = -.1
test[81:100] = -.05
test[96:100] = -.4


#form the figure
png(file = paste0("~/DUMP/interpsim-", Sys.Date(), ".png"), width = 4, height = 6, units = "in",
    res = 300)
par(mfrow=c(3,1))
par(mar=c(2,2,2,1))

#first plot (primal)
plot(y, col=rgb(.5,.5,.5), pch=16, cex=1.25, ylab = "Value")
.plot.primal(jump.mean, jump.location, y, fit, tol = 1e-6)

#second plot (residuals)
res = fit - truth
res.split = split.signal(res)
int.split = split.signal(interpolant)

plot(NA, xlim = c(0, n.vec[i]), ylim = c(min(res), max(res)), ylab = "Residual Value")
.plot.helper(res.split$location, res.split$mean, length(res), col = "green")
.plot.helper(int.split$location, int.split$mean, length(res), lty = 2)


