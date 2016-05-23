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
interpolant = lower.interpolant(fit, truth, include.demean = T)


#make an arbitrary example
n = floor(n.vec[i]/2)*2
truth2 = c(rep(0, n/2), rep(5, n/2))

#add non-iid noise
set.seed(10)
tmp = (1:(n/2)-n/4)^2
tmp = tmp/(max(tmp)/2)
res.mean = rep(tmp - max(tmp)/2, 2)
y2 = truth2 + rnorm(n, res.mean, 1)

res2 = fusedlasso1d(y2)
cv2 = cv.trendfilter(res2, verbose = F)

fit2 = coef(res2, lambda = cv2$lambda.1se)$beta

#construct the lowerinterpolant
interpolant2 = lower.interpolant(fit2, truth2, include.demean = T)

#form the figure
pdf(file = paste0("plots/interpsim-", Sys.Date(), ".pdf"), width = 9,
  height = 5)

par(mfcol = c(2,2))
par(mar = c(2,4,1,1))

#first plot (primal)
plot(y, col=rgb(.5,.5,.5), pch=16, cex=1.25, ylab = "Value")
.plot.primal(jump.mean, jump.location, y, fit, tol = 1e-6, verbose = F)

#second plot (residuals)
res.split = split.signal(interpolant$demeaned)
int.split = split.signal(interpolant$lower.int)

plot(NA, xlim = c(0, n.vec[i]), ylim = c(min(res.split$mean, int.split$mean),
 max(res.split$mean, int.split$mean)), ylab = "Residual Value")

truth.split = split.signal(truth)
for(i in 2:length(truth.split$location)){
  lines(x = rep(truth.split$location[i], 2), y = c(-100, 100), lty = 2,
   lwd = 2, col = 2)
}

.plot.helper(res.split$location, res.split$mean, n.vec[i], col = 3)
.plot.helper(int.split$location, int.split$mean, n.vec[i], lty = 2)

#third plot for the second simulation
plot(y2, col=rgb(.5,.5,.5), pch=16, cex=1.25, ylab = "Value")
.plot.primal(NA, NA, y2, fit2, truth = truth2,
  tol = 1e-6, verbose = F)

#fourth plot
res2.split = split.signal(interpolant2$demeaned)
int2.split = split.signal(interpolant2$lower.int)

plot(NA, xlim = c(0, n), ylim = c(min(res2.split$mean, int2.split$mean), 
 max(res2.split$mean, int2.split$mean)), ylab = "Residual Value")

truth2.split = split.signal(truth2)
for(i in 2:length(truth2.split$location)){
  lines(x = rep(truth2.split$location[i], 2), y = c(-100, 100), lty = 2,
   lwd = 2, col = 2)
}

.plot.helper(res2.split$location, res2.split$mean, n, col = 3)
.plot.helper(int2.split$location, int2.split$mean, n, lty = 2)


dev.off()
