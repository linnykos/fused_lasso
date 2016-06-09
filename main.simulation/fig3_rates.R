setwd("~/ryan/fused.git")

load("results/final-2016-05-22.RData")
library(assertthat)
library(grDevices)

.plot.rates <- function(mat, theoretical.rate, 
 n.vec = as.numeric(colnames(mat)), func = mean){
  assert_that(class(theoretical.rate) == "function")
  assert_that(length(n.vec) == ncol(mat))

  med.vec = apply(mat, 2, func, na.rm = T)
  sd.vec  = apply(mat, 2, sd, na.rm = T)

  #use linear regression to find best fit
  scaling = coef(lm(med.vec ~ theoretical.rate(n.vec) - 1))
 
  n.vec.full = seq(min(n.vec), max(n.vec), length.out = 200)
  lines(x = n.vec.full, y = scaling*theoretical.rate(n.vec.full), 
   lwd = 2, col = 2)

  #median values
  points(x = n.vec, y = med.vec, cex = 1.25, pch = 16)

  #error bars
  for(i in 1:length(n.vec)){
    lines(x = rep(n.vec[i], 2), y = c(med.vec[i]-sd.vec[i], 
     med.vec[i]+sd.vec[i]), lwd = 2)
  }

  invisible()
}


pdf(file = paste0("plots/rates-", Sys.Date(), ".pdf"), width = 7.5,
  height = 2.5)
par(mfrow = c(1,3))
par(mar = c(4,4,1,1))

plot(NA, xlab = "n", ylab = expression(lambda), ylim = c(0, median(res.list$lambda[,10]) + 
 sd(res.list$lambda[,10])), xlim = c(0, max(n.vec)), cex.lab = 1.25, 
 cex.axis = 1.25)
.plot.rates(res.list$lambda, theoretical.rate = function(x){sqrt(x)},
 n = n.vec)

#make the plot of the MSE
plot(NA, xlab = "n", ylab = "MSE", ylim = c(0, median(res.list$mse[,1]) + sd(res.list$mse[,1])), 
 xlim = c(0, max(n.vec)), cex.lab = 1.25, cex.axis = 1.25)
.plot.rates(res.list$mse, theoretical.rate = function(x){log(x)*log(log(x))/x})

plot(NA, xlab = "log n", ylab = "n MSE", ylim = c(0, (median(res.list$mse[,10]) + sd(res.list$mse[,10])) * n.vec[10]),
 xlim = c(min(log(n.vec)), max(log(n.vec))), cex.lab = 1.25, cex.axis = 1.25)

tmpmat = res.list$mse
for(i in 1:ncol(tmpmat)) tmpmat[,i] = tmpmat[,i]*n.vec[i]
.plot.rates(tmpmat, n.vec = log(n.vec), theoretical.rate = function(x){x*log(x)})

dev.off()
