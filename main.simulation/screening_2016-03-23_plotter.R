load("~/DUMP/screening_experiment2_2016-03-27.RData")

#make a 2-panel plot
lambda.mat = res$lambda.mat
n.vec = res$setup$n.vec


png(paste0("~/DUMP/screening-plot2_", DATE, ".png"), height = 5, width = 10, 
 units = "in", res = 300)

par(mfrow = c(1,2))

filter.opt.mean = apply(res$haus.filter.mat, 1, mean)
filter.opt.sd = apply(res$haus.filter.mat, 1, sd)
filter.est.opt.mean = apply(res$haus.filter.boot.mat, 1, mean)
filter.est.opt.sd = apply(res$haus.filter.boot.mat, 1, sd)
haar.opt.mean = apply(res$haus.haar.mat, 1, mean)
haar.opt.sd = apply(res$haus.haar.mat, 1, sd)

y.max = max(c(filter.opt.mean+filter.opt.sd, filter.est.opt.mean+filter.est.opt.sd, haar.opt.mean+haar.opt.sd))

#plot the hausdorff
plot(x = n.vec, y = filter.opt.mean, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylab = "Hausdorff Distance", ylim = c(0, y.max))
lines(x = n.vec, y = filter.opt.mean, lwd = 2)
for(i in 1:length(n.vec)) {
  lines(x = rep(n.vec[i], 2), y = c(filter.opt.mean[i] - filter.opt.sd[i], filter.opt.mean[i] + filter.opt.sd[i]),
   cex = 2)
}

points(x = n.vec, y = filter.est.opt.mean, col = "red", cex = 2, pch = 16, lty = 2)
lines(x = n.vec, y = filter.est.opt.mean, lwd = 2, col = "red", lty = 2)
for(i in 1:length(n.vec)) {
  lines(x = rep(n.vec[i], 2), y = c(filter.est.opt.mean[i] - filter.est.opt.sd[i], filter.est.opt.mean[i] + filter.est.opt.sd[i]),
   cex = 2, col = "red", lty = 2)
}

points(x = n.vec, y = haar.opt.mean, col = "blue", cex = 2, pch = 16, lty = 4)
lines(x = n.vec, y = haar.opt.mean, lwd = 2, col = "blue", lty = 4)
for(i in 1:length(n.vec)) {
  lines(x = rep(n.vec[i], 2), y = c(haar.opt.mean[i] - haar.opt.sd[i], haar.opt.mean[i] + haar.opt.sd[i]),
   cex = 2, col = "blue", lty = 4)
}

legend("topleft", c("Filter (Theoretical)", "Filter (Estimated)", "Haar (Theoretical)"), 
   lty = c(1,2,4), lwd = 2, bty = "n", col = c("black", "red", "blue"))


#second plot
par(mar = c(5,4,4,4))
plot(x = n.vec, y = filter.est.opt.mean, pch = 16, cex = 2, xlab = "Length of vector (n)",
     ylim = c(0, y.max), ylab = "", yaxt = 'n', col = "red")
lines(x = n.vec, y = filter.est.opt.mean, lwd = 2, col = "red", lty = 2)
for(i in 1:length(n.vec)) {
  lines(x = rep(n.vec[i], 2), y = c(filter.est.opt.mean[i] - filter.est.opt.sd[i], filter.est.opt.mean[i] + filter.est.opt.sd[i]),
   cex = 2, col = "red", lty = 2)
}
axis(2, ylim = c(0, y.max), col.axis = "red", col.ticks = "red")
mtext("Hausdorff Distance", side = 2, col = "red", line = 2.5)


level.mean = apply(res$level.mat, 1, mean)
level.sd = apply(res$level.mat, 1, sd)

par(new = TRUE)
plot(x = n.vec, y = level.mean, pch = 16, cex = 2, xlab = "", xaxt = 'n',
  ylim = c(0, max(level.mean + level.sd)), ylab = "", yaxt = 'n')
lines(x = n.vec, y = level.mean, lwd = 2)
for(i in 1:length(n.vec)) {
  lines(x = rep(n.vec[i], 2), y = c(level.mean[i] - level.sd[i], level.mean[i] + level.sd[i]),
   cex = 2)
}
axis(4, ylim = c(0,
  max(level.mean + level.sd)))
mtext("Estimated Filter Threshold", side = 4, line = 2.5)

dev.off()


