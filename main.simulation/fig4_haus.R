setwd("~/ryan/fused.git")
load("results/final-2016-06-15-0.95.RData")

compute.haus <- function(mat1, mat2){
  haus = pmax(mat1, mat2)
  med = apply(haus, 2, median, na.rm = T)
  upper = apply(haus, 2, quantile, prob = 0.75, na.rm = T)
  lower = apply(haus, 2, quantile, prob = 0.25, na.rm = T)

  list(median = med, upper = upper, lower = lower)
}

haus.list = vector("list", 3)
haus.list[[1]] = compute.haus(res.list$left.org, res.list$right.org)
haus.list[[2]] = compute.haus(res.list$left.oracle, res.list$right.oracle)
haus.list[[3]] = compute.haus(res.list$left.adapt, res.list$right.adapt)

idx = which( colnames(res.list$threshold.adapt[[1]]) == "0.95")
thres.list = vector("list", 2)
thres.list.upper = vector("list", 2)
thres.list.lower = vector("list", 2)

thres.list[[1]] = apply(res.list$threshold.oracle, 2, median, na.rm = T)
thres.list.upper[[1]] = apply(res.list$threshold.oracle, 2, quantile, 
  na.rm = T, prob = 0.75)
thres.list.lower[[1]] = apply(res.list$threshold.oracle, 2, quantile,
  na.rm = T, prob = 0.25)


thres.list[[2]] = sapply(1:length(res.list$threshold.adapt), function(x){
  median(res.list$threshold.adapt[[x]][,idx], na.rm = T)
})
thres.list.upper[[2]] = sapply(1:length(res.list$threshold.adapt), function(x){
  quantile(res.list$threshold.adapt[[x]][,idx], na.rm = T, prob = 0.75)
})
thres.list.lower[[2]] = sapply(1:length(res.list$threshold.adapt), function(x){
  quantile(res.list$threshold.adapt[[x]][,idx], na.rm = T, prob = 0.25)
})


pdf(file = paste0("plots/haus-dist-", Sys.Date(), ".pdf"), width = 9,
 height = 3)

par(mar = c(4.25,4.25,1,1))
par(mfrow = c(1,3))

#plot the hausdorff
plot(NA, xlim = c(0, max(n.vec)), ylim = c(0, max(haus.list[[1]]$upper)),
  xlab = "n", ylab = "Hausdorff distance", cex.axis = 1.25, cex.lab = 1.25)

col.vec = c(1, 2, 4)
lty.vec = c(3, 1, 2)

for(i in c(1,3,2)){
  points(x = n.vec, y = haus.list[[i]]$med, pch = 16, col = col.vec[i],
   cex = 2)
  lines(x = n.vec, y = haus.list[[i]]$med, lty = lty.vec[i], col = col.vec[i],
   lwd = 3)
  
  for(j in 1:length(n.vec)){
    lines(x = rep(n.vec[j], 2), y = c(haus.list[[i]]$lower[j], 
     haus.list[[i]]$upper[j]), lty = lty.vec[i], col = col.vec[i], lwd = 2)
  }
}

legend("topleft", c("Original", "Data-driven", "Oracle"),
 col = col.vec[c(1,3,2)], lty = lty.vec[c(1,3,2)], 
 lwd = 2, bty = "n", cex = 1.25)


#plot the hausdorff zoom-up
plot(NA, xlim = c(0, max(n.vec)), ylim = c(0, max(haus.list[[3]]$upper)),
  xlab = "n", ylab = "Hausdorff distance", cex.axis = 1.25, cex.lab = 1.25)

for(i in 3:2){
  points(x = n.vec, y = haus.list[[i]]$med, pch = 16, col = col.vec[i],
   cex = 2)
  lines(x = n.vec, y = haus.list[[i]]$med, lty = lty.vec[i], col = col.vec[i],
   lwd = 3)

  for(j in 1:length(n.vec)){
    lines(x = rep(n.vec[j], 2), y = c(haus.list[[i]]$lower[j],
     haus.list[[i]]$upper[j]), lty = lty.vec[i], col = col.vec[i], lwd = 2)
  }
}

legend("topright", c("Data-driven", "Oracle"),
 col = col.vec[3:2], lty = lty.vec[3:2], lwd = 2, bty = "n", cex = 1.25)


#plot the threshold values
col.vec = c(2, 4)
lty.vec = c(1,2)

plot(NA, xlim = c(0, max(n.vec)), ylim = c(min(unlist(thres.list.lower)), 
  max(unlist(thres.list.upper))+.2),
  xlab = "n", ylab = "Threshold", cex.lab = 1.25, cex.axis = 1.25)

for(i in 2:1){
  points(x = n.vec, y = thres.list[[i]], pch = 16, col = col.vec[i], cex = 2)
  lines(x = n.vec, y = thres.list[[i]], lwd = 3, col = col.vec[i], 
   lty = lty.vec[i])

  for(j in 1:length(n.vec)){
    lines(x = rep(n.vec[j], 2), y = c(thres.list.lower[[i]][j],
     thres.list.upper[[i]][j]), lwd = 2, col = col.vec[i], lty = lty.vec[i])
  }
}

legend("topright", c("Data-driven", "Oracle"),
 col = col.vec[2:1], lwd = 2, lty = c(2,1), bty = "n", cex = 1.25)

graphics.off()
quit()
