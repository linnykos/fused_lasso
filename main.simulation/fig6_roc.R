setwd("~/ryan/fused.git")
load("results/final-ROC2-2016-06-17.RData")
source("source_header.R")

n = n.vec[i]
filter.bandwidth = ceiling(0.25*(log(n))^2)

#left and right screening distance
haus.mat = vector("list", 4)
for(i in 1:4){
  haus.mat[[i]] = matrix(0, ncol = length(oracle.seq), nrow = trials)
}
names(haus.mat) = c("left.oracle.mat", "right.oracle.mat", "left.adapt.mat", "right.adapt.mat")

true.jumps = enumerate.jumps(truth)

for(i in 1:length(oracle.seq)){
  for(k in 1:trials){
    haus.mat[[1]][k,i] = compute.hausdorff(true.jumps, res.list[[2]][[i]][[k]], one.sided = T)
    haus.mat[[2]][k,i] = compute.hausdorff(res.list[[2]][[i]][[k]], true.jumps, one.sided = T)

    haus.mat[[3]][k,i] = compute.hausdorff(true.jumps, res.list[[3]][[i]][[k]], one.sided = T)
    haus.mat[[4]][k,i] = compute.hausdorff(res.list[[3]][[i]][[k]], true.jumps, one.sided = T)
  }
}

haus.vec = matrix(0, ncol = length(oracle.seq), nrow = 4)
for(i in 1:4){
  haus.vec[i,] = apply(haus.mat[[i]], 2, median, na.rm = T)
}

#classification quality
class.mat = vector("list", 4)
for(i in 1:4){
  class.mat[[i]] = matrix(0, ncol = length(oracle.seq), nrow = trials)
}
names(class.mat) = c("true.pos.oracle", "false.pos.oracle", "true.pos.adapt", "false.pos.oracle")

for(i in 1:length(oracle.seq)){
  for(k in 1:trials){
    res = classification.quality(true.jumps, res.list[[2]][[i]][[k]], filter.bandwidth, n) 

    class.mat[[1]][k,i] = res$true.pos
    class.mat[[2]][k,i] = res$false.pos

    res = classification.quality(true.jumps, res.list[[3]][[i]][[k]], filter.bandwidth, n)

    class.mat[[3]][k,i] = res$true.pos
    class.mat[[4]][k,i] = res$false.pos
  }
}

class.vec = matrix(0, ncol = length(oracle.seq), nrow = 4)
for(i in 1:4){
  class.vec[i,] = apply(class.mat[[i]], 2, mean, na.rm = T)
}

pdf(paste0("plots/ROC_", Sys.Date(), ".pdf"), height = 3, width = 9)
par(mar = c(4.25,4.25,1,1))
par(mfrow = c(1,3))

#plot 1: of hausdorff distance
plot(x = haus.vec[1,], y = haus.vec[2,], pch = 16,
  xlab = "Screening distance", 
  ylab = "Precision distance", cex.axis = 1.25, cex.lab = 1.25)
lines(x = haus.vec[1,], y = haus.vec[2,], lwd = 2)

#lines(x = haus.vec[3,], y = haus.vec[4,], lwd = 2, lty = 2, col = 2)

idx = which(colnames(filter.mat) == "0.95")
points(x = haus.vec[3,idx], y = haus.vec[4,idx], pch = 16, cex = 2, col = 2)

#plot 2: of classification
plot(x = class.vec[2,], y = class.vec[1,], pch = 16,
  ylab = "TPR", 
  xlab = "FPR", cex.axis = 1.25,
  cex.lab = 1.25)
lines(x = class.vec[2,], y = class.vec[1,], lwd = 2)

#lines(x = class.vec[4,], y = class.vec[3,], lwd = 2, lty = 2, col = 2)

idx = which(colnames(filter.mat) == "0.95")
points(x = class.vec[4,idx], y = class.vec[3,idx], pch = 16, cex = 2, col = 2)

#plot 3: of the proportion
plot(x = 1 - as.numeric(colnames(filter.mat)), y = class.vec[4,], pch = 16,  
 xlab = "1 - q", ylab = "FDR",
 cex.lab = 1.25, cex.axis = 1.25)
lines(x = 1 - as.numeric(colnames(filter.mat)), y = class.vec[4,], lwd = 2)
lines(x = c(0,1), y = c(0,1), lty = 2, lwd = 2)

#lines(x = rep(0.95, 2), y = c(-1,2), lty = 2, lwd = 2, col = 2)

graphics.off()
quit()

