plot.hausdorff <- function(mat1, mat2){
  n = ncol(mat1)
  assert_that(ncol(mat2) == n)

  mean1 = apply(mat1, 2, mean, na.rm = T)
  mean2 = apply(mat2, 2, mean, na.rm = T)
  
  sd1 = apply(mat1, 2, sd, na.rm = T)
  sd2 = apply(mat2, 2, sd, na.rm = T)

  top1 = mean1 + sd1
  bot1 = mean1 - sd1
  top2 = mean2 + sd2
  bot2 = mean2 - sd2
  
  #lower bound by 0
  bot1 = sapply(bot1, function(x){max(x,0)})
  bot2 = sapply(bot2, function(x){max(x,0)})

  x.axs = as.numeric(colnames(mat1))
  ylim1 = c(min(bot1, na.rm = T), max(top1, na.rm = T))
  ylim2 = c(min(bot2, na.rm = T), max(top2, na.rm = T))

  par(mar = c(5,4,4,4))
  plot(x = x.axs, y = mean2, ylim = ylim2, pch = 16, cex = 1.5, 
    col = "red", xlab = "Filter Threshold Value",
    ylab = "", yaxt = 'n')
  for(i in 1:n){
    lines(x = rep(x.axs[i], 2), y = c(bot2[i], top2[i]), lwd = 1, lty = 1,
    col = "red")
  }
  axis(2, ylim = ylim2, col.axis = "red", col.ticks = "red")
  mtext("Distance from Estimate to Truth", side = 2, col = "red", line = 2.5)

  par(new = T)

  plot(x = x.axs, y = mean1, ylim = ylim1, pch = 16, cex = 1.5,
   xlab = "", xaxt = 'n', ylab = "", yaxt = 'n')
  for(i in 1:n){
    lines(x = rep(x.axs[i], 2), y = c(bot1[i], top1[i]), lwd = 1, lty = 1)
  }
  axis(4, ylim = ylim1)
  mtext("Distance from Truth to Estimate", side = 4, line = 2.5)	

  invisible()
}

plot.roc <- function(res, n.level = 10){
  oracle.1 = apply(res$oracle.left.list[[n.level]], 2, mean, na.rm = T)
  oracle.2 = apply(res$oracle.right.list[[n.level]], 2, mean, na.rm = T) 

  #remove any indices with NA's
  idx = unique(c(which(is.na(oracle.1)), which(is.na(oracle.2))))
  oracle.1 = oracle.1[-idx]
  oracle.2 = oracle.2[-idx]

  #keep track of the normalizing constants
  max1 = max(oracle.1)
  max2 = max(oracle.2)
  oracle.1 = oracle.1/max1
  oracle.2 = oracle.2/max2

  #plot the ROC curve
  par(mar = c(5,6,2,2))
  plot(x = oracle.1, y = 1 - oracle.2, pch = 16, xlab = "Distance from
   Truth to Estimate (Normalized)", ylab = "Distance from Estimate to
   Truth (Normalized)", yaxt = 'n')
  axis(2, seq(0, 1, by = 0.2), seq(1, 0, by = -0.2))

  #now compute the same for the actual estimated methods
  j = 1
  meth.1 = apply(res$left.list[[j]][[n.level]], 2, mean, na.rm = T)
  meth.2 = apply(res$right.list[[j]][[n.level]], 2, mean, na.rm = T)
  meth.1 = meth.1/max1
  meth.2 = meth.2/max2

  points(x = meth.1, y = 1 - meth.2, pch = 16, col = "red")

}


#WARNING: MESSY CODE!!
library(assertthat)
load("~/ryan/fused.git/results/filterExperiment_hausdorff_2016-05-04.RData")
mat1 = res$oracle.left.list[[10]]
mat2 = res$oracle.right.list[[10]]

png(paste0("~/DUMP/Haus_", Sys.Date(), ".png"), height = 5, width = 5, 
 units = "in", res = 300)
plot.hausdorff(mat1, mat2)
dev.off()


haus.res = res
load("~/ryan/fused.git/results/filterExperiment_2016-05-02.RData")
level.res = res
