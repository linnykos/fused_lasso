setwd("~/ryan/fused.git/")

load("results/final-ROC2-2016-05-21.RData")
source("source_header.R")

left.mat = matrix(0, ncol = length(oracle.seq), nrow = trials)
colnames(left.mat) = oracle.seq
right.mat = left.mat

true.jumps = enumerate.jumps(truth)

for(i in 1:length(oracle.seq)){
  for(k in 1:trials){
    left.mat[k,i] = compute.hausdorff(true.jumps, res.list[[2]][[i]][[k]], one.sided = T)
    right.mat[k,i] = compute.hausdorff(res.list[[2]][[i]][[k]], true.jumps, one.sided = T)
  }
}

pdf(paste0("plots/Haus_", Sys.Date(), ".pdf"), height = 5, width = 5)

plot.hausdorff(left.mat, right.mat)

dev.off()
