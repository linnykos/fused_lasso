rm(list=ls())
library(genlasso)

dim1 = dim2 = 100
beta0 = matrix(0, dim1, dim2)

x = -(dim1/2):(dim1/2)
y = x^3-(dim1*3/8)^2*x
rescale.coef = max(y)/(dim1/2)

for(i in 1:dim1){
 
  x = i-dim1/2
  y = (x^3-(dim1*3/8)^2*x)/rescale.coef
  
  j = 1:dim2 
  idx = which(j - dim2/2 > y) 
  beta0[i,idx] = 1
}

set.seed(10)
y = beta0 + matrix(rnorm(dim1*dim2), dim1, dim2)

out = fusedlasso2d(y, dim1 = dim1, dim2 = dim2, maxsteps = 80000)
save.image(file = paste0("~/ryan/fused.git/results/graph-dim_",
 dim1, "_", Sys.Date(), ".RData"))
