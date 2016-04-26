rm(list=ls())
library("genlasso")
source("source.R")
source("source_plotter.R")
date = Sys.Date()
add.info = "Dmatrix"


sigma = 1
len = 20
jump.mean.org =     c(0, 2,  4, 1, 4)
jump.location.org= c(0, .2, .4, .6, .8)

y = generate.problem(len, jump.mean.org, jump.location.org,  sigma)

#first changepoint
D = form.diff.matrix(len)
Z = solve(D%*%t(D))%*%D
y2 = Z%*%y
y.sum = sum(y)

y3 = numeric(len)
for(i in 1:len){
  y3[i] = abs(sum(y[1:i])-i/len*y.sum)
}

dev.off()
res = dualplot_suite(y, jump.mean.org, jump.location.org, num.est.jumps=1)
idx = which(abs(diff(res$coef))>1e-5)

par(pty="s")
plot(y=cumsum(y),x=(1:len)/len*y.sum)
larg = 500
lines(x=c(-larg,larg), y=c(-larg,larg),col="red",lwd=2)
tmp = sum(y[1:idx])
tmp2 = idx/len*y.sum

lines(x=c(tmp2,tmp2), y=c(tmp,tmp2), col="blue", lwd=2)
abline(a=tmp2-tmp,b=1,col="blue",lwd=1)
abline(a=tmp-tmp2,b=1,col="blue",lwd=1)

lambda1 = max(Z%*%y)
s = sign(max(Z%*%y))

D2 = D[-idx,]
Z2 = solve(D2%*%t(D2))%*%D2
num = Z2%*%y
tmp = Z2%*%t(D[idx,,drop=FALSE])*s
denom = apply(cbind(tmp+1,tmp-1),1,max)
lambda2 = max(num/denom)