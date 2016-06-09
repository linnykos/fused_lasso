#specific to 2d graph
plot.2dgraph <- function(mat, est.edge.mat = NA, true.edge.mat = NA, 
 color.vec = gray.colors(30), zlim = c(min(as.numeric(mat)),
 max(as.numeric(mat)))) {

  m = nrow(mat)
  n = ncol(mat)

  image(mat, col = color.vec, zlim = zlim, asp = T, xlab = "",
   ylab = "", yaxt = "n", xaxt = "n", bty = "n")

  s1 = 1:(m-1)*1/(m-1)-1/((m-1)*2)
  s2 = 1:(n-1)*1/(n-1)-1/((n-1)*2)

  s1t = 0:m*1/(m-1)-1/((m-1)*2)
  s2t = 0:m*1/(m-1)-1/((m-1)*2)

  if(!all(is.na(est.edge.mat))){
    for(i in 1:nrow(est.edge.mat)){
      tmp = est.edge.mat[i,]
      if(tmp[1] != tmp[3]){
        lines(x = rep(s2[min(tmp[1],tmp[3])],2), 
         y = c(s1t[tmp[2]],s1t[tmp[2]+1]), col = "blue2", lwd = 2)
      } else {
        lines(x = c(s2t[tmp[1]],s1t[tmp[1]+1]), 
         y = rep(s1[min(tmp[2],tmp[4])],2), col = "blue2", lwd = 2)
      }
    }
  }

  if(!all(is.na(true.edge.mat))){
    for(i in 1:nrow(true.edge.mat)){
      tmp = true.edge.mat[i,]
      if(tmp[1] != tmp[3]){
        lines(x = rep(s2[min(tmp[1],tmp[3])],2),
         y = c(s1t[tmp[2]],s1t[tmp[2]+1]), col = 2, lwd = 2)
      } else {
        lines(x = c(s2t[tmp[1]],s1t[tmp[1]+1]), 
         y = rep(s1[min(tmp[2],tmp[4])],2), col = 2, lwd = 2)
      }
    }
  }

  invisible()
}

