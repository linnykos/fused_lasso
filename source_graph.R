#specific to 2d graphs
enumerate.jumps.2dgraph <- function(mat){
  m = nrow(mat)
  n = ncol(mat)
  
  edge.mat = matrix(NA,1,4)
  
  #populate edge.mat
  for(i in 1:m){
    for(j in 1:(n-1)){
      if(mat[i,j]!=mat[i,j+1]) edge.mat = rbind(edge.mat,c(i,j,i,j+1))
    }
  }
  
  for(j in 1:n){
    for(i in 1:(m-1)){
      if(mat[i,j]!=mat[i+1,j]) edge.mat = rbind(edge.mat,c(i,j,i+1,j))
    }
  }
  
  #remove the first edge
  edge.mat = edge.mat[-1,]
  
  edge.mat
}

#m is num of pixels along x, n along y
generate.2dgraph <- function(m,n){
  edges = matrix(NA,1,2)
  
  for(j in 1:n){
    for(i in 1:(m-1)){
      edges = rbind(edges,c((j-1)*m+i,(j-1)*m+i+1))
    }
  }
  
  for(i in 1:m){
    for(j in 1:(n-1)){
      edges = rbind(edges,c(i+(j-1)*m,i+j*m))
    }
  }
  
  edges = edges[-1,]
  
  ggraph = add.edges(graph.empty(m*n), as.numeric(t(edges)))
  
  ggraph
}

#specific to 2d graphs
.convert.edge2graph <- function(edge.mat, m, n){
  if(nrow(edge.mat)==0) return()
  
  conversion <- function(vec){
    if(vec[1] != vec[3]){
      tmp = min(vec[1],vec[3])
      tmp1 = c((vec[2]-1)*m+tmp,(vec[2]-1)*m+tmp+1)
    } else {
      tmp = min(vec[2],vec[4])
      tmp1 = c((tmp-1)*m+vec[1],(tmp)*m+vec[1])
    }
    
    tmp1
  }
  
  node.vec = apply(edge.mat,1,conversion)

  t(node.vec)
}

#for any set of nodes
screening_dist <- function(ggraph, base.node.vec, cover.node.vec){
  if(length(base.node.vec)==0) return(0)
  if(length(cover.node.vec)==0) return(NA)
  
  base.vec = unique(as.numeric(base.node.vec))
  cover.vec = unique(as.numeric(cover.node.vec))
  
  dist = 0
  shortpath = shortest.paths(ggraph)
  
  for(i in 1:nrow(base.node.vec)){
    tmp1 = min(shortpath[base.node.vec[i,1],cover.vec])
    tmp2 = min(shortpath[base.node.vec[i,2],cover.vec])
    
    tmp = min(tmp1,tmp2)
    #special case only when tmp = 0
    if(tmp1 == tmp2){
      if(tmp1 == 0){
        idx = which(cover.node.vec==base.node.vec[i,1],arr.ind=TRUE)[,1]
        if(length(which(cover.node.vec[idx,]==base.node.vec[i,2]))==0) tmp = tmp+1
      } else {
        tmp = tmp+1
      }
    } else {
      if(tmp1 == 0 | tmp2 == 0) tmp = 1
    }
    
    #print(paste("Row ",i," got distance of ",tmp1," and ",tmp2,sep=""))
    
    if(tmp>dist) dist = tmp
  }
  
  dist
}

