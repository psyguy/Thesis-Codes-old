library(dplyr)

func.1 <- function(x, a = 1.7){1 - a*(x^2)}

make.random.graph <- function(size = 5,
                              num.links = 10,
                              distribution = "binary",
                              parameters = c(0,1), seed = 1){
  set.seed(seed)
  edges <- rep(0,size*(size-1)/2)
  
  if(distribution == "binary"){
    ones <- sample.int(length(edges), num.links)
    edges[ones] <- 1
  }
  
  g <- matrix(0, size, size)
  g[base::lower.tri(g, diag=FALSE)] <- edges
  
  (adj.matrix <- g + t(g)) %>% return()
}

GongvLeeuwen2004.logistic <- function(x.input,
                                      connectivity.matrix,
                                      a = 1.7, eps = 0.8,
                                      order = 0){
  # eps: coupling strength
  # x.input <- x.lastrow
  #   t(x.out) %>% rbind(t(x.out))
  # # connectivity.matrix <- conn
  # 
  
  x.now <- x.input # %>% t()
  # if(dim(x.input)[1]=1) x.now <- x.input %>% tail(1)
  # if(is.null(dim(x.input))) x.now <- x.input %>% t()# as.matrix() %>% t()
  
  n.nodes <- x.now %>% ncol()
  
  # x.now <- x.input %>% t()
  
  # Higher orders are not implemented yet
  #if(order) x <- x.input[,len-order:len]
  
  # unit.vector allows to calculate M_i by multiplying it the connectivity matrix
  # unit.vector <- matrix(1, height, 1)
  
  M <- connectivity.matrix %*% rep(1,n.nodes) %>% t()
  fx <- x.now %>% func.1() %>% as.matrix()# %>% t()
  
  x.next <- (1 - eps)*fx +  fx %*% connectivity.matrix*eps / M
  
  x.next %>% return()
  
}

func.1 <- function(x, a = 1.7){1 - a*(x^2)}