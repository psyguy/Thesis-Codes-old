my.rewiring <- function(inp.x, inp.conn, a = 1.7, eps = 0.8){
  
  x.temp <- inp.x
  conn <- inp.conn
  
  for (i in 1:t.transient) x.temp <- x.temp %>% Hellrigel2019.logistic(conn, a = a, eps = eps)
  
  distances <- Hellrigel2019.syncerror(inp.x, connectivity.matrix = conn)
  
  num.nodes <- dim(inp.conn)[1]
  
  i_ <- sample.int(num.nodes,1)
  d_ <- distances[,i_]
  k_ <- which.min(d_)
  l_ <- (d_ * inp.x[i_]) %>% which.max()
  
  if(!conn[i_,k_]) conn <- conn %>% swap.edge()

  # Skipped step 5 probabilistic approach
  
  g <- igraph::graph_from_adjacency_matrix(conn, mode = "undirected")
  clustering.coef <- transitivity(g)
  
  output <- list(x = x.temp, connectivity.matrix = conn, clustering.coef = clustering.coef)
  
  output %>% return()

}
