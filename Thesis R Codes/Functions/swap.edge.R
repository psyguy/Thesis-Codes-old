library(dplyr)
swap.edge <- function(connectivity.matrix, i1, j1, i2, j2){
  
  edge.1 <- connectivity.matrix[i1,j1]
  edge.2 <- connectivity.matrix[i2,j2]
  
  connectivity.matrix[i1,j1] <- edge.2
  connectivity.matrix[j1,i1] <- edge.2

  connectivity.matrix[i2,j2] <- edge.1
  connectivity.matrix[j2,i2] <- edge.1
  
  connectivity.matrix %>% return()
  
}
