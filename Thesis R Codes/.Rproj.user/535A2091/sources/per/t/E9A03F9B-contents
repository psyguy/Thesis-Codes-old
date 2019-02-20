library(dplyr)
library("Hmisc")

make.random.graph <- function(size = 20, num.links = 100){
  
  vals <- sample.int(size ^ 2, num.links + size)
  diagonals <- c(1:size)#cbind(1:size,1:size)
  
  edge.list <- cbind(vals %/% size + 0, vals %% size + 1)
  
  for(i in  1:(num.links+size)){
    if(edge.list[i,1]==edge.list[i,2]) edge.list[i,] <- NA
  }
  
  edge.list <- edge.list %>% na.omit()
  
  # edge.list <- edge.list[1:num.links,]
  # 
  # edge.list <- unique(edge.list[ , 1:2 ] )
  
  edge.list.rev <- edge.list[,2] %>% cbind(edge.list[,1])
  edge.list.symmetric <- edge.list %>% rbind(edge.list.rev)
  edge.list.symmetric <- unique(edge.list.symmetric[ , 1:2 ] )
  
  adj.matrix <- matrix(0, size, size)
  adj.matrix[edge.list.symmetric] <- 1
  
  # adj.matrix %>% sum() %>% sum()
  # 
  # diag(adj.matrix)
  
  adj.matrix %>% return()
  
}