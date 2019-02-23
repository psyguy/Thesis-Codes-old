library(dplyr)
library("Hmisc")

make.random.graph <- function(size = 20, num.links = 100, distribution = "binary",
                              parameters = c(0,1), seed = 1){
  
  set.seed(seed)
  
  if(distribution == "binary") vals <- sample.int(size ^ 2, num.links + size)

  diagonals <- c(1:size)#cbind(1:size,1:size)
  
  edge.list <- cbind(vals %/% size + 0, vals %% size + 1)
  
  for(i in  1:(num.links+size)){
    if(edge.list[i,1]==edge.list[i,2]) edge.list[i,] <- NA
  }
  
  edge.list <- edge.list %>% na.omit()#)[1:num.links,]
  n.e <- dim(edge.list)[1]
  # edge.list <- edge.list[1:num.links,]
  # 
  # edge.list <- unique(edge.list[ , 1:2 ] )
  
  edge.list.rev <- edge.list[n.e:1,2] %>% cbind(edge.list[n.e:1,1])
  edge.list.symmetric <- edge.list %>% rbind(edge.list.rev)
  edge.list.symmetric <- unique(edge.list.symmetric[ , 1:2 ] )

  # l.e <- dim(edge.list.symmetric)[1]/2
  # 
  # while(num.links < l.e){
  #   
  #   
  #   
  # }
  # 
  # 
  # length.diff <- l.e - num.links
  # 
  # edge.list.final <- edge.list.symmetric[(1+length.diff):(l.e-length.diff),]
  
  adj.matrix <- matrix(0, size, size)
  adj.matrix[edge.list.symmetric] <- 1
  
  # adj.matrix %>% sum() %>% sum()
  # 
  # diag(adj.matrix)
  
  adj.matrix %>% return()
  
}
