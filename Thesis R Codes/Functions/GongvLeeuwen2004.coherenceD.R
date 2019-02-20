library(dplyr)

GongvLeeuwen2004.coherenceD <- function(x.input, connectivity.matrix){

  inp.dim <- x.input %>% dim()
  hight <- inp.dim[1]
  leng <- inp.dim[2]
  
  x.now <- x.input[,leng]
  
  x.duplicated <- x.now
  for (i in 1:(hight-1)){x.duplicated <- x.now %>% cbind(x.duplicated)}
  
  x.duplicated <- x.duplicated %>% as.matrix()
  
  d <- (t(x.duplicated) - x.duplicated) %>% abs()
  
  d %>% return()
  
}
