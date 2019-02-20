library(dplyr)

L_c <- 5200
N <- 300
T_ <- 6000 

conn <- make.random.graph(size = N, num.links = L_c)

x.init <- N %>% runif(-1,1)
x.out <- x.init

for (i in 1:T_) {
  x.temp <- x.out %>% GongvLeeuwen2004.logistic(conn)
  x.out <- x.out %>% cbind(x.temp)
}

