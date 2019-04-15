library(dplyr)
rm(list=ls())
source("./My Thesis Functions/April.GongvLeeuwen2004.R")

seed <- 1
# number of links
L_c <- 5 #520
# number of nodes
N <- 5 #3000
# number of iterations
T_ <- 60 #600 
# tol <- 0.001

set.seed(seed)

conn <- make.random.graph(size = N, num.links = L_c, seed = seed)

x.init <- N %>% runif(-1,1)
x.out <- x.init %>% as.matrix() %>% t()

for (i in 1:T_) {
  x.lastrow <- x.out %>% tail(1)
  x.temp <- x.lastrow %>% GongvLeeuwen2004.logistic(conn)
  x.out <- x.temp %>% rbind(x.out)
}

# save(x.out, file = "x.out_5200.300.6000.Rdata")

distances <- x.out %>% tail(1) %>% GongvLeeuwen2004.coherenceD()#,connectivity.matrix = conn)

i_ <- sample.int(N,1)

d_ <- distances[,i_]
# d_ <- d_[i_]

j_1 <- which.min(d_)
j_2 <- which.max(d_)

if(!conn[i_,j_1]) conn <- conn %>% swap.edge(i1 = i_, j1 = j_1, i2 = i_, j2 = j_2)

g <- graph_from_adjacency_matrix(conn, mode = "undirected")

ClCoef <- transitivity(g)


transitivity
