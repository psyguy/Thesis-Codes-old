source("./My Thesis Functions/my.library.loader.R")

seed <- 1

set.seed(seed)
num.edges <- 5200
num.nodes <- 300
eps <- c(0.3, 0.4, 0.5)
a <- 1.7
t.transient <- 20

tol <- 0.001

conn <- make.random.graph(size = num.nodes, num.links = num.edges,
                          dist = "binary", seed = seed)


x.init <- num.nodes %>% runif(-1,1)
x.out <- x.init

# This for loop keeps the history of the x's, which would be better if implemented separately.
for (i in 1:t.transient) {
  x.temp <- x.out %>% Hellrigel2019.logistic(conn)
  x.out <- x.out %>% cbind(x.temp)
}


# save(x.out, file = "x.out_5200.300.6000.Rdata")

distances <- Hellrigel2019.syncerror(x.out,connectivity.matrix = conn)

i_ <- sample.int(num.nodes,1)

d_ <- distances[,i_]

k_ <- which.min(d_)



l_ <- (d_ * x.out[i_,]) %>% which.max()

if(!conn[i_,j_1]) conn <- conn %>% swap.edge()

# Skipped step 5 probabilistic approach

conn <- conn %>% swap.edge()

g <- graph_from_adjacency_matrix(conn, mode = "undirected")

