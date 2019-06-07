a <- matrix(
  c(
    0,1,0,0,0,0,0,
    1,0,1,0,1,0,0,
    0,1,0,1,1,0,0,
    0,0,1,0,0,1,0,
    0,1,1,0,0,1,0,
    0,0,0,1,1,0,1,
    0,0,0,0,0,1,0
  ),
  nrow = 7
)

tr(a)

a %*% a %*% a %>% tr()


a %>% graph_from_adjacency_matrix() %>% transitivity()

a %*% a %>% tr()

a



rm(list = ls())
source("./functions/functions_GongvLeeuwen2004.R")
seed <- 1
# number of nodes
num_nodes <- 30
# number of links
num_edges <- round(2 * log(num_nodes) * (num_nodes - 1))
# number of iterations
T_ <- 10
# tol <- 0.001


# start of the script -----------------------------------------------------

set.seed(seed)
conn <- make_random_graph(size = num_nodes,
                          num.links = num_edges,
                          seed = rnorm(1))
conn %>% my_clustceof()

conn %>% pimage()
o <- conn %>% seriate()
pimage(conn,o)
