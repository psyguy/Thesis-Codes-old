# Notes -------------------------------------------------------------------

#
#   
#   at each heartbeat nodes are updated (logistic functions)
#   and every minute the rewiring takes place
#
#
#
#
#
#


# init --------------------------------------------------------------------

rm(list = ls())
source("./functions/functions_heartbeat.R")
seed <- 1

# number of nodes
num_nodes <- 300
# number of links
num_edges <- round(2 * log(num_nodes) * (num_nodes - 1))

# number of heartbeats per minute
num_hr <- 10
# number of rewirings
num_minutes <- 500


# start of the script -----------------------------------------------------

set.seed(seed)

# making a random matrix
mat_adj <- make_random_graph(size = num_nodes,
                          num.links = num_edges,
                          seed = seed)

# making initial values of the node activations
activities_init <- num_nodes %>% runif(-1, 1) %>% t()

# Cl
c <-  mat_adj %>% my_clustceof()


GvL.start <- list(
  x.tot = x.init,
  clustering.coefficient = ClCoef,
  connectivity.matrix = conn
)

GvL.finished <- GvL.start

time_start <- Sys.time()
# T_ <- 4000 #since more is redundant
for (i in 1:T_) {
  print(i)
  GvL.finished <- GvL.finished %>%
    GongvLeeuwen2004_adaptive_rewire()
  g <- GvL.finished$connectivity.matrix %>%
    igraph::graph_from_adjacency_matrix(mode = "undirected")
  (ClCoef <- igraph::transitivity(g)) %>% print()
  if (!(i %% (T_ / 40))) {
    title <- paste0("t = ", i, ", Clustering Coefficient = ", ClCoef)
    d <- GvL.finished$connectivity.matrix
    o <- seriate(d)
    pimage(d, o, main = title)
    # GvL.finished$connectivity.matrix %>% corrplot(method = "square")
  }
  
}


time_taken <- Sys.time - time_start
GvL.finished$clustering.coefficient %>% plot()

save_vars(ls())

# GvL.finished$connectivity.matrix %>% corrplot(method = "square")

# GvL.finished$x.tot[, ]
#
# (GvL.finished$x.tot[2, ] - GvL.finished$x.tot[3, ]) %>% sum()


# for(k in c(1:11)){

# densityplot(GvL.finished$x.tot[50,])
#
# GvL.finished$x.tot %>% corrplot()




