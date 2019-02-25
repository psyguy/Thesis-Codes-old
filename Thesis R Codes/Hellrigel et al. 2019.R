source("./My Thesis Functions/my.library.loader.R")

seed <- 1

set.seed(seed)
num.edges <- 7#5200
num.nodes <- 5#300
eps <- c(0.3, 0.4, 0.5)
a <- 1.7
t.transient <- 3#20 %>% as.vector()
num.epochs <- 10#10000
tol <- 0.001

conn <- make.random.graph(size = num.nodes, num.links = num.edges,
                          dist = "binary", seed = seed)


x.init <- num.nodes %>% runif(-1,1)
inp.x <- x.init

x <- x.init
connectivity.matrix <- conn
clustering.coef <- NULL

time.start <- Sys.time()
rewired <- NULL
for(iterations in 1:num.epochs){
  
  if(!(iterations %% 100)) paste("System time at iteration", iterations, "is", Sys.time()) %>% print()

  rewired <- my.rewiring(inp.x = x, inp.conn = connectivity.matrix,
                         t.transient = t.transient, eps = eps[3])
  x <- rewired$x
  (connectivity.matrix <- rewired$connectivity.matrix) %>% print()
  (clustering.coef[iterations] <- rewired$clustering.coef) %>% print()
}


time.end <- Sys.time()
time.taken <- time.end - time.start
paste("It took", time.taken, "to run for", num.epochs, "epochs.") %>% print()

current.status <- list(rewired = rewired, clustering.coef = clustering.coef, num.epochs = num.epochs, eps = eps[3], seed = 1)

save(current.status, file = "current.status.20180225.RData")
