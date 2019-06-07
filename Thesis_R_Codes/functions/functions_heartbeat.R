# Intro -------------------------------------------------------------------

## This script keeps all functions needed for GongvLeeuwen2004.
## It's standalone and loads functions_my automatically

# rm(list = ls())


# loading the basic functions and packages --------------------------------

if(!exists("path.to.functions_my")) path.to.functions_my <- "functions"
source(paste0(path.to.functions_my,"/functions_my.R"))


# logistic function -------------------------------------------------------

heartbeat_logistic <- function(activities_,
                               m, #connectiviti matrix
                               a = 1.7,
                               eps = 0.8) {
  # eps: coupling strength
  # eps <- 0.8
  # x.input <- x.init
  # connectivity.matrix <- conn
  
  activities_now <- activities_
  num_nodes <- activities_now %>% ncol()
  
  
  M <- m %*% rep(1, num_nodes) %>% t()
  fx <- activities_now %>% mini_logistic() %>% as.matrix()
  
  neighbors <- fx %*% m * eps / M
  neighbors[is.nan(neighbors)] <- 0
  
  activities_next <- (1 - eps) * fx + neighbors
  
  activities_next %>% return()
  
}


# making a S4 class -------------------------------------------------------

brain <- setClass("brain",
                    slots = list(birthday = "character",
                                 age = "list",
                                 starting_values = "list", # activities_init
                                 parameters = "list", # num_nodes, num_edges
                                 history = "list", # activities_, coef.clustering, connmatrix
                                 now = "list" # activities, coef.clustering, mat.connectivity
                                 ))


# functions fir brain class --------------------------------------------------

setMethod("brain", signature(object = "brain"), function(object) {
  object@sides
})

setMethod("update",
          "brain",
          function(object) {
            beat <- object@age$beat
            m <- object@now$mat.connectivity
            a <- object@now$activities
            c <- object@now$coef.clustering
            
            a_next <- a %>% GongvLeeuwen2004_logistic(m)
            
            object@history$activities_all <-
              object@history$activities_all %>%
              rbind(a_next)
            
            object@now$activities <- a_next
            object@age$beat <- beat + 1
            # c_next <- m %>% my_clustceof()
          })

setMethod("rewiring",
          # signature(object = "brain"),
          "brain",
          function(object) {
            minute <- object@age$minute
            m <- object@now$mat.connectivity
            a <- object@now$activities
            c <- object@now$coef.clustering
            
            m_next <- my_rewire(a, m)
            
            object@history$mat.connectivity[beat] <- m_next
            
            object@now$mat.connectivity <- m_next
            object@age$minute <- minute + 1
            # c_next <- m %>% my_clustceof()
          })


x.start <- input$x.tot %>% tail(1)
conn.start <- input$connectivity.matrix

x.end <- x.start %>% GongvLeeuwen2004_logistic(conn.start)
conn.end <- my_rewire(x.end, conn.start)

g <-
  igraph::graph_from_adjacency_matrix(conn.end, mode = "undirected")
ClCoef <- igraph::transitivity(g)

output <- list(
  x.tot = rbind(input$x.tot, x.end),
  clustering.coefficient = c(input$clustering.coefficient, ClCoef),
  connectivity.matrix = conn.end
)



# GongvLeeuwen2004 adaptive rewiring  ------------------------------------------------------

heartbeat_rewire <- function(input) {
  x.start <- input$x.tot %>% tail(1)
  conn.start <- input$connectivity.matrix
  
  x.end <- x.start %>% GongvLeeuwen2004_logistic(conn.start)
  conn.end <- my_rewire(x.end, conn.start)
  
  g <-
    igraph::graph_from_adjacency_matrix(conn.end, mode = "undirected")
  ClCoef <- igraph::transitivity(g)
  
  output <- list(
    x.tot = rbind(input$x.tot, x.end),
    clustering.coefficient = c(input$clustering.coefficient, ClCoef),
    connectivity.matrix = conn.end
  )
  class(output) <- append(class(output), "GongvLeeuwen2004")
  output %>% return()
}
