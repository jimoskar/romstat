### Problem 2 ----
library(spatial)
library(tidyverse)

## a) ----

a.mat = read.table("obsprob.txt", header = TRUE)
a.vec = a.mat$alpha
p.mat = read.table("obspines.txt", header = TRUE)
p.vec = p.mat$N_obs
p.mat

x.center <- y.center <- (1:30)*10 - 5 
grid <- expand.grid(y.center, x.center)

ggplot(data = a.mat, aes(x, y)) +
  geom_tile(aes(fill = alpha)) +
  scale_fill_viridis_c(name= expression(alpha)) +
  coord_fixed() + xlab("x") + ylab("y") + 
  ggtitle("Detection probability") + theme_minimal()

ggplot(p.mat, aes(x, y)) +
  geom_tile(aes(fill = N_obs)) +
  scale_fill_viridis(expression(N[obs])) +
  coord_fixed() + xlab("x") + ylab("y") + 
  ggtitle("Detected pine counts") + theme_minimal()

## c) ----
# Calculate estimator
C <- 1/(100 * sum(a.vec))
Lambda2 <- C*sum(p.vec)
Lambda2
# Simulate counts within each cell
sim.points <- function(lambda, plot=TRUE){
  dx <- dy <- 10 # width and height of the cells
  lambda.cell <- 100*lambda
  N.mat <- data.frame(x = p.mat$x, y = p.mat$y)
  N.vec <- rpois(900, lambda.cell)
  N.mat$N <- N.vec
  p <- ggplot()
  for(i in 1:nrow(N.mat)){
    n <- N.mat$N[i] # number of points in current cell
    x <- runif(n, N.mat$x[i] - dx/2, N.mat$x[i] + dx/2)
    y <- runif(n, N.mat$y[i] - dy/2, N.mat$y[i] + dy/2)
    p <- p + geom_point(data = data.frame(x = x, y = y), aes(x,y))
  }
  if(plot == TRUE) p + theme_minimal()
}
# Generate 3 realizations
set.seed(2)
sim.points(Lambda2)
sim.points(Lambda2)
sim.points(Lambda2)
