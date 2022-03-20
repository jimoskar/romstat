### Problem 2 ----
library(spatial)
library(tidyverse)
library(viridis)

## a) ----

a.mat = read.table("obsprob.txt", header = TRUE)
a.vec = a.mat$alpha
p.mat = read.table("obspines.txt", header = TRUE)
p.vec = p.mat$N_obs

x.center <- y.center <- (1:30)*10 - 5 
grid <- expand.grid(y.center, x.center)

# Plot observation probs
ggplot(data = a.mat, aes(x, y)) +
  geom_tile(aes(fill = alpha)) +
  scale_fill_viridis(name= expression(alpha)) +
  coord_fixed() + xlab("x") + ylab("y") + 
  ggtitle("Detection probability") + theme_minimal()

# Plot counts
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
# Function for plotting point pattern or just generate realization
sim.points <- function(lambda, p.mat, plot.points=TRUE){
  dx <- dy <- 10 # width and height of the cells
  lambda.cell <- 100*lambda
  N.mat <- data.frame(x = p.mat$x, y = p.mat$y)
  N.vec <- rpois(900, lambda.cell)
  N.mat$N <- N.vec
  if(plot.points){
    p <- ggplot()
    for(i in 1:nrow(N.mat)){
      n <- N.mat$N[i] # number of points in current cell
      x <- runif(n, N.mat$x[i] - dx/2, N.mat$x[i] + dx/2)
      y <- runif(n, N.mat$y[i] - dy/2, N.mat$y[i] + dy/2)
      p <- p + geom_point(data = data.frame(x = x, y = y), aes(x,y))
    }
    p + theme_minimal()
  }
  else{
    return(N.mat)
  }
}
# Generate 3 realizations
set.seed(2)
sim.points(Lambda2, p.mat)
sim.points(Lambda2, p.mat)
sim.points(Lambda2, p.mat)

## d) ----
# Function for plotting point pattern or generate realization
sim.cond <- function(lambda, p.mat, a.mat, plot.points=TRUE){
  a.vec <- a.mat$alpha
  dx <- dy <- 10 # width and height of the cells
  lambda.undetected <- 100*lambda*(1 - a.vec)
  undetected.vec <- rpois(900, lambda = lambda.undetected) 
  N.mat <- data.frame(x = a.mat$x, y = a.mat$y)
  N.mat$N <- undetected.vec + p.mat$N_obs
  if(plot.points){
    p <- ggplot()
    for(i in 1:nrow(N.mat)){
      n <- N.mat$N[i] # number of points in current cell
      x <- runif(n, N.mat$x[i] - dx/2, N.mat$x[i] + dx/2)
      y <- runif(n, N.mat$y[i] - dy/2, N.mat$y[i] + dy/2)
      p <- p + geom_point(data = data.frame(x = x, y = y), aes(x,y))
    }
    p + theme_minimal()
  }
  else{
    return(N.mat)
  }
}

# Create 3 realizations:
set.seed(2)
sim.cond(Lambda2, p.mat, a.mat)
sim.cond(Lambda2, p.mat, a.mat)
sim.cond(Lambda2, p.mat, a.mat)

## e) ----
# Calculate the means of the conditional  N|M and marginal N over 500 realizations:
N.mat <- NM.mat <- matrix(rep(NA,500*900), nrow = 500)
for(i in 1:500){
  N.mat[i, ] = sim.points(Lambda2, p.mat, plot.point=FALSE)$N
  NM.mat[i, ] = sim.cond(Lambda2, p.mat, a.mat, plot.point=FALSE)$N
}
N.mean <- apply(N.mat, 2, mean)
NM.mean <- apply(NM.mat, 2, mean)
means.df <- data.frame(mean = c(N.mean, NM.mean), 
                       x = rep(p.mat$x, 2), y = rep(p.mat$y, 2),
                       id = c(rep("widehat(E(N))", 900), rep("widehat(E(N*'|'*M == m))", 900)))

ggplot(means.df, aes(x,y)) +
  geom_tile(aes(fill = mean)) + facet_grid(~id, labeller = label_parsed) + 
  scale_fill_viridis(" ") +
  coord_fixed() + xlab("x") + ylab("y")  + theme_minimal()

# Do the same for standard deviation:
N.sd <- apply(N.mat, 2, sd)
NM.sd <- apply(NM.mat, 2, sd)
sd.df <- data.frame(sd = c(N.sd, NM.sd), 
                       x = rep(p.mat$x, 2), y = rep(p.mat$y, 2),
                       id = c(rep("widehat(SD(N))", 900), rep("widehat(SD(N*'|'*M == m))", 900)))

ggplot(sd.df, aes(x,y)) +
  geom_tile(aes(fill = sd)) + facet_grid(~id, labeller = label_parsed) + 
  scale_fill_viridis(" ") +
  coord_fixed() + xlab("x") + ylab("y")  + theme_minimal()
