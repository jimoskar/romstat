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
