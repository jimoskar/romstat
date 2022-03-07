### Problem 1 ----
library(spatial)
library(MASS)
library(tidyverse)

## a) ----
# Plot cells point data.
cells <- ppinit("cells.dat")
cells.df <- data.frame(x = cells$x, y = cells$y)
cells$area
ggplot(cells.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), alpha = 0, color = 1) + theme_minimal()

# Plot redwood point data.
redwood <- ppinit("redwood.dat")
redwood.df <- data.frame(x = redwood$x, y = redwood$y)
redwood$area
ggplot(redwood.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -1, ymax = 0), alpha = 0, color = 1) + theme_minimal()

# Plot pines point data.
pines <- ppinit("pines.dat")
pines.df <- data.frame(x = pines$x, y = pines$y)
ggplot(pines.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = pines$area[1], xmax = pines$area[2], ymin = pines$area[3], ymax = pines$area[4]),
            alpha = 0, color = 1) + theme_minimal()
ggsave("figures/pines.pdf")

## b) ----
# L for cells
L.cells <- Kfn(ppinit("cells.dat"), fs = 0.5, k = 100)
L.cells.df <- data.frame(x = L.cells$x, y = L.cells$y)
ggplot(L.cells.df) + geom_line(aes(x, y)) + geom_abline(intercept = 0, slope = 1, color = "red") + 
  theme_minimal()

# L for redwood
L.redwood <- Kfn(ppinit("redwood.dat"), fs = 0.5, k = 100)
L.redwood.df <- data.frame(x = L.redwood$x, y = L.redwood$y)
ggplot(L.redwood.df) + geom_line(aes(x, y)) + geom_abline(intercept = 0, slope = 1, color = "red") + 
  theme_minimal()

# L for pines
L.pines <- Kfn(ppinit("pines.dat"), fs = 5, k = 100)
L.pines.df <- data.frame(x = L.pines$x, y = L.pines$y)
ggplot(L.pines.df) + geom_line(aes(x, y)) + geom_abline(intercept = 0, slope = 1, color = "red") + 
  theme_minimal()

## c) ----

sim.ppp <- function(xmin, xmax, ymin, ymax, n, fs){
  N <- 100 # Number of simulations
  k <- 100 # Number of grid points
  L.mat <- matrix(rep(NA, N*k), nrow = N)
  for(i in 1:N){
    x.sim <- runif(n, xmin, xmax)
    y.sim <- runif(n, ymin, ymax)
    L <- Kfn(list(x = x.sim, y = y.sim), fs = fs, k = k)
    L.mat[i, ] = L$y
  }
  return(list(x = L$x, y = L.mat))
}

# Plot CI for cells
n.cells <- length(cells$x)
cells.pois <- sim.ppp(0,1,0,1, n.cells, 0.5)

upper <- apply(cells.pois$y, 2,  quantile, probs = c(0.95))
lower <- apply(cells.pois$y, 2,  quantile, probs = c(0.05))
ggplot(data = data.frame(l = lower, u = upper, x = cells.pois$x)) + 
  geom_ribbon(aes(x = x, ymin = l, ymax = u), alpha = 0.2) +
  geom_line(data = L.cells.df, aes(x, y)) + theme_minimal()
  
# Plot CI for redwood
n.redwood <- length(redwood$x)
redwood.pois <- sim.ppp(0,1,-1,0, n.redwood, 0.5)

upper <- apply(redwood.pois$y, 2,  quantile, probs = c(0.95))
lower <- apply(redwood.pois$y, 2,  quantile, probs = c(0.05))
ggplot(data = data.frame(l = lower, u = upper, x = redwood.pois$x)) + 
  geom_ribbon(aes(x = x, ymin = l, ymax = u), alpha = 0.2) +
  geom_line(data = L.redwood.df, aes(x, y)) + theme_minimal()

# Plot CI for pines
n.pines <- length(pines$x)
pines.pois <- sim.ppp(0,9.6,0,10, n.pines, 5)

upper <- apply(pines.pois$y, 2,  quantile, probs = c(0.95))
lower <- apply(pines.pois$y, 2,  quantile, probs = c(0.05))
ggplot(data = data.frame(l = lower, u = upper, x = pines.pois$x)) + 
  geom_ribbon(aes(x = x, ymin = l, ymax = u), alpha = 0.2) +
  geom_line(data = L.pines.df, aes(x, y)) + theme_minimal()
