### Problem 4 ----
library(spatial)
cells <- ppinit("cells.dat")
x <- cells$x
y <- cells$y
n <- length(x) # Number of points

# Calculate avg. distance between points / 2
point.mat <- matrix(c(x,y), nrow = length(x))
dist.mat <- dist(point.mat)
dist.mat[1]
mean(dist.mat)/2

# Simulate with Strauss
r0 <- 0.1 # Initial guess: r0 = 0.06
beta <- 6 # Initial guess: beta = 5
S <- Strauss(n, exp(-beta), r0)
ggplot(data.frame(x = S$x, y = S$y)) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), alpha = 0, color = 1) + theme_minimal()

# L for dataset
L.cells <- Kfn(cells, fs = 0.5, k = 100)
L.cells.df <- data.frame(x = L.cells$x, y = L.cells$y)

## Compute L 100 times for Strauss

# For ggplot
plot.df <- data.frame(x = c(x, rep(NA, 3*n)), y = c(y, rep(NA, 3*n)), 
                      idx = c(rep("Cells dataset", n), rep("Realization 1", n), 
                            rep("Realization 2", n), rep("Realization 4", n)))
sim.Strauss <- function(beta, r0, plot.df){
  set.seed(4250)
  L.mat <- matrix(NA, nrow = 100, ncol = 100)
  for(i in 1:100){
    S <- Strauss(n, exp(-beta), r0)
    if(i < 4){
      plot.df$x[(n + 1):((i + 1)*n)] = S$x
      plot.df$y[(n + 1):((i + 1)*n)] = S$y
    }
    L.Strauss <- Kfn(S, fs = 0.5, k = 100)
    L.mat[i, ] <- L.Strauss$y
  }
  return(list(L.mat = L.mat, plot.df = plot.df))
}

sim <- sim.Strauss(beta, r0, plot.df)
L.mat <- sim$L.mat
plot.df <- sim$plot.df

# Create prediction intervals
pred.int <- apply(L.mat, 2, quantile, probs = c(0.05, 0.95))

# Plot predints and  L function
ggplot(data.frame(x = L.cells$x, y = L.cells$y, lower = pred.int[1,], upper = pred.int[2,])) + 
  geom_ribbon(aes( x = x, ymin = lower, ymax = upper, fill = "Prediction interval"), alpha = 0.2) + 
  geom_line(aes(x, y, color = "Empirical L-function")) + 
  scale_color_manual(name = " ", values = c("Empirical L-function" = "slategray")) +
  scale_fill_manual(name = " ", values = c("Prediction interval" = "navyblue"))+
  ylab(" ") + ggtitle(expression(beta~" = 5, "~r[0]~" = 0.06")) + theme_minimal()

# Plot 3 realizations and cells dataset
ggplot(plot.df) + 
  geom_point(aes(x = x, y = y)) + 
  facet_wrap(~idx, nrow = 2) + 
  theme_minimal()
