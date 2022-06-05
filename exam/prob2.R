library(geoR)
library(MASS)
library(ggplot2)
library(wesanderson)
library(dplyr)

curve(1/5 * cov.spatial(x, cov.mod = "exponential", cov.pars = c(5, 20)), 
      from = 0, to = 100, col = 2, xlab = "distance", ylab = expression(rho(h)))
curve(1/5 *cov.spatial(x, cov.mod = "exponential", cov.pars = c(5, 10)), 
      from = 0, to = 100, col = 2, lwd = 2, add = TRUE)

n <- 50 # Number of grid points
D.tilde <- 1:n # Grid
sigma2 = 10
a = 10

mu <- rep(0, n)
Sigma <- cov.spatial(as.matrix(dist(expand.grid(D.tilde))), 
                       cov.mod = "exponential", 
                       cov.pars = c(sigma2, a))
X <- mvrnorm(1, mu, Sigma)
obs <- c(10, 25, 30)
X.obs <- X[obs]
Sig <-  cov.spatial(as.matrix(dist(expand.grid(obs))), 
                            cov.mod = "powered.exponential", 
                            cov.pars = c(sigma2, a), kappa = 1.9)
c <-  cov.spatial(as.matrix(dist(expand.grid(D.tilde)))[,obs], 
                  cov.mod = "powered.exponential", 
                  cov.pars = c(sigma2, a), kappa = 1.9)
X.krig <- c %*% solve(Sig) %*% X.obs

plot(X.krig, type = "l")

