setwd("~/Github/romstat/project3")
source("resources/functions.R")
load("resources/Admin1Geography.RData")

# Inverse of logit
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

v1 <- read.table("resources/Admin1Graph.txt")
d1 <- as.matrix(v1) %*% as.vector(rep(1, 37))
R1 <- diag(as.numeric(d1)) - unlist(v1)

direct.est.df <- read.table("resources/DirectEstimates.txt", header = TRUE, sep = "", dec = ".")
V <- direct.est.df$StdDev^2
y <- direct.est.df$Observation

## c) ----

# Sample from multivariate normal with precision matrix
sample.mvn <- function(n, mu, prec){
  X <- matrix(NA, nrow = n, ncol = length(mu))
  for(i in 1:n){
    L <- chol(prec)
    z <- rnorm(nrow(L))
    v <- solve(L, z)
    X[i, ] = v + mu
  }
  return(X)
}

# sample P_A|Y
V.mat <- diag(V)
inv.V.mat <- diag(1/V)
mu.cond <- -solve(inv.V.mat + R1) %*% R1 %*% y
Q.cond <- inv.V.mat + R1
cond <- sample.mvn(100, mu.cond, Q.cond)

# Find and plot median
P.cond <- expit(cond)
P.median <- apply(P.cond, 2, median)
plotAreaCol("figures/P_median.pdf", 20, 20, P.median, nigeriaAdm1, "Median")

# Find and plot CV
P.sd <- apply(P.cond, 2, sd)
P.mean <- apply(P.cond, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient of Variation")

## d) ----
y38 <- 0.5
y.tilde <- c(y, y38)
V.tilde <- c(V, 0.1^2)
Q1 <- R1
Q2 <- diag(1/V.tilde)
M <- matrix(0, nrow = 38, ncol = 37)
M[col(M) == row(M)] = 1
M[38, 19] = 1 # Kaduna is area number 19

mu <- solve(Q1 + t(M)%*%Q2%*%M) %*% t(M) %*% Q2 %*% y.tilde
Q <- Q1 + t(M)%*%Q2%*%M

sim.P <- function(n, Q, mu){
  P.mat <- matrix(NA, nrow = n, ncol = 37)
  x.var <- 1/diag(Q)
  for(i in 1:n){
    x <- rnorm(37, mu, sqrt(x.var))
    P.mat[i, ] <- expit(x)
  }
  return(P.mat)
}

# Find and plot Median
P.mat <- sim.P(100, Q, mu)
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_d.pdf", 20, 20, P.median, nigeriaAdm1, "Median")

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_d.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation")
