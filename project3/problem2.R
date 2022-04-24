setwd("~/Github/romstat/project3")
source("resources/functions.R")
load("resources/Admin1Geography.RData")

# Inverse of logit
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# Sample P_a|Y=y
sim.P <- function(n, Q, mu){
  P.mat <- matrix(NA, nrow = n, ncol = 37)
  x.var <- 1/diag(Q)
  for(i in 1:n){
    x <- rnorm(37, mu, sqrt(x.var))
    P.mat[i, ] <- expit(x)
  }
  return(P.mat)
}

v1 <- read.table("resources/Admin1Graph.txt")
d1 <- as.matrix(v1) %*% as.vector(rep(1, 37))
R1 <- diag(as.numeric(d1)) - unlist(v1)

direct.est.df <- read.table("resources/DirectEstimates.txt", header = TRUE, sep = "", dec = ".")
V <- direct.est.df$StdDev^2
y <- direct.est.df$Observation

cs.CV <- c(0, 0.3) # Color scale coefficient of variation
cs.median <- c(0, 1) # Color scale median

## b) ----

# sample P_A|Y
P.mat <- sim.P(100, diag(1/V), y)

# Find and plot median
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_b.pdf", 20, 20, P.median, nigeriaAdm1, "Median", cs.median)

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_b.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation", cs.CV)


## c) ----

# sample P_A|Y
inv.V.mat <- diag(1/V)
mu.cond <- solve(inv.V.mat + R1) %*% inv.V.mat %*% y
Q.cond <- inv.V.mat + R1
P.mat <- sim.P(100, Q.cond, mu.cond)

# Find and plot median
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_c.pdf", 20, 20, P.median, nigeriaAdm1, "Median", cs.median)

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_c.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation", cs.CV)

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

# Find and plot Median
P.mat <- sim.P(100, Q, mu)
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_d.pdf", 20, 20, P.median, nigeriaAdm1, "Median", cs.median)

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_d.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation", cs.CV)

## e) ----

## tau = 0.1
tau1 <- 0.1
mu1 <- solve(inv.V.mat + tau1*R1) %*% inv.V.mat %*% y
Q1 <- inv.V.mat + tau1*R1

# Find and plot Median
P.mat <- sim.P(100, Q1, mu1)
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_e1.pdf", 20, 20, P.median, nigeriaAdm1, "Median", cs.median)

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_e1.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation", cs.CV)

## tau = 10
tau2 <- 10
mu2 <- solve(inv.V.mat + tau2*R1) %*% inv.V.mat %*% y
Q2 <- inv.V.mat + tau2*R1

# Find and plot Median
P.mat <- sim.P(100, Q2, mu2)
P.median <- apply(P.mat, 2, median)
plotAreaCol("figures/P_median_e2.pdf", 20, 20, P.median, nigeriaAdm1, "Median", cs.median)

# Find and plot CV
P.sd <- apply(P.mat, 2, sd)
P.mean <- apply(P.mat, 2, mean)
P.CV <- P.sd/P.mean # Cofficient of variation
plotAreaCol("figures/P_CV_e2.pdf", 20, 20, P.CV, nigeriaAdm1, "Coefficient\nof Variation", cs.CV)


## f) ----

# rewrite to fit f) code
sim.P <- function(n, Q, mu){
  P.mat <- matrix(NA, nrow = n, ncol = 37)
  x.var <- 1/diag(Q)
  for(i in 1:n){
    x <- rnorm(37, mu, sqrt(x.var))
    # P.mat[i, ] <- expit(x)
    P.mat[i, ] <- x
  }
  return(P.mat)
}

sample.Pa = function(tau, R1, V.inv.mat, n=100){
  mu <- solve(tau*R1 + V.inv.mat) %*% V.inv.mat %*% y
  Q <- tau*R1 + V.inv.mat
  return(sim.P(n, Q, mu))
}

V <- direct.est.df$StdDev^2
y <- direct.est.df$Observation
V.mat <- diag(V)
inv.V.mat <- diag(1/V)
cond = sample.Pa(1, R1, inv.V.mat)

# relabel to fit newly introduced variable names
D = V
D.inv = inv.V.mat
R = R1
x.sim = apply(cond, 2, median)

loglikeY = function(tau){
  mu_c <- solve(tau*R + D.inv) %*% D.inv %*% y
  Qc <- tau*R + D.inv 
  
  
  return(36/2 * log(tau) 
         - tau/2 *t(x.sim) %*% R %*% x.sim
         - 1/2 * t(y - x.sim)%*%D.inv%*% (y-x.sim)
         - 1/2 * log(abs(det(Qc)))
         + 1/2 * t(x.sim - mu_c) %*% Qc %*% (x.sim - mu_c))
}

tau.test = seq(0.001, 5, length.out = 1e4)

tau.opt = optimize(loglikeY, tau.test, maximum = TRUE)
tau.opt

condOpt = sample.Pa(tau.opt$maximum, R, D.inv)

P.condOpt <- expit(condOpt)
P.medianOpt <- apply(P.condOpt, 2, median)
P.sdOpt <- apply(P.condOpt, 2, sd)
P.meanOpt <- apply(P.condOpt, 2, mean)
P.CVOpt <- P.sdOpt/P.meanOpt

plotAreaCol("figures/2f_P_medianOpt.pdf", 20, 20, P.medianOpt, nigeriaAdm1, "Median", cs.p)
plotAreaCol("figures/2f_P_CVOpt.pdf", 20, 20, P.CVOpt, nigeriaAdm1, "Coefficient\nof Variation")
