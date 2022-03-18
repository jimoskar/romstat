# Project 2, 18.03.2022
# Problem 3 ----
library(spatial)
library(MASS)
library(tidyverse)
library(spatstat)

# Config
setwd('~/Fysmat/8 Semester V2022/RomStat/romstat/project2/')
figpath = './figures'
# getwd()

## Read,plot ----
redwood    <- ppinit("redwood.dat")
redwood.df <- data.frame(x = redwood$x, y = redwood$y)

# Plot

t1 = c(xl=1,xu=2)
redwood.plot = ggplot(redwood.df) + 
  geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -1, ymax = 0),
            alpha = 0, color = 1) + 
  theme_minimal()
# redwood.plot


## i) Empirical fit ----
"Make empirical fit"
sim.NS = function(lm, lc, sc, d = 0, w=c(xl=0, xu=1, yl=-1, yu=0)){
  result = matrix(NA, nrow=0, ncol=2)
  lme = lm*((w[2]-w[1]+2*d)*(w[4]-w[3]+2*d)) # intensity on extended window
  k=0
  km = rpois(1,lm)
  for (j in 1:km){
    xm = c(runif(1,w[1]-d,w[2]+d), runif(1,w[3]-d,w[4]+d)) # parent location
    # xm = c(runif(1,0-d,1+d), runif(1,-1-d,0+d)) # parent location
    kc = rpois(1,lc) # #children for this parent
    for (i in 1:kc) {
      xc = mvrnorm(mu = xm, Sigma = sc*diag(2))
      if (xc[1] > w[1] && xc[1] < w[2] && xc[2] > w[3] && xc[2] < w[4]){
      # if (xc[1] > 0 && xc[1] < 1 && xc[2] > -1 && xc[2] < 0){
        result = rbind(result, xc)
        k = k + kc # Total #points
      }
    }
  }
  return(result)
}


sim.NS.seq = function(lambda_m, lambda_c, sigma_c, N = 100, k=100, 
                      fs=0.7, d=0, w=c(xl=0, xu=1, yl=-1, yu=0)){
  L.mat = matrix(rep(NA, N*k), nrow = N)
  lambda_m_ext = lambda_m*((w[2]-w[1]+2*d)*(w[4]-w[3]+2*d))
  for (i in 1:N){
    x = sim.NS(lambda_m_ext, lambda_c, sigma_c, d, w)
    L = Kfn(list(x = x[,1], y = x[,2]), fs = fs, k = k)
    # Store to matrix?
    L.mat[i, ] = L$y # all L-fn values
  }
  return(list(x = L$x, y = L.mat))
}


plot.pi = function(sim, L){
  upper <- apply(sim$y, 2,  quantile, probs = c(0.95))
  lower <- apply(sim$y, 2,  quantile, probs = c(0.05))
  
  gg.NS = ggplot(data = data.frame(l = lower, u = upper, x = sim$x)) + 
    geom_ribbon(aes(x = x, ymin = l, ymax = u), alpha = 0.2) +
    geom_line(data = L, aes(x, y), color="cyan3") + 
    geom_abline(slope=1, intercept = 0, linetype=2) +
    xlim(0,0.5) + ylim(0,.6) +
    theme_minimal()
  gg.NS
}






# lambda_m = 20, lambda_c = 3, sigma_c = 0.002
L.redwood <- Kfn(ppinit("redwood.dat"), fs = 0.7, k = 100)
L.redwood.df <- data.frame(x = L.redwood$x, y = L.redwood$y)

# Initial guess (6, 10, 0.01)
Nc = data.frame(n=c(5, 9, 8, 13, 19, 6)) # Points per cluster
lambda_m = 6
lambda_c = mean(Nc$n) # Daughter intensity
sigma_c = 0.01 # simply guessed

set.seed(321)
redwood.ns1 = sim.NS.seq(lambda_m,lambda_c, sigma_c)
plot.pi(redwood.ns1, L.redwood.df)

# second guess with fnc(6, 10, 0.001)
set.seed(321)
redwood.ns2 = sim.NS.seq(lambda_m=6,
                         lambda_c=10, 
                         sigma_c =0.001)
plot.pi(redwood.ns2, L.redwood.df)

temp = kppm(as.ppp(redwood), "Thomas")






## ii) Improved fit ----
"Iterate guestimate to improve fit"
"Try to make empirical L-fnc consistent wit 90% pi."
"List final guestimates of model params."

## iii) Display ----
"Display data next to three realizations from guestimated NS model"


## saveFigs ----
ggsave("3_NS1.pdf", plot = ns.p1, path = figpath, width=4, height = 4)


