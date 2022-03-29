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
  sim.plot = matrix(NA, nrow=0, ncol = 2)
  len = c()
  for (i in 1:N){
    x = sim.NS(lambda_m_ext, lambda_c, sigma_c, d, w)
    L = Kfn(list(x = x[,1], y = x[,2]), fs = fs, k = k)
    # Store to matrix?
    L.mat[i, ] = L$y # all L-fn values
    if (i < 4){
      # sim.x.plot[i] = ( x[,1])
      # sim.y.plot[i] = ( x[,2])
      sim.plot = rbind(sim.plot,x)
      len = cbind(len, length(x[,2]))
    }
  }
  return(list(x = L$x, y = L.mat, sim.plot = sim.plot, len = len))
}


plot.pi = function(sim, L){
  upper <- apply(sim$y, 2,  quantile, probs = c(0.95))
  lower <- apply(sim$y, 2,  quantile, probs = c(0.05))
  
  gg.NS = ggplot(data = data.frame(l = lower, u = upper, x = sim$x)) + 
    geom_ribbon(aes(x = x, ymin = l, ymax = u), alpha = 0.2) +
    geom_line(data = L, aes(x, y), color="Cyan3") + 
    geom_abline(slope=1, intercept = 0, linetype=2) +
    xlim(0,0.5) + ylim(0,.6) +
    labs(color="Legend")
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
ns.p1 = plot.pi(redwood.ns1, L.redwood.df)

## ii) Improved fit ----
# second guess with kppm thomas process
thomasFit = kppm(as.ppp(redwood), clusters = "Thomas")
lambda_m = thomasFit$par[1]
lambda_c = thomasFit$mu
sigma_c = thomasFit$par[2]

c(lambda_m = lambda_m, lambda_c = lambda_c, sigma_c)

set.seed(321)
redwood.ns2 = sim.NS.seq(lambda_m,lambda_c, sigma_c)
ns.p2 = plot.pi(redwood.ns2, L.redwood.df)
ns2.u <- apply(redwood.ns2$y, 2,  quantile, probs = c(0.95))
ns2.l <- apply(redwood.ns2$y, 2,  quantile, probs = c(0.05))

# Third guess with extended window

set.seed(321)
redwood.ns3 = sim.NS.seq(lambda_m,lambda_c, sigma_c, d=0.5)
ns.p3 = plot.pi(redwood.ns3, L.redwood.df)
ns3.u <- apply(redwood.ns3$y, 2,  quantile, probs = c(0.95))
ns3.l <- apply(redwood.ns3$y, 2,  quantile, probs = c(0.05))

t=  (sum(ns3.u - ns3.l)/sum(ns2.u- ns2.l))
(1-t)*sum(ns2.u- ns2.l)
sum(ns2.u- ns2.l) - sum(ns3.u - ns3.l)
# show plots
ns.p1
ns.p2
ns.p3

"List final guestimates of model params."

## iii) Display ----
"Display data next to three realizations from guestimated NS model"
x = redwood.df$x
y = redwood.df$y
n = length(x)
n1 = (redwood.ns3$len[1])
n2 = (redwood.ns3$len[2])
n3 = (redwood.ns3$len[3])



head(redwood.ns3$sim.plot)

plot.df <- data.frame(x = c(x, redwood.ns3$sim.plot[,1]), 
                      y = c(y, redwood.ns3$sim.plot[,2]), 
                      idx = c(rep("Redwood dataset", n),
                              rep("Realization 1", n1),
                              rep("Realization 2", n2),
                              rep("Realization 3", n3)
                              ))
# for (i in 2:4){
#   plot.df$x[((i-1)*n+1):(i*n)] = x
#   # you need this shit
#   plot.df$y[((i-1)*n+1):(i*n)] = x
# }

ns.realizations = ggplot(plot.df) + 
  geom_point(aes(x = x, y = y)) + 
  facet_wrap(~idx, nrow = 2) + 
  theme_minimal()

length(redwood.ns3$sim.x.plot)


## saveFigs ----
ggsave("3_NSinit.pdf", plot = ns.p1, path = figpath, width=4, height = 4)
ggsave("3_NSkppm.pdf", plot = ns.p2, path = figpath, width=4, height = 4)
ggsave("3_NSext.pdf", plot = ns.p3, path = figpath, width=4, height = 4)
ggsave("3_NSrealizations.pdf", plot = ns.realizations, path = figpath, 
       width = 7, height = 7)



