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
redwood$area
redwood.plot = ggplot(redwood.df) + 
  geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -1, ymax = 0),
            alpha = 0, color = 1) + 
  theme_minimal()
# redwood.plot


## i) Empirical fit ----
"Make empirical fit"
N = data.frame(n=c(5, 9, 8, 13, 19, 6))
a = mean(N$n)
lambda = 6

sim.ns <- function(w, lambda, d=0.5, lambda_c=10, sigma_c=c(1,1), 
                   N=100, k=100, fs=0.7){
  # w: Window
  # lambda: Parent intensity
  # d: Distance to extend window
  # pc: Pmf #daughters (removed)
  # fr: pdf parent-daughter distance (removed)
  # cMean: mean number of daughter points per parent
  # N: #sims
  # k: #grid points

  xmin = w[1,1]; xmax = w[1,2]
  ymin = w[2,1]; ymax = w[2,2]
  wExt = w + matrix(c(-d, -d, d, d), ncol = 2) # Extended window
  nu   = (wExt[1,2] - wExt[1,1])*(wExt[2,2] - wExt[2,1]) # nu(W)
  n    = rpois(N, lambda*nu) # Create #parents for each sim
  L.mat = matrix(rep(NA, N*k), nrow = N) # Allocate memory for L estimates
  
  for(i in 1:N){
    # parent points
    y1 = runif(n[i], w[1,1], w[1,2]) # Not extended window
    y2 = runif(n[i], w[2,1], w[2,2])
    # y1 = runif(n[i], wExt[1,1], wExt[1,2])
    # y2 = runif(n[i], wExt[2,1], wExt[2,2])
    # Daughter points
    x1 = c() # x daughter val
    x2 = c() # y daughter val
    nt = n[i]
    for (p in 1:nt){
      # browser()
      c = rpois(nt, lambda=lambda_c) # #daughter points for this parent (p_c)
      x1 = c(x1, rnorm(c, mean = y1[p], sd = sigma_c[1]^2)) # f_R & \hat r
      x2 = c(x2, rnorm(c, mean = y2[p], sd = sigma_c[2]^2)) # f_R & \hat r
    }
    # browser()
    # plot(x1[1:(10*n[1])], x2[1:(10*n[1])], main="NS sim",
    #      xlim=c(xmin-d,xmax+d), ylim=c(ymin-d, xmin-d))
    # tp = 3
    # plot(x1[(10*sum(n[1:tp])):(10*sum(n[1:(tp+1)]))], 
    #      x2[(10*sum(n[1:tp])):(10*sum(n[1:(tp+1)]))], 
    #      main="NS sim", xlim=c(xmin-d,xmax+d), ylim=c(ymin-d, xmin-d))
    
    # Remove daughter points outside original window {x: x \in w}
    where = x1>=rep(xmin,length(x1)) && x1<=rep(xmax,length(x1)) &&
      x2>=rep(ymin,length(x2)) && x2<=rep(ymax,length(x2))
    sum(where)
    x1 = x1[where] 
    x2 = x2[where] 
    L = Kfn(list(x = x1, y = x2), fs = fs, k = k)
    # Store to data frame (matrix?)
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




"Evaluate by simulating 100 realizations from NS to produce 90% pi. of L-fnc"
# Plot PI for redwood sim by NS process
# redwood.pois <- sim.ppp(0,9.6,0,10, n.pines, 5) # ur l.mat
# n.pines <- length(pines$x)
w = matrix(c(0, 1,
             -1, 0), ncol = 2, byrow = T)
L.redwood <- Kfn(ppinit("redwood.dat"), fs = 0.5, k = 100)
L.redwood.df <- data.frame(x = L.redwood$x, y = L.redwood$y)

set.seed(321)
redwood.ns <- sim.ns(w, 10, sigma_c = c(0.01,0.01))

ns.p1 = plot.pi(redwood.ns, L.redwood.df)
ns.p1

set.seed(321)
sc = 0.17
redwood.ns <- sim.ns(w, 23, sigma_c = c(sc,sc), lambda_c = 2.63)

ns.p2 = plot.pi(redwood.ns, L.redwood.df)
ns.p2

temp = kppm(redwood.df, "Thomas")






## ii) Improved fit ----
"Iterate guestimate to improve fit"
"Try to make empirical L-fnc consistent wit 90% pi."
"List final guestimates of model params."

## iii) Display ----
"Display data next to three realizations from guestimated NS model"


## saveFigs ----
ggsave("3_NS1.pdf", plot = ns.p1, path = figpath, width=4, height = 4)



sim.NS = function(lm, lc, sc){
  result = matrix(NA, nrow=0, ncol=2)
  k=0
  km = rpois(1,lm)
  for (j in 1:km){
    xm = c(runif(1,0,1), runif(1,-1,0))
    kc = rpois(1,lc)
    for (i in 1:kc) {
      xc = mvrnorm(mu = xm, Sigma = sc*diag(2))
      if (xc[1] > 0 && xc[1] < 1 && xc[2] > -1 && xc[2] < 0){
        result = rbind(result, xc)
        k = k + kc
      }
    }
  }
  return(result)
}