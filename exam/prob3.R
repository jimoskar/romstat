
a <- 1
b <- 200
lam <- function(s1){
  return((a + b*s1)/10)
}

ippp <- function(W = c(0,10,0,10)){
  X <- matrix(NA, nrow = 1, ncol = 2)
  max.lam <- lam(W[2])
  v <- 100 # Volume
  n.star <- rpois(1, max.lam*v)
  n <- 0
  for(i in 1:n.star){
    x1.star <- runif(1,W[1],W[2])
    x2.star <- runif(1,W[3],W[4])
    u <- runif(1)
    if(u <= lam(x1.star)/max.lam){
      X <- rbind(X, c(x1.star, x2.star))
      n <- n + 1
    }
  }
  return(list(X = X[2:nrow(X), ], n = n))
  
}

x <- ippp()
plot(x$X)
