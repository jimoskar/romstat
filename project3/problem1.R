source("resources/functions.R")

### Problem 1 ----
## a) ----
v1 <- read.table("resources/Admin1Graph.txt")
v2 <- read.table("resources/Admin2Graph.txt")

make.R <- function(dim, ne){
  R <- matrix(NA, nrow = dim, ncol = dim)
  for(i in 1:dim){
    for(j in 1:dim){
      if(i == j){
        R[i,j] <- sum(ne[i, 2:dim])
      } else{
        R[i,j] <- ifelse(ne[i,j], -1, 0)
      }
    }
  }
  return(R)
}

# Construct R1 and R2
R1 <- make.R(37, v1)
R2 <- make.R(775, v2)

# Compute proportion of non-zero elements
sum(R1 > 0)/(37*37)
sum(R2 > 0)/(775*775)

# Create plots of sparsity patterns
plot.spars <- function(mat){
  spars <- (mat != 0)
  dim(spars) <- NULL
  spars.df <- data.frame(fill = as.numeric(spars), coord = expand.grid(1:nrow(mat), 1:nrow(mat)))
  
  p <- ggplot(spars.df, aes(coord.Var1, coord.Var2)) +
    geom_tile(aes(fill = fill)) + 
    coord_fixed() + xlab("Column") + ylab("Row")  + theme_minimal() + theme(legend.position="none")
  return(p)
}

p1 <- plot.spars(R1)
ggsave("figures/spars1.pdf", plot = p1)

p2 <- plot.spars(R2)
ggsave("figures/spars2.png", plot = p2)

## b) ----
load("resources/Admin1Geography.RData")
load("resources/Admin2Geography.RData")

sim.1o.intrinsic <- function(Q, eps){
  n <- nrow(Q)
  Q.hat <- Q + eps*diag(n)
  L.hat <- chol(Q.hat)
  Z <- rnorm(n)
  v <- solve(L.hat, Z)
  x <- v - mean(v)
  return(x)
}

set.seed(4250)
cs = c(-3, 3) # color scale
leg = " "
x.besag1 <- sim.1o.intrinsic(1*R1, 2)
x.besag2 <- sim.1o.intrinsic(1*R1, 2)
plotAreaCol("figures/besag1.pdf", 20, 20, x.besag1, nigeriaAdm1, leg, cs)
plotAreaCol("figures/besag2.pdf", 20, 20, x.besag2, nigeriaAdm1, leg, cs)

x.norm1 <- rnorm(37)
x.norm2 <- rnorm(37)
plotAreaCol("figures/norm1.pdf", 20, 20, x.norm1, nigeriaAdm1, leg, cs)
plotAreaCol("figures/norm2.pdf", 20, 20, x.norm2, nigeriaAdm1, leg, cs)


  
  