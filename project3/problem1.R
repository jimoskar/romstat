source("resources/functions.R")

### Problem 1 ----
## a) ----
v1 <- read.table("resources/Admin1Graph.txt")
v2 <- read.table("resources/Admin2Graph.txt")

# Construct R1 and R2
d1 <- as.matrix(v1) %*% as.vector(rep(1, 37))
R1 <- diag(as.numeric(d1)) - unlist(v1)

d2 <- as.matrix(v2) %*% as.vector(rep(1, 775))
R2 <- diag(as.numeric(d2)) - unlist(v2)

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
    coord_fixed() + scale_y_reverse() + 
    xlab("Column") + ylab("Row")  + theme_minimal() + theme(legend.position="none")
  return(p)
}

p1 <- plot.spars(R1)
ggsave("figures/spars1.pdf", plot = p1, height = 5, width = 5)

p2 <- plot.spars(R2)
ggsave("figures/spars2.png", plot = p2, height = 5, width = 5)

## b) ----
load("resources/Admin1Geography.RData")
load("resources/Admin2Geography.RData")

sim.1o.intrinsic <- function(Q, eps){
  n <- nrow(Q)
  Q.hat <- Q + eps*diag(n)
  L.hat <- t(chol(Q.hat))
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

## c) ----

set.seed(4250)
cs = c(-3, 3) # color scale
leg = " "
x.besag1 <- sim.1o.intrinsic(1*R2, 2)
x.besag2 <- sim.1o.intrinsic(1*R2, 2)
plotAreaCol("figures/besag_admin2_1.pdf", 20, 20, x.besag1, nigeriaAdm2, leg, cs)
plotAreaCol("figures/besag_admin2_2.pdf", 20, 20, x.besag2, nigeriaAdm2, leg, cs)

x.norm1 <- rnorm(775)
x.norm2 <- rnorm(775)
plotAreaCol("figures/norm_admin2_1.pdf", 20, 20, x.norm1, nigeriaAdm2, leg, cs)
plotAreaCol("figures/norm_admin2_2.pdf", 20, 20, x.norm2, nigeriaAdm2, leg, cs)


## d) ----

# Create  matrix and realizations
besag.mat <- matrix(NA, nrow = 100, ncol = 775)
for(i in 1:100){
  besag.mat[i, ] = sim.1o.intrinsic(1*R2, 2)
}

# Caluclate empirical variance and plot it
besag.var <- apply(besag.mat, 2, var)
plotAreaCol("figures/besag_var.pdf", 20, 20, besag.var, nigeriaAdm2, "Empirical variance", c(0,1))

# Calculate cor of Gubio (150) with all others
besag.cor <- cor(besag.mat[,150], besag.mat)
plotAreaCol("figures/besag_cor.pdf", 20, 20, t(besag.cor), nigeriaAdm2, "Empircal correlation", c(-0.3, 1))
besag.cor.neg <- as.numeric(besag.cor > 0)
plotAreaCol("figures/besag_cor_neg.pdf", 20, 20, besag.cor.neg, nigeriaAdm2, "Empircal correlation", c(0, 1))

  