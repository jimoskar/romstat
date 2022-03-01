
### Problem 1

library(geoR)
library(MASS)
library(ggplot2)
library(wesanderson)
library(dplyr)

## a)

# Correlation functions
curve(1/5 * cov.spatial(x, cov.mod = "powered.exponential", cov.pars = c(5, 10), kappa = 1), 
      from = 0, to = 100, col = 2, xlab = "distance", ylab = expression(rho(h)))
curve(1/5 *cov.spatial(x, cov.mod = "powered.exponential", cov.pars = c(5, 10), kappa = 1.9), 
      from = 0, to = 100, col = 2, lwd = 2, add = TRUE)
curve(1/5 * cov.spatial(x, cov.mod = "matern", cov.pars = c(5, 20), kappa = 1), 
      from = 0, to = 100, col = 3, add = TRUE)
curve(1/5 *cov.spatial(x, cov.mod = "matern", cov.pars = c(5, 20), kappa = 3), 
      from = 0, to = 100, col = 3, lwd = 2, add = TRUE)
legend("topright", c(expression("Powered exponential, " ~alpha~ " = 1"), 
                     expression("Powered exponential, " ~alpha~ " = 1.9"),        
                     expression("Matérn, " ~nu~ " = 1"), expression("Matérn, " ~nu~ " = 3")),
       col = c(2,2,3,3), lwd = c(1,2,1,2))

# Semi-variograms
sigma2 <- 5
semi.variogram <- function(x, ...){
  return(cov.spatial(0, ...) - cov.spatial(x, ...))
}

curve(semi.variogram(x, cov.mod = "powered.exponential", cov.pars = c(sigma2, 10), kappa = 1), 
      from = 0, to = 100, col = 4, xlab = "distance", ylab = expression(gamma(h)))
curve(semi.variogram(x, cov.mod = "powered.exponential", cov.pars = c(sigma2, 10), kappa = 1.9), 
      from = 0, to = 100, col = 4, lwd = 2, add = TRUE)
curve(semi.variogram(x, cov.mod = "matern", cov.pars = c(sigma2, 20), kappa = 1), 
      from = 0, to = 100, col = 5, add = TRUE)
curve(semi.variogram(x, cov.mod = "matern", cov.pars = c(sigma2, 20), kappa = 3), 
      from = 0, to = 100, col = 5, lwd = 2, add = TRUE)
legend("bottomright", c(expression("Powered exponential, " ~alpha~ " = 1"), 
                        expression("Powered exponential, " ~alpha~ " = 1.9"),
                        expression("Matérn, " ~nu~ " = 1"), 
                        expression("Matérn, " ~nu~ " = 3")),
       col = c(4,4,5,5), lwd = c(1,2,1,2))
title(main = expression(~sigma~ "= 5"))

## b)

# Exp. cov.mod.
n <- 50 # Number of grid points
D.tilde <- 1:n # Grid

df.exp <- data.frame(x = rep(D.tilde, 4), 
                     y1 = rep(NA, 4), y2 = rep(NA, 4), y3 = rep(NA, 4), y4 = rep(NA, 4),
                     combination = rep(NA, 4 * n))
for(i in 1:nrow(params.exp)){
  mu <- rep(0, n)
  Sigma <- cov.spatial(as.matrix(dist(expand.grid(D.tilde))), 
                       cov.mod = "powered.exponential", 
                       cov.pars = c(params.exp$sigma2[i], params.exp$a.exp[i]), 
                       kappa = params.exp$alpha[i])
  X <- mvrnorm(4, mu, Sigma)
  df.exp$y1[((i-1)*n + 1):(i*50)] = X[1, ]
  df.exp$y2[((i-1)*n + 1):(i*50)] = X[2, ]
  df.exp$y3[((i-1)*n + 1):(i*50)] = X[3, ]
  df.exp$y4[((i-1)*n + 1):(i*50)] = X[4, ]
  df.exp$combination[((i-1)*n + 1):(i*50)] = rep(i, 50)
}

palette <- wes_palette("GrandBudapest2", n = 4)
string1 = rep("list(sigma^2 == 1, alpha ==1, a == 10)", 50)
string2 = rep("list(sigma^2 == 5, alpha ==1, a == 10)", 50)
string3 = rep("list(sigma^2 == 1, alpha ==1.9, a == 10)", 50)
string4 = rep("list(sigma^2 == 5, alpha ==1.9, a == 10)", 50)

df.exp$combination[df.exp$combination == 1] = string1
df.exp$combination[df.exp$combination == 2] = string2
df.exp$combination[df.exp$combination == 3] = string3
df.exp$combination[df.exp$combination == 4] = string4

ggplot(df.exp, aes(x = x)) + geom_line(aes(y = y1), color = palette[1]) + 
  geom_line(aes(y = y2), color = palette[2]) +
  geom_line(aes(y = y3), color = palette[3]) + 
  geom_line(aes( y = y4), color = palette[4]) +
  facet_wrap( ~combination, nrow = 2, labeller = label_parsed) + xlab("Discretized grid") + ylab("Realization") + theme_minimal()

# Matern cov.mod
df.matern <- data.frame(x = rep(D.tilde, 4), 
                        y1 = rep(NA, 4), y2 = rep(NA, 4), y3 = rep(NA, 4), y4 = rep(NA, 4),
                        combination = rep(NA, 4 * n))

for(i in 1:nrow(params.matern)){
  mu <- rep(0, n)
  Sigma <- cov.spatial(as.matrix(dist(expand.grid(D.tilde))), 
                       cov.mod = "matern", 
                       cov.pars = c(params.matern$sigma2[i], params.matern$a.matern[i]), 
                       kappa = params.matern$nu[i])
  X <- mvrnorm(4, mu, Sigma)
  df.matern$y1[((i-1)*n + 1):(i*50)] = X[1, ]
  df.matern$y2[((i-1)*n + 1):(i*50)] = X[2, ]
  df.matern$y3[((i-1)*n + 1):(i*50)] = X[3, ]
  df.matern$y4[((i-1)*n + 1):(i*50)] = X[4, ]
  df.matern$combination[((i-1)*n + 1):(i*50)] = rep(i, 50)
}

palette <- wes_palette("GrandBudapest1", n = 4)
string1 = rep("list(sigma^2 == 1, nu ==1, a == 20)", 50)
string2 = rep("list(sigma^2 == 5, nu ==1, a == 20)", 50)
string3 = rep("list(sigma^2 == 1, nu ==3, a == 20)", 50)
string4 = rep("list(sigma^2 == 5, nu ==3, a == 20)", 50)

df.matern$combination[df.matern$combination == 1] = string1
df.matern$combination[df.matern$combination == 2] = string2
df.matern$combination[df.matern$combination == 3] = string3
df.matern$combination[df.matern$combination == 4] = string4

ggplot(df.matern, aes(x = x)) + geom_line(aes(y = y1), color = palette[1]) + 
  geom_line(aes(y = y2), color = palette[2]) +
  geom_line(aes(y = y3), color = palette[3]) + 
  geom_line(aes( y = y4), color = palette[4]) +
  facet_wrap( ~combination, nrow = 2, labeller = label_parsed) + xlab("Discretized grid") + ylab("Realization") + theme_minimal()

## d)

s = c(10, 25, 30) # Positions of interest
realization = filter(df.matern, 
                     combination == "list(sigma^2 == 1, nu ==1, a == 20)")[, 2] # First realization over entire grid
y = realization[s]

D.tilde <- 1:50
Sigma_XX <- cov.spatial(as.matrix(dist(expand.grid(D.tilde))), 
                        cov.mod = "matern", 
                        cov.pars = c(5, 20), 
                        kappa = 1)
Sigma_XX
Sigma_XY <- Sigma_XX[ , s]
# Construct Sigma_YY
H <- matrix(rep(0, 3*50), nrow = 3, ncol = 50)
H[1, 10] = 1
H[2, 25] = 1
H[3, 30] = 1

sigma2_N <- 0
Sigma_YY <- H %*% Sigma_XX %*% t(H) + sigma2_N*diag(3) # sigma_N^2 = 0.25
inv.YY <- solve(Sigma_YY)
BLUP <- Sigma_XY %*% inv.YY %*% y
var.BLUP <- Sigma_XX  - Sigma_XY %*% inv.YY %*% t(Sigma_XY)
point.variance <- diag(var.BLUP)
point.variance[point.variance < 0] = 0 # Remove small negative values due to numerical imprecision.
diag(var.BLUP) <- point.variance
z <- qnorm(0.95)
conf.int.upper <- BLUP + sqrt(point.variance)*z
conf.int.lower <- BLUP - sqrt(point.variance)*z
ggplot(data.frame(x = D.tilde, y = BLUP, lower = conf.int.lower, upper = conf.int.upper), aes(x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = "blue", fill = "blue", alpha = 0.1) + 
  geom_line(aes(y = y), color = "red") + xlab("Discretized grid") + ylab("Realization") + theme_minimal()

## e
XX <- matrix(rep(NA, 50 *100), nrow = 100) # matrix of realizations
for(i in 1:100){
  XX[i, ] <-  mvrnorm(1, BLUP, var.BLUP)
}

variance <- sapply(as.data.frame(XX), var)
mean <- colMeans(XX)
plot.df <- data.frame(x = D.tilde, y = mean, lower = mean - z*sqrt(variance), 
                      upper = mean + z*sqrt(variance))
p1 <- ggplot(data = as.data.frame(t(XX)), aes(x = D.tilde))
for (i in 1:100) {
  p1 <- p1 + geom_line(aes_string(y = paste0("V", i)), alpha = 0.2, color = "blue")
}
p1 +  geom_ribbon(data = plot.df, aes(ymin = lower, ymax = upper), alpha = 0.1, fill = "blue", color = "blue") +
  geom_line(data = plot.df, aes(y = y, color = "red"), color = "red") + 
  ylab("Empirical prediction") + xlab("Discretized grid") + theme_minimal()

## f)
A <- function(x){
  return(sum((x>2)*(x-2)))
}

# Use realizations
temp <- apply(XX, 1, A)
A.hat <- mean(temp)
A.hat
temp
var(temp)
# Use simple Kriging
A.tilde <- A(BLUP)
A.tilde
BLUP

