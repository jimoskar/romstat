### Problem 3
library(MASS)
library(fields)
library(geoR)
library(gridExtra)
library(latex2exp)
library(tidyverse)
library(reshape2)
set.seed(3)

## a)
sigma2 <- 2
a <- 3

D.tilde <- expand.grid(1:30,1:30)
field1 <- grf(1, grid = D.tilde, cov.model = "exponential", cov.pars = c(sigma2, a))

df <- data.frame(field = field1$data, x = D.tilde[,1], y = D.tilde[, 2])
ggplot(data = df, aes(x, y)) +
  geom_tile(aes(fill = field)) +
  scale_fill_viridis_c(name="Realization") +
  coord_fixed() + xlab("x") + ylab("y") + theme_minimal()

## b) 
semi.variogram <- function(x, ...){
  return(cov.spatial(0, ...) - cov.spatial(x, ...))
}

vario <- variog(coord = D.tilde, data = field1$data)
plot(vario, type = "l", lwd = 2, xlab = "h", ylab = expression(gamma(h)))
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("Empirical semi-variogram", "Analytic semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))

## c)
fields <- grf(nsim=3, grid = D.tilde, cov.model = "exponential", cov.pars = c(sigma2, a))
df <- data.frame(fields = c(fields$data[,1], fields$data[,2], fields$data[,3]),
                            x = rep(D.tilde[,1], 3), y = rep(D.tilde[, 2],3),
                 r = c(rep(1, 900), rep(2, 900), rep(3, 900)))

ggplot(data = df, aes(x, y)) +
  geom_tile(aes(fill = fields)) +
  scale_fill_viridis_c(name="Realization") + facet_wrap(~r) +
  coord_fixed() + xlab("x") + ylab("y")  + theme_minimal()

df1.vario <- list(data = fields$data[, 1], coords = D.tilde)
vario1 <- variog(df1.vario)
df2.vario <- list(data = fields$data[, 2], coords = D.tilde)
vario2 <- variog(df2.vario)
df3.vario <- list(data = fields$data[, 3], coords = D.tilde)
vario3 <- variog(df3.vario)


par(mfrow = c(1, 3))
plot(vario1, type = "l", lwd = 2, xlab = "h", ylab = expression(gamma(h)), ylim = c(0, 4))
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(vario2, type = "l", lwd = 2, xlab = "h", ylab = " ", ylim = c(0, 2.5))
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(vario3, type = "l", lwd = 2, xlab = "h", ylab = " ", ylim = c(0, 3))
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("Empirical semi-variogram", "Analytic semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))
dev.off()

## d)
set.seed(3)
idx <- sample(1:900, 36)

# Use the realizarion from a).
df.vario.d <- list(data = field1$data[idx], coords = D.tilde[idx, ])
vario.d <- variog(df.vario.d)

plot(vario.d, type = "l", lwd = 2, xlab = "h", ylab = expression(gamma(h)), ylim = c(0, 3.5))
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("Empirical semi-variogram", "True semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))

# All points on grid
ml.full <- likfit(data = field1$data, coords = D.tilde, cov.model = "exponential", ini.cov.pars = c(1,1))
# 36 points on grid
ml.reduced <- likfit(data = field1$data[idx], coords = D.tilde[idx,], 
                     cov.model = "exponential", ini.cov.pars = c(1,1))
par(mfrow = c(1, 2))
plot(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = ml.full$cov.pars), type = "l", 
     xlab = "h", ylab = expression(gamma(h)), ylim = c(0, 2.5),  main = "All grid points ")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = ml.reduced$cov.pars), type = "l", 
     xlab = "h", ylab = expression(rho(h)), ylim = c(0, 2.5), main = "36 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("ML-based semi-variogram", "True semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))

ml.reduced$cov.pars
ml.full$cov.pars
 expression("All grid points, "~hat(sigma) == ml.full$cov.pars[1])

## e)
set.seed(3)
idx9 <- sample(1:900, 9)
idx64 <- sample(1:900, 64)
idx100 <- sample(1:900, 100)

# Empirical
vario9 <- variog(coord = D.tilde[idx9,], data = field1$data[idx9])
vario64 <- variog(coord = D.tilde[idx64,], data = field1$data[idx64])
vario100 <- variog(coord = D.tilde[idx100,], data = field1$data[idx100])

par(mfrow = c(1, 3))
plot(vario9, type = "l", lwd = 2, xlab = "h", ylab = expression(gamma(h)), ylim = c(0, 4), main = "9 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(vario64, type = "l", lwd = 2, xlab = "h", ylab = " ", ylim = c(0, 4), main = "64 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(vario100, type = "l", lwd = 2, xlab = "h", ylab = " ", ylim = c(0, 4), main = "100 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("Empirical semi-variogram", "True semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))


# ML-based
ml9 <- likfit(data = field1$data[idx9], coords = D.tilde[idx9,], 
              cov.model = "exponential", ini.cov.pars = c(1,1))
ml64 <- likfit(data = field1$data[idx64], coords = D.tilde[idx64,], 
               cov.model = "exponential", ini.cov.pars = c(1,1))
ml100 <- likfit(data = field1$data[idx100], coords = D.tilde[idx100,], 
                cov.model = "exponential", ini.cov.pars = c(1,1))

par(mfrow = c(1, 3))
plot(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = ml9$cov.pars), type = "l", 
     xlab = "h", ylab = expression(gamma(h)), ylim = c(0, 2.5), main = "9 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = ml64$cov.pars), type = "l", 
     xlab = "h", ylab = " ", ylim = c(0, 2.5), main = "64 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")

plot(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = ml100$cov.pars), type = "l", 
     xlab = "h", ylab = " ", ylim = c(0, 3.5), main = "100 grid points")
lines(0:40, semi.variogram(0:40, cov.mod = "exponential", cov.pars = c(sigma2, a)), lty = "dashed")
legend("bottomright", c("ML-based semi-variogram", "True semi-variogram"), lwd = c(2,1), lty = c("solid", "dashed"))
ml9$cov.pars
