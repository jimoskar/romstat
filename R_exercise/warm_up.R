# 1.1.1
m1 <- matrix(nrow = 6, ncol = 10)

for(i in 1:ncol(m)){
  m1[,i] = rnorm(nrow(m))
}


# Starting out with NA's is useful if you do not know the values beforehand.

# 1.1.2
x <- 1:10
y <- 11:20
z <- runif(10)
m2 <- matrix(data = c(x,y,z), nrow = 3, ncol = 10, byrow = TRUE)


# 1.1.3
df1 <- data.frame(m1)
df2 <- data.frame(m2)


# 1.1.4
d <- diag(x = 1:10)

# 1.1.5
d %*% x

# 1.2.1
library(ggplot2)
ggplot(df1) + geom_point(aes(x = X1, y = X2, color = 1:6)) + 
  ggtitle("hei")

# 1.2.2
library(fields)
image.plot(x = 1:6, y = 1:10, z = m1)

# 1.2.3
image.plot(x = m2[1, ], y = m2[2,], z = diag(nrow = ncol(m2)) * m2[3, ])

# 

