## Example 1
simBrownian = function(xMin, xMax, nRes, sigma = 1){
  h = (xMax-xMin)/(nRes-1)
  return(list(x = seq(xMin, xMax, length.out = nRes),
              B = cumsum(c(0, rnorm(nRes-1, sd = sqrt(h)*sigma)))))
}

# a) Brownian motion
set.seed(100)
B1 = simBrownian(0, 10, 200)
B2 = simBrownian(0, 5, 200)
plot(B1$x, B1$B, type = "l", ylim = c(-2,2), xlim = c(0,5), main = "Brownian motion", xlab = "s", ylab = "X")
lines(B2$x, B2$B, col = "red")

# b) Coordinate-wise Brownian motion
library(ggplot2)
library(viridis)
df = data.frame(s1 = rep(B1$x, length(B2$x)), 
                s2 = rep(B2$x, each = length(B1$x)),
                X = rep(B1$B, length(B2$B))+rep(B2$B, each = length(B1$B)))

fig = ggplot(data = df, aes(s1, s2)) + 
  geom_raster(aes(fill = X)) +
  scale_fill_viridis_c() + coord_fixed() + 
  ggtitle("Coordinate-wise Brownian motion")
fig

# c) Brownian sheet
simSheet = function(xRan, yRan, dim, sigma = 1){
  dA = diff(xRan)*diff(yRan)/prod(dim-1)
  x = seq(xRan[1], xRan[2], length.out = dim[1])
  y = seq(yRan[1], yRan[2], length.out = dim[2])
  z = matrix(0, nrow = dim[1], ncol = dim[2])
  for(i in 1:dim[1]){
    for(j in 1:dim[2]){
      if(!(i == 1 && j == 1)){
        z[i,j] = rnorm(1, sd = sigma*sqrt(dA))
      }
      if(i > 1){
        z[i,j] = z[i,j] + z[i-1,j]
      }
      if(j > 1){
        z[i,j] = z[i,j] + z[i,j-1]
      }
      if(i > 1 && j > 1){
        z[i,j] = z[i,j] - z[i-1, j-1]
      }
    }
  }
  return(data.frame(s1 = rep(x, length(y)),
                    s2 = rep(y, each = length(x)),
                    X  = as.vector(z)))
}

df = simSheet(c(0, 10), c(0, 5), dim = c(200, 200))  
fig = ggplot(data = df, aes(s1, s2)) + 
  geom_raster(aes(fill = X)) +
  scale_fill_viridis_c() + 
  coord_fixed() +
  ggtitle("Brownian sheet")
fig
