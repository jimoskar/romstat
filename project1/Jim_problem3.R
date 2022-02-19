### Problem 2 (Jim)
library(akima)
library(fields)
library(ggplot2)
library(viridis)

topo <- read.csv("~/Github/romstat/project1/topo.dat", sep="")
head(topo)

## a)
topo.interp <- with(topo, interp(x, y, z, nx = 100, ny = 100))

# Set-up for color in `persp()` plot
nbcol <- 100
color <- hcl.colors(nbcol, palette = "Terrain 2")
zfacet <- with(topo.interp, z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])
facetcol <- cut(zfacet, nbcol)

# 2D plot
image.plot(topo.interp, col = color)
contour(topo.interp, nlevels = 11, add = TRUE)
with(topo, points(x, y, pch = 4, col = "white"))
with(topo, points(x, y, pch = 3, col = "white"))


persp(topo.interp, scale = FALSE,
      theta = 145, phi = 20,
      col = color[facetcol],
      lwd = 0.1, border = rgb(0, 0, 0, alpha = 0.8), box = TRUE,
      xlab = "X", zlab = "Z"
)

