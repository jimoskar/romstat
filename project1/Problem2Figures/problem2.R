library(geoR)    
library(akima)   
library(fields)   
library(MASS)     
library(viridis)  
library(tidyverse)
library(spatstat)
library(scatterplot3d)
library(plotly)

# Surface 3D plot function
s3dPlot = function(grid, fn="surfacePlot.png",z=1.8, x=1, y=1){
  fig = plot_ly(x=grid$x, y=grid$y,z = grid$z) %>% 
    add_surface(
      colors = terrain.colors(5),
      contours = list(
        z = list(
          show=TRUE,
          usecolormap=F,
          highlightcolor="#ff1111",
          project=list(z=TRUE)
        ),
        start = 0,
        end=1000,
        size = 0
      )
    ) 
  m <- list(
    l = 1,
    r = 1,
    b = 1,
    t = 1,
    pad = 100
  )
  fig <- fig %>% layout(
    margin = m,
    scene = list(
      xaxis = list(title = "y"), 
      yaxis = list(title = "x"),
      zaxis = list(title = "Elevation"),
      camera=list(
        eye = list(x=x, y=y, z=z)
      )
    )
  )%>% 
    config(
      toImageButtonOptions = list(
        format = "png",
        filename = "myplot",
        width = 800,
        height = 700
      )
    )
  save_image(fig, file=fn)
  return(fig)
}


# Misc
figPath = "./Fysmat/8 Semester V2022/RomStat/romstat/project1/Problem2Figures/"
dt <- read.table(file="~/Fysmat/8 Semester V2022/RomStat/Project1/Prj1Code_romstat/topo.dat")
setwd("C:/Users/chris/OneDrive/Documents/Fysmat/8 Semester V2022/RomStat/romstat/project1/Problem2Figures/")

############################## a)
# Setup
gd = as.geodata(dt)
grid <- interp(dt$x, dt$y, dt$z)

# Simple visualization by interpolation
# pdf("2a_contour.pdf")
image.plot(grid, col=terrain.colors(200), asp=1,
           xlab="x", ylab="y"
)
contour(grid, add=T)
points(dt, pch=4, col="red", )
# dev.off()

# surface2a = s3dPlot(grid, fn="2a_surface.pdf")

# Linear regression to see linear data trend
# pdf("2a_linearPlane.pdf")
plane = lm(dt$z ~ dt$x + dt$y)
planeplot = scatterplot3d(dt$x,dt$y,dt$z, angle=-20,
                          xlab = "x", ylab = "y", zlab = "z")
planeplot$plane3d(plane)
# dev.off()

# pdf("2a_obsDensity.pdf")
heat <- as.ppp(cbind(dt$x,dt$y), c(0, 315, 0, 315))
plot(density(heat, adjust=0.5), col=viridis(200, direction=-1), main="")
# dev.off()


############################## b)

############################## c)
D = seq(1,315,1)
gridOrd = expand.grid(x = D, y = D)
# Spatial Prediction - Ordinary Kriging
spOK = krige.conv(
  geodata=gd,
  locations=gridOrd,
  krige = krige.control(
    # defaults to ordinary Kriging and constant mean
    cov.model = "powered.exponential",
    cov.pars = c(50^2, 100), # (partial sill, range)
    kappa = 1.5 # Smoothness
  )
)

predOK = matrix(spOK$predict, ncol = 315) # to fit image.plot structure
varOK = matrix(spOK$krige.var, ncol = 315)

# Ordinary Kriging predictor visualization
# pdf("2c_ordKrigPred.pdf")
predImageOK = list(x = D, y = D, z = predOK)
image.plot(predImageOK, col=terrain.colors(200),asp=1)
contour(predImageOK, add = T)
# dev.off()

# Associated prediction variance
# pdf("2c_ordKrigVar.pdf")
varImageOK = list(x = D, y = D, z = varOK)
image.plot(varImageOK, col=viridis(200),asp=1)
points(dt$x, dt$y, col="white", pch=4)
# dev.off()

# surface3d = s3dPlot(predImageOK, fn = "2c_ordKrigSurface.pdf")


############################## d)
gridUniv = expand.grid(x = D, y = D)
spUK     = krige.conv(
  geodata=gd,
  # coords    = cbind(x=x,y=y), 
  # data      = z, 
  locations = gridOrd,
  krige     = krige.control(
    type.krige = "ok",
    #"2nd": beta0 + beta1*x1 + beta2*x2 + beta3*(x1)^2 + beta4*(x2)^2 + beta5*x1*x2
    trend.d    = "2nd", 
    trend.l    = "2nd",
    cov.model  = "powered.exponential",
    cov.pars   = c(50^2, 100),
    kappa      = 1.5
  )
)

# Estimated beta by Universal Kriging
betaHat = (spUK$beta.est)
o2 = c(1,2,3,6,4,5) # order of betas corresponding to task description
betaHat[o2]

# Universal Kriging predictor
predUK = matrix(spUK$predict, ncol=315)
imageUK = list(x=D, y=D, z=predUK)
# pdf("2d_univKrigPred.pdf")
image.plot(imageUK,
           col=terrain.colors(200), 
           asp=1)
contour(list(x=D, y=D, z=predUK), add=T)
# points(dt$x,dt$y, pch=4, col = "red")
# dev.off()

# Associated variance predictions
varUK = matrix(spUK$krige.var,ncol = 315) # Predicted variance
# pdf("2d_univKrigVar.pdf")
image.plot(list(x=D,y=D,z=varUK), col=viridis(200), asp=1)
points(dt$x,dt$y, pch=4, col = "white")
# dev.off()
# surface3d2 = s3dPlot(imageUK, fn="2d_univKrigSurface.pdf")

# Difference between ordinary and universal
predDiff = (predUK - predOK)
# pdf("2d_predDiff.pdf")
image.plot(list(x=D,y=D,z=predDiff), col=turbo(200), asp=1)
# points(dt$x,dt$y, pch=4, col = "white")
# dev.off()

# Logarithmic color gradient for difference
predDiffLog = abs(predUK - predOK)
breaksPred = 1.1^(-190:10)/1.1^10*max(predDiff)
# pdf("2d_predDiffLog.pdf")
image.plot(list(x=D,y=D,z=predDiffLog), col=viridis(200), asp=1,breaks = breaksPred)
# points(dt$x,dt$y, pch=4, col = "white")
# dev.off()

varDiff = varUK - varOK
breaksVar = 1.1^(-189:10)/1.1^10*(max(varDiff)+1)
breaksVar = c(-0.001,breaksVar)
# pdf("2d_varDiff.pdf")
image.plot(list(x=D,y=D,z=varDiff), col=viridis(200), asp=1, breaks = breaksVar)
# points(dt$x,dt$y, pch=4, col = "white")
# dev.off()


############################## e)
vars0 = matrix(spOK$krige.var, ncol=315)[100,100]
x0 = matrix(spOK$predict, ncol = 315)[100,100]
c(krigPred = x0, predVar = vars0)
# P(X((100,100))>850)
pnorm(850, mean = x0, sd=sqrt(vars0), lower.tail = F)
# critical value for which the true val is below
# with a prob of 0.9
qnorm(0.9, x0, sqrt(vars0))