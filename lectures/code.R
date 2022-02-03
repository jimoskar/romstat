####################
### Package list ### 
####################

# mvrnorm
library(MASS)

# image.plot
library(fields)
# interp
library(akima)
# kriging
library(geoR)
# Plotting package
library(ggplot2)
# Color package
library(viridis)

# seed
set.seed(8080)
# same seed across operative systems
# set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion")

# ?function_name - gives function details in help window, useful to see what arguments function takes
# variable name autofill tab-button

### Variable list ###



# Starting out with NA, could be useful to indiciate if any operations went wrong and not all
# of the matrix has been filled out with values for example, since some calculations might produce
# zeros as valid values, so if we start out with all zeros, this might be tricky
data_struct <- matrix(NA, nrow = 6, ncol = 10)

for(i in 1:ncol(data_struct)){
  
  # Filling out matrices column-wise
  data_struct[,i] <- rnorm(nrow(data_struct),0,1)
  # Similarly row-wise
  # data_struct[j,] <- rnorm(ncol(data_struct),0,1), where j 1:nrow(data_struct)
  
}

x_coord <- seq(1,10,1)
# x_coord <- 1:10
# x_coord <- as.vector(1:10)
y_coord <- seq(1,10,length.out =10)
z_value <- rnorm(10,0,1)
z_value_new <- matrix(rnorm(100,0,1), nrow=10, ncol=10)

coord_matrix <- cbind(x_coord,y_coord,z_value)


data_struc_df <- as.data.frame(data_struct)
coord_df <- as.data.frame(coord_matrix)

diag_mx <- diag(1:10)

# multiplication test - two ways two multiple two structures, f ex two matrices, a vector and a matrix etc
# elementwise

elem_multiplic <- diag_mx*z_value


# matrix multiplication

multip_mx <- diag_mx%*%z_value

# inverting the diagonal matrix
inv_diag <- solve(diag_mx)

# transpose

mx_transp <- t(elem_multiplic)


# Visualizatio

## base R and image.plot
plot(data_struct[1,], xlab="Entry number", ylab="Value", main="Scatterplot of the first matrix", pch=20, col="magenta")
points(data_struct[2,],xlab="Entry number", ylab="Value", main="Scatterplot of the first matrix", pch=20, col="blue")
points(data_struct[3,],xlab="Entry number", ylab="Value", main="Scatterplot of the first matrix", pch=20, col="orange")

# bty - points, pch also for points, lty for line type
# locations either given as topright, topleft etc or coordinates
legend("topright", legend=c("Row 1", "Row 2", "Row 3"), pch = c(20,20,20),
              col=c("magenta", "blue", "orange"), bty="o", cex=0.8,
              title="Matrix rows", text.font=4)

image.plot(elem_multiplic, legend.lab = "Level", main="Your favourite matrix", xlab = "index1", ylab="index2")

## ggplot



# Loading data into R

a <- read.table('filename.fileext')

# generating mvnormal smaple with built in function and using cholesky factorization
# need a mean and a covariance matrix
grid_1d <- 1:10
# making covariance matrices

# expand.grid - making a grid structure a bit like coordinates

# matern() can be used for the matern function, or specifying matern instead of powered.exponential; cov.pars=c(sigma^2,phi)
cov_mx <- cov.spatial(as.matrix(dist(expand.grid(grid_1d))), cov.model="powered.exponential", cov.pars=c(1,5), kappa=1.5)


sample_1 <- mvrnorm(1,z_value,cov_mx)

# using cholesky factorization
L <- chol(cov_mx)
sample_ch <- t(L)%*%rnorm(length(z_value),0,1)+z_value



data(akima)


image.plot(interp(akima$x, akima$y, akima$z), asp=1)
interp_surface <-interp(akima$x, akima$y, akima$z)
contour(interp_surface$x, interp_surface$y, interp_surface$z)

akima_data <- as.data.frame(akima)
# kriging
temp <- as.geodata(akima_data)

# here you can use this method, or just use rdist
locs <- expand.grid(1:50,1:50)
#krige.conv()

#krige.control()

krigPred <- krige.conv(temp,locations = locs, krige = 
                        krige.control(type.krige = "ok", trend.d = "2nd", trend.l = "2nd",
                                      cov.pars = c(750, 150)))
prediction <- krigPred$predict

image.plot(interp(locs$Var1,locs$Var2,prediction))


# Plotting with ggplot

interp_df <- as.data.frame(cbind(locs,prediction))

fig = ggplot(data = interp_df, aes(Var1, Var2)) +
  geom_raster(aes(fill = prediction)) +
  scale_fill_viridis_c(name="Prediction") +
  coord_fixed() + xlab("Variable 1") + ylab("Variable 2")+
  ggtitle("Predicted surface")
fig

# Scatter plot
# picking out every third element, just to have less grid points in order to illustrate the scatterplot
interp_df_sc <- as.data.frame(cbind(locs$Var1[seq(1,length(locs$Var1),3)],locs$Var2[seq(1,length(locs$Var1),3)],prediction[seq(1,length(locs$Var1),3)]))

# similar can be done with fields package and quilt.plot

ggplot(data = interp_df_sc) + geom_point(aes(V1,V2,color = V3)) +
  scale_color_viridis_c(name="scale") + ggtitle("Scatter plot")

# for latex-expressions in titles or labels look up combination of paste and expression
ggplot(data = interp_df_sc) + geom_point(aes(V1,V2,color = V3)) +
  scale_color_viridis_c(name="scale") + ggtitle(expression(paste("With symbols in title","\n", sigma==0.05,", ", phi==35)))


# for base R a similar scatterplot for irregular spatial data
quilt.plot(interp_df_sc)

