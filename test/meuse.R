# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)
# Load the data
data(meuse)
d <- cbind(meuse$x, meuse$y)
y <- 5*meuse$elev*meuse$x/max(meuse$x)
y2 <- meuse$cadmium*meuse$y/meuse$x
y3 <- meuse$lead/30

# Find anchorpoints
a <- find_anchorpoints.lsm(d,12)
# Build the empiric variogram
vario <- variogram.lsm(cbind(y,y2,y3,meuse$y/meuse$x*y3/y2),d,a$anchorpoints,570,4,15,dim = 1,kernel_id = "gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario , "exponential", c(570/2,571/2,0.1,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000),lower.bound = c(1e-8,1e-8,1e-8,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf))
# Plot of the solutions
solu$solutions

 vario1 <- variogram.lsm(y,d,a$anchorpoints,570,4,15,dim = 1,kernel_id = "gaussian")

solu <- findsolutions.lsm(vario1, "maternNuNugget_1", c(200,200,0.01,100,0.1))
##
x11()
mypoints<-plot.lsm(model = solu, a = a, z = c(y), d = d, n_points = 15, points_arrangement = "straight", kriging = TRUE, 
                   ellipse_scale = 2.5, arrow_scale = 1.5)

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d, plot_output = T)
max(previsions$zpredicted - y)

# Test the performace of our model via cross-validation
cv.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))
