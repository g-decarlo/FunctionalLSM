# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)

data(SwissData)

d <- SwissData[,1:2]
z <- SwissData[,3]

# Find anchorpoints
a <- find_anchorpoints.lsm(d,12)
# Build the empiric variogram
vario <- variogram.lsm(z,d,a$anchorpoints,46,6,15,"gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario, "exponential", c(30,30,0.5,100),upper.bound = c(46,46,pi/2,250))
solu
plotvario(vario,3)
# Plot of the solutions
x11()
mypoints<-plot.lsm(model = solu, a = a, z = y, d = d, n_points = 5, points_arrangement = "straight",ellipse_scale = 1.5, arrow_scale =  2)

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d)
max(previsions$zpredicted - z)
