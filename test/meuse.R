# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)


d <- cbind(simulated_process[,1],simulated_process[,2])
a <- find_anchorpoints.lsm(d,n = 6)
z <- as.matrix(simulated_process[,3])

vario <- variogram.lsm(z = z,d = d,a = a$anchorpoints,epsilon=0.2,n_angles=6,n_intervals=12,dim=1,kernel_id="gaussian")
plotvario(variogram = vario, 6)
solu <- findsolutions.lsm(vario ,lower.delta = 2, "exponentialnugget", c(0.2,0.2,0.1,0.3,0.03), upper.bound = c(0.4,0.4,pi/2,10,1))
solu$solutions
mypoints<-plot.lsm(model = solu, a=a, n_points = 5, points_arrangement = "straight", kriging = FALSE, 
                   ellipse_scale = 2.5 , arrow_scale = 1.5)



# Load the data
data(meuse)
d <- cbind(meuse$x, meuse$y)
y <- 5*meuse$elev*meuse$x/max(meuse$x)
y2 <- meuse$cadmium*meuse$y/meuse$x
y3 <- meuse$lead/30

# Find anchorpoints
a <- find_anchorpoints.lsm(d,12)
# Build the empiric variogram
vario <- variogram.lsm(cbind(y,y2),d,a$anchorpoints,570,4,15,dim = 2,kernel_id = "gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario ,lower.delta = 1, c("nugget","exp"), c(250,250,0.1,1000,1000,1000,250,250,0.1,100,100,100),lower.bound = c(1e-8,1e-8,1e-8,-Inf,-Inf,-Inf,1e-8,1e-8,1e-8,-Inf,-Inf,-Inf))
# Plot of the solutions
solu$solutions

 vario1 <- variogram.lsm(cbind(y,y2),d,a$anchorpoints,570,4,15,dim = 1,kernel_id = "gaussian")

solu <- findsolutions.lsm(vario1, "exponential", c(300,200,0.01,100,10), remove_not_convergent = T)
solu$solutions
plotvario(vario1,6)
##
x11()
mypoints<-plot.lsm(model = solu, a = a, n_points = 10, points_arrangement = "straight", kriging = TRUE, 
                   ellipse_scale = 2.5 , arrow_scale = 1.5)

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, a$anchorpoints, plot_output = F, predict_y = F)
max(abs(previsions$zpredicted - cbind(y,y2)))



