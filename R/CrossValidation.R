#' Cross-Validation MSE
#' 
#' @description calculate the mean squared error via cross-validation
#' @param z the vector containing z(d)
#' @param d the matrix containing the coordinates in which we know the value of z
#' @param anchorpoints a matrix with the coordinates of the anchor points which can be obtained calling find_anchorpoints.lsm
#' @param epsilon the value of the bandwidth parameter epsilon
#' @param n_angles the number of angles for the grid
#' @param n_intervals the number of intervals for the grid
#' @param kernel_id the type of kernel to be used
#' @param id the type of variogram to be used
#' @param initial.position the starting position to be given to the optimizer
#' @param lower.bound the lower bound for the optimization, by default (1e-8, 1e-8, ...)
#' @param upper.bound the upper bound for the optimization, by default (Inf, Inf, pi/2, Inf, Inf, ...)
#' @param lower.delta set the minimum value for Cross-Validation search for optimal delta in smoothing equal to lowerdelta*epsilon
#' @param upper.delta set the maximum value for Cross-Validation search for optimal delta in smoothing equal to upperdelta*epsilon
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @examples 
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' cv.loo.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))
cv.loo.lsm <- function(z, d, anchorpoints, epsilon, n_angles, n_intervals, kernel_id, id, initial.position, lower.bound = rep(1e-8,length(initial.position)), upper.bound = c(c(Inf,Inf,pi/2), rep(Inf, length(initial.position)-3)), lower.delta = 0.1, upper.delta = 10, n_threads = -1){
  # set the MSE to 0
  MSE=0
  
  # create the progress bar
  pb <- txtProgressBar(min = 1,         # Minimum value of the progress bar
                       max = length(z), # Maximum value of the progress bar
                       style = 3,       # Progress bar style
                       width = 50,      # Progress bar width
                       char = "=")      # Character used to create the bar
  
  for(i in 1:length(z)){
    # create a new couple (znew, dnew) deleting the i-th element from (z, d)
    znew <- z[-i]
    dnew <- d[-i,]
    
    # predict the value of f(d[i, ]) and update the MSE
    vario <- variogram.lsm(znew,dnew,anchorpoints,epsilon,n_angles,n_intervals,kernel_id,FALSE, n_threads = n_threads)
    solu <- findsolutions.lsm(vario, id, initial.position, lower.bound, upper.bound, lower.delta, upper.delta, print=FALSE, n_threads = n_threads)
    previsions <- predict.lsm(solu, rbind(d[i,]),FALSE,FALSE, n_threads = n_threads)
    MSE <- MSE + (previsions$zpredicted - z[i])^2
    
    # update the progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  MSE=MSE/length(z)
  return(MSE)
}


#' Cross-Validation MSE for epsilon
#' 
#' @description calculate the best epsilon and mean squared error via cross-validation
#' @param z the vector (or matrix) containing z(d)
#' @param d the matrix containing the coordinates in which we know the value of z
#' @param anchorpoints a matrix with the coordinates of the anchor points which can be obtained calling find_anchorpoints.lsm
#' @param epsilons the values of the bandwidth parameter epsilon to test
#' @param n_angles the number of angles for the grid
#' @param n_intervals the number of intervals for the grid
#' @param kernel_id the type of kernel to be used
#' @param id the type of variogram to be used
#' @param initial.position the starting position to be given to the optimizer
#' @param lower.bound the lower bound for the optimization, by default (1e-8, 1e-8, ...)
#' @param upper.bound the upper bound for the optimization, by default (Inf, Inf, pi/2, Inf, Inf, ...)
#' @param lower.delta set the minimum value for Cross-Validation search for optimal delta in smoothing equal to lowerdelta*epsilon
#' @param upper.delta set the maximum value for Cross-Validation search for optimal delta in smoothing equal to upperdelta*epsilon
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @param sample_perc The percentage of points to be used as test set.
#' @examples 
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' cv.loo.lsm(y,d,a$anchorpoints,epsilons = c(350,360,370),8,8,"gaussian","exponential", c(200,200,0.01,100), sample_perc = 0.1)

cv.lsm <- function(z, d, anchorpoints, epsilons, n_angles, n_intervals, kernel_id, id, initial.position, lower.bound = rep(1e-8,length(initial.position)), upper.bound = c(c(Inf,Inf,pi/2) , rep(Inf, length(initial.position)-3)), dim = 1, sample_perc = 0.1, lower.delta = 0.1, upper.delta = 10, n_threads = -1){
  # set the MSE to 0
  MSE=numeric(length(epsilons))
  sample_idx = sample(1:dim(z)[1],sample_perc*dim(z)[1])
  
  # create the progress bar
  pb <- txtProgressBar(min = 1,         # Minimum value of the progress bar
                       max = length(epsilons), # Maximum value of the progress bar
                       style = 3,       # Progress bar style
                       width = 50,      # Progress bar width
                       char = "=")      # Character used to create the bar

  # create a new couple (znew, dnew) deleting the i-th element from (z, d)
  znew <- z[-sample_idx,]
  dnew <- d[-sample_idx,]
  min.MSE = 0
  for(i in 1:length(epsilons)){

    # predict the value
    vario <- variogram.lsm(z = znew,d = dnew,anchorpoints = a$anchorpoints,epsilon = epsilons[i],n_angles = n_angles,n_intervals = n_intervals,kernel_id = kernel_id,print_output = F, n_threads = n_threads, dim = dim)
    solu <- findsolutions.lsm(vario = vario, id = id, initial.position = initial.position, lower.bound = lower.bound, upper.bound = upper.bound, lower.delta = lower.delta, upper.delta = upper.delta,remove_not_convergent = T, print=FALSE, n_threads = n_threads)
    previsions <- predict.lsm(sol = solu, newpos = d[sample_idx,],plot_output = FALSE,print_output = FALSE, n_threads = n_threads, predict_y = T)
    
    # update the progress bar
    setTxtProgressBar(pb, i)
    MSE[i] <- mean((previsions$zpredicted - z[sample_idx,])^2)
    
    if(i == 1 || MSE[i] < min.MSE){
      min.MSE <- MSE[i]
      min.epsilon <- epsilons[i]
    }
    
  }
  close(pb)
  list <- list("MSE"=MSE,"min.epsilon" = min.epsilon, "min.MSE" = min.MSE)

  return(list)
}
