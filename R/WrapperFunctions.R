#' Find Solutions (Fit Variogram)
#' 
#' @description for each anchorpoints solves a problem of nonlinear optimization and returns the results
#' @param vario a "sample_variogram" object obtained using variogram.lsm()
#' @param id the type of variogram to be used. Can be one of the following: "exponential", "gaussian", "matern", "maternNuFixed N", where N must
#' be replaced with a double of your choice. If the id is not correct, "exponential" is used by default. Remember that the matern has an extra
#' parameter nu that must be passed to the function.
#' @param initial.position the starting position to be given to the optimizer
#' @param lower.bound the lower bound for the optimization, by default (1e-8, 1e-8, ...)
#' @param upper.bound the upper bound for the optimizaion, by default (Inf, Inf, pi/2, Inf, Inf, ...)
#' @param lower.delta set the minimum value for Cross-Validation search for optimal delta in smoothing equal to lowerdelta*epsilon
#' @param upper.delta set the maximum value for Cross-Validation search for optimal delta in smoothing equal to upperdelta*epsilon
#' @param remove_not_convergent if set to TRUE removes the anchorpoints which cause troubles to the optimizer, by default is TRUE
#' @param print_output if set to FALSE suppress the console output, by default is TRUE
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @return an object containing the matrix with the optimal parameters, the optimal value of delta for smoothing, the value of the bandwidth
#' parameter epsilon used, the matrix with the coordinates of the anchor points used, the id of the variogram chosen and the id of the kernel
#' you used to generate the sample variogram
#' @details given an object of type "sample_variogram" returned by variogram.lsm, this function solves a problem of non linear
#' optimization in order to find the parameters that better fit the variogram function chosen via id in each anchor point. The initial position
#' to find the optimum must be provided by the user, which can also provide the upper and lower bounds for the solutions. Always remember that
#' the order of the parameters is always lambda1, lambda2, phi and sigma followed by additional ones required by the chosen variogram function
#' @examples
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' vario <- variogram.lsm(y,d,a$anchorpoints,370,8,8,"gaussian")
#' solu <- findsolutions.lsm(vario, "exponential", c(200,200,0.01,100))
findsolutions.lsm<-function(vario, id, initial.position, lower.bound = rep(1e-8,length(initial.position)), upper.bound = c(c(Inf,Inf,pi/2), rep(Inf, length(initial.position)-3)), lower.delta = 0.1, upper.delta = 10, remove_not_convergent = FALSE, print_output = TRUE, n_threads = -1)
{ 
  if(grepl("maternNuFixed", id, fixed = TRUE))
  {
    id_check <- "maternNuFixed"
  }
  else if(grepl("maternNuNugget", id, fixed = TRUE))
  {
    id_check <- "maternNuNugget"
  }
  else
  {
    id_check <- id
  }
  #if(vario$dim == 1 &( length(initial.position) != variogramfunctions$n_parameters[which(variogramfunctions$name == id_check)] || length(lower.bound) != variogramfunctions$n_parameters[which(variogramfunctions$name == id_check)] || length(upper.bound) != variogramfunctions$n_parameters[which(variogramfunctions$name == id_check)]))
  #{
   # stop("wrong number of initial parameters")#change to stop
  #}
  result <- findsolutionslsm(vario$anchorpoints, vario$empiricvariogram, vario$squaredweigths, vario$dim, vario$mean.x, vario$mean.y, id, vario$kernel_id, initial.position, lower.bound, upper.bound, vario$epsilon, lower.delta, upper.delta, print_output, n_threads)
  if (remove_not_convergent)
  {
    for (i in 1:dim(result$solutions)[1])
    {
      if (norm((result$solutions[i, ]-initial.position), type="2") < 1e-12)
      {
        result$solutions <- result$solutions[-i,]
        result$anchorpoints <- result$anchorpoints[-i, ]
      }
    }
  }
  result$id <- id
  result$kernel_id <- vario$kernel_id
  result$initial_coordinates <- vario$initial_coordinates
  result$initial_z <- vario$initial_z
  result$dim <- vario$dim
  class(result) <- "lsm"
  return(result)
}

#' Predict LSM (Kriging)
#' 
#' @description for each couple of coordinates in newpos predict the mean and punctual value of z
#' @param sol an object of type lsm obtained by calling findsolutions.lsm
#' @param newpos a matrix with the coordinates of the points where to evaluate z
#' @param plot_output if set to TRUE plot the solutions, by default is TRUE
#' @param print_output if set to FALSE suppress the console output, by default is TRUE
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @return an object containing the vector with the means, the vector with the punctual predictions and the vector with the kriging variance
#' in newpos
#' @details given an object of type "lsm" returned by findsolutions.lsm, this function performs kriging on the coordinates provided by newpos
#' and possibly plot the results found
#' @examples
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' vario <- variogram.lsm(y,d,a$anchorpoints,370,8,8,"gaussian")
#' solu <- findsolutions.lsm(vario, "exponential", c(200,200,0.01,100))
#' previsions <- predict.lsm(solu, d)
predict.lsm<-function(sol, newpos, plot_output = TRUE, print_output = TRUE, n_threads = -1)
{
  d <- sol$initial_coordinates
  z <- sol$initial_z
  predictedvalues <- predikt(as.matrix(z),d,sol$anchorpoints,sol$epsilon,sol$delta,sol$dim,sol$solutions,newpos,sol$id,sol$kernel_id,print_output,n_threads)
  if (plot_output)
  {
    par(ask=TRUE)
    for(i in 1:dim(predictedvalues$zpredicted)[2]){
      newpos <- as.data.frame(newpos)
      colnames(newpos) <- c("X", "Y")
      means <- ggplot2::ggplot(newpos, ggplot2::aes(x=X, y=Y, color=predictedvalues$predictedmean[,i])) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
      ys <- ggplot2::ggplot(newpos, ggplot2::aes(x=X, y=Y, color=predictedvalues$zpredicted[,i])) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
      means<-means+ggplot2::labs(color="mean") + ggplot2::theme_light()
      ys<-ys+ggplot2::labs(color="z") + ggplot2::theme_light()
      title <- cowplot::ggdraw() + cowplot::draw_label(paste("Predicted mean and z -",as.character(i)), fontface='bold')
      p <- cowplot::plot_grid(means, ys)
      print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
    }
    par(ask=FALSE)
  }
  return(predictedvalues)
}

#' Find Anchor Points
#' 
#' @description given a dataset find the corresponding equally spaced anchorpoints
#' @param dataset a dataset full of coordinates
#' @param n a parameter proportional to the number of anchorpoints
#' @param plot_output if set to true plot the original dataset and the anchorpoints
#' @return the coordinates of the anchor points found and the specification of the grid used to build them
#' @details given a set of points this function builds n*n grid covering all the locations. Then takes as anchor points all the
#' centers of the cells such that at least one of the points of d belongs to the same cell 
#' @examples
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12)
find_anchorpoints.lsm<-function(dataset, n, plot_output = TRUE)
{
  min.c <- rep(0., dim(dataset)[2])
  for (i in 1:dim(dataset)[2])
  {
    min.c[i] <- abs(min(dataset[,i]))
    dataset[, i] <- dataset[, i] + min.c[i] + 1
  }
  result <- find_anchorpoints(dataset, n)
  for (i in 1:dim(dataset)[2])
  {
    dataset[, i] <- dataset[, i] - min.c[i] -1
    result$anchorpoints[, i] <- result$anchorpoints[, i] - min.c[i] - 1
  }
  if (plot_output)
  {
    D <- as.data.frame(dataset)
    A <- as.data.frame(result)
    colnames(D) <- c("X", "Y")
    colnames(A) <- c("X", "Y")
    p <- ggplot2::ggplot(data = D, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_point() + ggplot2::geom_point(data = A, ggplot2::aes(x=X, y=Y, color="red"), shape = 3) + ggplot2::theme_light() + ggplot2::theme(legend.position = "none") + ggplot2::coord_fixed()
    print(p)
  }
  return(result)
}

#' Empiric Variogram LSM
#' 
#' @description compute the sample variogram in the anchorpoints
#' @param z the vector contatining f(d)
#' @param d the matrix contatining the coordinates in which we know the value of z
#' @param anchorpoints a matrix with the coordinates of the anchorpoints which can be obtained calling find_anchorpoints.lsm
#' @param epsilon the value of epsilon regulating the kernel
#' @param n_angles the number of angles for the grid
#' @param n_intervals the number of intervals for the grid
#' @param kernel_id the type of kernel to be used. At the moment the only possibility is "gaussian".
#' @param print_output if set to FALSE suppress the console output, by default is TRUE
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @return an object of type "sample_variogram" containing the kernel matrix, the grid matrix, the vecors with the value of x, y and norm of
#' every tile of the grid, the matrix of the squaredweights, the matrix with the sample variogam, the matrix with the anchor points used,
#' the value of the bandwidth parameter epsilon used, the id of the kernel function, the number of angles and of intervals used to build
#' the grid, the matrix with the coordinates of the initial points and the vector with the function z evaluated in these points.
#' @details the purpose of this function is to calculate the value of the sample variogram in every anchor point. To do so 
#' the function requires to be given as input all the information about the construction of the grid and of the kernel as in the paper by
#' Fouedjio.
#' @examples
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' vario <- variogram.lsm(y,d,a$anchorpoints,370,8,8,"gaussian")
variogram.lsm <- function(z, d, anchorpoints, epsilon, n_angles, n_intervals, dim = 1, kernel_id, print_output=TRUE, n_threads = -1)
{ 
  z=as.matrix(z)
  if(dim(z)[1] != dim(d)[1])
  {
    print("The length of z and the number or rows of d do not coincide")
  }
  if(dim(z)[2] != dim && dim != 1)
  {
    print("Dimension automatically set to 1")
    dim <- 1
  }
  vario <- variogramlsm(z, d, anchorpoints, epsilon, n_angles, n_intervals, dim, kernel_id, print_output, n_threads)
  vario$kernel_id <- kernel_id
  vario$n_angles <- n_angles
  vario$n_intervals <- n_intervals
  vario$initial_coordinates <- d
  vario$initial_z <- z
  vario$dim <- dim
  class(vario) <- "sample_variogram"
  return(vario)
}

#' Smooth LSM
#' 
#' @description compute the value of the parameters of the variogram function of model in the points contained in newpoints
#' @param model a "lsm" object generated via findsolutions.lsm
#' @param newpoints a matrix with the coordinates of the points the knowledge of the parameters is needed
#' @param n_threads the number of threads for OpenMP, by default is equal to -1, which means that OpenMP will use all the available threads.
#' @return a matrix with the values of the paramters smoothed in newpoints
#' @details given model, this function exploits model$solutions and model$delta to perform smoothing and find the value of the 
#' parameters regulating the variogram function in other points beyond the anchor ones. model$delta already contains the optimal value of 
#' delta which does not need to be evaluated again.
#' @examples 
#' data(meuse)
#' d <- cbind(meuse$x, meuse$y)
#' y <- meuse$elev
#' a <- find_anchorpoints.lsm(d,12,FALSE)
#' vario <- variogram.lsm(y,d,a$anchorpoints,370,8,8,"gaussian")
#' solu <- findsolutions.lsm(vario, "exponential", c(200,200,0.01,100))
#' newparams <- smooth.lsm(solu, d)
smooth.lsm <- function(model, newpoints, n_threads = -1)
{
  result <- smoothing(model$solutions,model$anchorpoints,model$delta,newpoints,model$kernel_id,n_threads)
  return(result)
}
