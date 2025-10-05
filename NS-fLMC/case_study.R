# 1. Setup ----
# Install and load necessary packages
library(devtools)
install_github("g-decarlo/FunctionalLSM")
library(LocallyStationaryModels)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  compositions,
  LocallyStationaryModels,
  data.table,
  robCompositions,
  usmap,
  ggplot2,
  cowplot,
  MASS
)

# Set a seed for reproducibility
set.seed(2023)

# Create a directory for the plots
if (!dir.exists("plots")) {
  dir.create("plots")
}

# 2. Data Loading and Preprocessing ----
# Load the datasets
coords <- as.matrix(read.csv("~/FunctionalLSM/NS-fLMC/data/coordinatesrain.csv"))
rain_obs <- read.csv("~/FunctionalLSM/NS-fLMC/data/rainobservations.csv")
density_mat <- as.matrix(read.table("~/FunctionalLSM/NS-fLMC/data/density_matrix.prn", as.is = TRUE))
colnames(density_mat) <- NULL
density_mat <- unname(density_mat)
density_mat <- matrix(density_mat, nrow = dim(coords)[1], ncol = 60)
mode(density_mat) <- "numeric"

# Probability of rain - CLR transformation
prob_rain <- rowMeans(rain_obs > 0)
prob_rain_comp <- cbind(prob_rain, 1 - prob_rain)
clr_rain <- clr(prob_rain_comp)[, 1, drop = FALSE]

# Smooth the density data
knots <- 10
xcp <- 1:60
data_range <- range(xcp)
knot_range <- c(data_range[1] - 0.1, data_range[2] + 0.1)
knots_vector <- seq(from = knot_range[1], to = knot_range[2], length.out = knots)
smoothed_splines <- smoothSplines(
  k = 3,
  l = 2,
  alpha = 0.98,
  xcp = as.numeric(xcp),
  data = density_mat,
  knots = knots_vector,
  num_points = 256,
  prior = "sq"
)

# Perform PCA and select components that explain most of the variance
pca_result <- prcomp(smoothed_splines$Y_clr, scale. = TRUE, center = TRUE)
variance_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
num_components <- which(variance_explained > 0.995)[1]
z_scores <- pca_result$x[, 1:num_components]

# Plot some density realizations
sample_idx <- sample(1:nrow(smoothed_splines$Y), 10)
matplot(
  seq(0, 60, length.out = 256),
  t(smoothed_splines$Y[sample_idx, ]),
  type = "l",
  lwd = 2,
  lty = 1,
  ylab = "Density",
  xlab = "Daily Precipitation (sqrt(mm))",
  main = "Sample of Smoothed Precipitation Densities"
)
ggsave("plots/sample_densities.png", width = 8, height = 6)


# 3. Non-Stationary Model Fitting for Densities ----
# Find anchor points
anchor_points <- find_anchorpoints.lsm(coords, 12, TRUE)
ggsave("plots/anchor_points.png", width = 8, height = 6)

# Compute the empirical variogram for the functional data
vario_trace <- variogram.lsm(
  z = z_scores,
  d = coords,
  a = anchor_points$anchorpoints,
  epsilon = 10,
  n_angles = 6,
  n_intervals = 24,
  dim = 1,
  kernel_id = "gaussian"
)

# Fit the non-stationary model
solu_trace <- findsolutions.lsm(
  vario_trace,
  remove_not_convergent = TRUE,
  lower.delta = 0.5,
  upper.bound = c(25, 25, pi / 2, 20, 200),
  lower.bound = c(2, 2, 0, 1e-8, 1e-8),
  initial.position = c(10, 10, pi / 3, 10, 50),
  id = "exponentialnugget",
)
solu_trace$solutions

# 4. Non-Stationary Model Fitting for Rain Probability ----
# Compute the empirical variogram for the CLR-transformed rain probability
vario_clr <- variogram.lsm(
  z = clr_rain,
  d = coords,
  a = anchor_points$anchorpoints,
  epsilon = 5,
  n_angles = 8,
  n_intervals = 16,
  dim = 1,
  kernel_id = "gaussian"
)

# Fit the non-stationary model
solu_clr <- findsolutions.lsm(
  vario_clr,
  remove_not_convergent = TRUE,
  lower.delta = 1,
  upper.bound = c(200, 200, pi / 2, 8, 200),
  lower.bound = c(2, 2, 0, 0, 1e-8),
  initial.position = c(2, 2, pi / 4, 6, 0.1),
  id = "exponentialnugget"
)
solu_clr$solutions

# 5. Spatially Varying Parameter Plots ----
# Generate plots of the estimated parameters
# These functions are assumed to be in your 'PlotFunctions.R'
# If they are not, you might need to write custom plotting functions
# using ggplot2 and the results from 'solu_trace' and 'solu_clr'

# For densities
smoothed_params_densities <- smooth.lsm(solu_trace, coords)
# For rain probability
smoothed_params_prob <- smooth.lsm(solu_clr, coords)

# Example of plotting one parameter (e.g., sigma)
# You'll need to adapt this for all your parameters
plot_df_densities <- data.frame(
  x = coords[, 1],
  y = coords[, 2],
  sigma = smoothed_params_densities$parameters[, 4]
)

p_sigma_densities <- ggplot(plot_df_densities, aes(x = x, y = y, color = sigma)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Spatially Varying Sigma (Densities)", color = "Sigma")

ggsave("plots/sigma_densities.png", plot = p_sigma_densities, width = 8, height = 6)

# Repeat for other parameters and for the rain probability model

# 6. Kriging and Prediction ----
# Define a grid for prediction
grid_to_predict <- as.matrix(expand.grid(
  x = seq(min(coords[, 1]), max(coords[, 1]), length.out = 100),
  y = seq(min(coords[, 2]), max(coords[, 2]), length.out = 100)
))


# Kriging for the functional data (z_scores)
predictions_densities <- predict.lsm(solu_trace, grid_to_predict, plot_output = FALSE)

# Kriging for the rain probability (clr_rain)
predictions_prob <- predict.lsm(solu_clr, grid_to_predict, plot_output = FALSE)

# 7. Final Hazard Map ----
# Here you would combine the predictions to create the final hazard map.
# This part is highly specific to your research question and may require
# back-transforming the CLR and PCA results.

# Example of plotting the predicted mean for the first principal component
plot_df_pred_densities <- data.frame(
  x = grid_to_predict[, 1],
  y = grid_to_predict[, 2],
  predicted_mean = predictions_densities$smoothed_means[, 1]
)

p_pred_densities <- ggplot(plot_df_pred_densities, aes(x = x, y = y, fill = predicted_mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Predicted Mean of 1st PC (Densities)", fill = "Mean")

ggsave("plots/predicted_mean_pc1.png", plot = p_pred_densities, width = 8, height = 6)


# Example of plotting the predicted rain probability
# Back-transform from CLR space
predicted_prob <- ilrInv(predictions_prob$predictedmean, orig = prob_rain_comp)[, 1]

plot_df_pred_prob <- data.frame(
  x = grid_to_predict[, 1],
  y = grid_to_predict[, 2],
  predicted_prob = predicted_prob
)

p_pred_prob <- ggplot(plot_df_pred_prob, aes(x = x, y = y, fill = predicted_prob)) +
  geom_raster() +
  scale_fill_viridis_c(labels = scales::percent) +
  theme_minimal() +
  labs(title = "Predicted Probability of Rain", fill = "Probability")

ggsave("plots/predicted_prob_rain.png", plot = p_pred_prob, width = 8, height = 6)