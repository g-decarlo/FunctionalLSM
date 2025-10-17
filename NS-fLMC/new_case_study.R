# 1. Setup ----
# Install and load necessary packages
# Ensure the custom package is installed from GitHub
# library(devtools)
# install_github("g-decarlo/FunctionalLSM")
library(LocallyStationaryModels)
# Use pacman to load/install packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  compositions,
  data.table,
  robCompositions,
  ggplot2,
  dplyr
)

# Set a seed for reproducibility
set.seed(42)

# ---
# 2. Data Loading and Preprocessing ----
# ---

# Load the datasets (adjust paths as needed)
# For this script to be reproducible, please replace the placeholder paths
# with the actual paths to your data files.
tryCatch({
  coords <- as.matrix(read.csv("~/FunctionalLSM/NS-fLMC/data/coordinatesrain.csv"))
  rain_obs <- read.csv("~/FunctionalLSM/NS-fLMC/data/rainobservations.csv")
  density_mat <- as.matrix(read.table("~/FunctionalLSM/NS-fLMC/data/density_matrix.prn", as.is = TRUE))
}, error = function(e) {
  stop("Data files not found. Please ensure 'coordinatesrain.csv', 'rainobservations.csv', and 'density_matrix.prn' are in a 'data/' subdirectory or adjust the paths.")
})


# Reshape and format the density matrix
colnames(density_mat) <- NULL
density_mat <- unname(density_mat)
density_mat <- matrix(density_mat, nrow = dim(coords)[1], ncol = 60)
mode(density_mat) <- "numeric"

# Smooth the density data using B-splines
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

# Perform PCA on the centered log-ratio (clr) transformed densities
pca_result <- prcomp(smoothed_splines$Y_clr, scale. = TRUE, center = TRUE)
z_scores <- pca_result$x

# ---
# 3. Data Splitting ----
# ---

# Split data into training and test sets
n_total <- nrow(coords)
n_test <- 4500
n_train <- n_total - n_test

cat(sprintf("Total observations: %d\n", n_total))
cat(sprintf("Training set size: %d\n", n_train))
cat(sprintf("Test set size: %d\n", n_test))

test_idx <- sample(1:n_total, size = n_test)
train_idx <- setdiff(1:n_total, test_idx)

coords_train <- coords[train_idx, ]
coords_test <- coords[test_idx, ]

# ---
# 4. Trace-NS Model Analysis (20 PCs) ----
# ---
cat("\n--- Starting Trace-NS Analysis (Top 20 PCs) ---\n")

# Select the top 20 PCA scores
n_pcs_trace <- 20
z_scores_trace_train <- z_scores[train_idx, 0:n_pcs_trace]
z_scores_trace_test <- z_scores[test_idx, 0:n_pcs_trace]

# Find anchor points using the training data coordinates
anchor_points <- find_anchorpoints.lsm(coords_train, 12, TRUE)

# Compute the empirical trace-variogram
variogram_trace <- variogram.lsm(
  z = z_scores_trace_train,
  d = coords_train,
  a = anchor_points$anchorpoints,
  epsilon = 10,
  n_angles = 6,
  n_intervals = 24,
  dim = 1, # dim=1 for trace-variogram
  kernel_id = "gaussian"
)

# Fit the non-stationary model to the trace-variogram
solution_trace <- findsolutions.lsm(
  variogram_trace,
  remove_not_convergent = TRUE,
  lower.delta = 0.5,
  upper.bound = c(10, 10, pi / 2, 20, 50),
  lower.bound = c(2, 2, 0, 1e-8, 1e-8),
  initial.position = c(10, 10, pi / 3, 10, 1),
  id = "exponentialnugget"
)

# Predict scores at the test locations
predictions_trace <- predict.lsm(solution_trace, coords_test, plot_output = FALSE)

# Calculate Mean Squared Prediction Error (MSPE) for Trace-NS
# MSPE is the mean of the sum of squared errors across all components
error_matrix_trace <- z_scores_trace_test - predictions_trace$zpredicted
squared_errors_per_location_trace <- rowSums(error_matrix_trace^2)
mspe_trace <- mean(squared_errors_per_location_trace)

cat(sprintf("MSPE for Trace-NS (20 PCs): %.4f\n", mspe_trace))

# ---
# 5. Simplified Op-NS Model Analysis (6 PCs) ----
# ---
cat("\n--- Starting Simplified Op-NS Analysis (Top 6 PCs, independently) ---\n")

# Select the top 6 PCA scores
n_pcs_op <- 6
z_scores_op_train <- z_scores[train_idx, 1:n_pcs_op]
z_scores_op_test <- z_scores[test_idx, 1:n_pcs_op]

# Matrix to store predictions for each PC
predictions_op_matrix <- matrix(NA, nrow = n_test, ncol = n_pcs_op)
colnames(predictions_op_matrix) <- paste0("PC", 1:n_pcs_op)

# Loop through each of the 6 principal components
for (k in 1:n_pcs_op) {
  cat(sprintf("Fitting model for PC %d...\n", k))
  
  # Isolate the k-th PC
  z_scalar_train <- z_scores_op_train[, k, drop = FALSE]
  
  # Compute the empirical variogram for this single component
  variogram_scalar <- variogram.lsm(
    z = z_scalar_train,
    d = coords_train,
    a = anchor_points$anchorpoints,
    epsilon = 10,
    n_angles = 6,
    n_intervals = 24,
    dim = 1,
    kernel_id = "gaussian"
  )
  
  # Fit the non-stationary model
  solution_scalar <- findsolutions.lsm(
    variogram_scalar,
    remove_not_convergent = TRUE,
    lower.delta = 0.5,
    upper.bound = c(10, 10, pi / 2, 20, 50),
    lower.bound = c(2, 2, 0, 1e-8, 1e-8),
    initial.position = c(10, 10, pi / 3, 10, 1),
    id = "exponentialnugget"
  )
  
  # Predict the k-th score at test locations
  predictions_scalar <- predict.lsm(solution_scalar, coords_test, plot_output = FALSE)
  
  # Store the predictions
  predictions_op_matrix[, k] <- predictions_scalar$zpredicted
}

# Calculate the combined MSPE for the simplified Op-NS
error_matrix_op <- z_scores_op_test - predictions_op_matrix
squared_errors_per_location_op <- rowSums(error_matrix_op^2)
mspe_op <- mean(squared_errors_per_location_op)

cat(sprintf("\nMSPE for Simplified Op-NS (6 PCs): %.4f\n", mspe_op))


# ---
# 6. Final Comparison ----
# ---
cat("\n--- Final Results ---\n")
cat(sprintf("Trace-NS MSPE (20 PCs): %.4f\n", mspe_trace))
cat(sprintf("Simplified Op-NS MSPE (6 PCs): %.4f\n", mspe_op))

if (mspe_trace < mspe_op) {
  cat("Conclusion: The Trace-NS model performed better.\n")
} else if (mspe_op < mspe_trace) {
  cat("Conclusion: The simplified Op-NS model performed better.\n")
} else {
  cat("Conclusion: Both models performed equally well.\n")
}