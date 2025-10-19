# SCRIPT 1: MODEL COMPARISON AND SELECTION
#
# PURPOSE:
# This script performs a systematic comparison of four spatial models
# (Trace-Stat, Op-Stat, Trace-NS, Op-NS)
# to identify the best-performing model for the U.S. rainfall dataset.
# The best model's configuration is then saved for the final analysis.

# ---
# 1. SETUP
# ---
# Required packages
library(devtools)
install_github("g-decarlo/FunctionalLSM", ref = "gdecarlo/case-study")
library(LocallyStationaryModels)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(compositions, data.table, robCompositions, ggplot2, dplyr)

# Ensure reproducibility
set.seed(42)

# ---
# 2. DATA LOADING AND PREPROCESSING
# ---
# Load data from the 'data' subdirectory.
tryCatch({
  coords <- as.matrix(read.csv("NS-fLMC/data/coordinatesrain.csv"))
  rain_obs <- read.csv("NS-fLMC/data/rainobservations.csv")
  density_mat <- as.matrix(read.table("NS-fLMC/data/density_matrix.prn", as.is = TRUE))
}, error = function(e) {
  stop("Data files not found. Please ensure they are in a 'data/' subdirectory.")
})

# Reshape and smooth density data
colnames(density_mat) <- NULL
density_mat <- matrix(unname(density_mat), nrow = dim(coords)[1], ncol = 60)
mode(density_mat) <- "numeric"

knots <- 10
xcp <- 1:60
knots_vector <- seq(from = range(xcp)[1] - 0.1, to = range(xcp)[2] + 0.1, length.out = knots)
smoothed_splines <- smoothSplines(
  k = 3, l = 2, alpha = 0.98, xcp = as.numeric(xcp), data = density_mat,
  knots = knots_vector, num_points = 256, prior = "sq"
)

# Perform PCA and get scores
pca_result <- prcomp(smoothed_splines$Y_clr, scale. = TRUE, center = TRUE)
z_scores <- pca_result$x

# ---
# 3. STATION FILTERING
# ---
station_filter_epsilon <- 0.2
n_initial <- nrow(coords)
indices_to_process <- 1:n_initial
indices_to_keep <- c()

cat(sprintf("Starting station filtering with epsilon = %.2f...\n", station_filter_epsilon))

while(length(indices_to_process) > 0) {
  current_index <- indices_to_process[1]
  indices_to_keep <- c(indices_to_keep, current_index)
  if (length(indices_to_process) > 1) {
    current_coords <- coords[current_index, , drop = FALSE]
    remaining_indices <- indices_to_process[-1]
    remaining_coords <- coords[remaining_indices, , drop = FALSE]
    distances <- sqrt(rowSums(sweep(remaining_coords, 2, current_coords, "-")^2))
    close_indices <- remaining_indices[which(distances <= station_filter_epsilon)]
    indices_to_remove_from_pool <- c(current_index, close_indices)
  } else {
    indices_to_remove_from_pool <- current_index
  }
  indices_to_process <- setdiff(indices_to_process, indices_to_remove_from_pool)
}

cat(sprintf("Finished filtering. %d stations remaining.\n", length(indices_to_keep)))

coords <- coords[indices_to_keep, ]
z_scores <- z_scores[indices_to_keep, ]

# ---
# 4. MODEL COMPARISON SETUP
# ---
# Dynamic selection of PCs to explain >99.999% of variance
variance_explained <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
n_pcs_trace <- which(variance_explained > 0.99999)[1]
cat(sprintf("Using %d components for Trace models (explaining >99.999%% of variance).\n", n_pcs_trace))

# Define train/test split
n_total <- nrow(coords)
test_percentage <- 0.30
n_test <- floor(n_total * test_percentage)
n_train <- n_total - n_test

cat(sprintf("Splitting data: %d training, %d test observations.\n", n_train, n_test))

test_idx <- sample(1:n_total, size = n_test)
train_idx <- setdiff(1:n_total, test_idx)

coords_train <- coords[train_idx, ]
coords_test <- coords[test_idx, ]

z_scores_trace_train <- z_scores[train_idx, 1:n_pcs_trace]
z_scores_trace_test <- z_scores[test_idx, 1:n_pcs_trace]

# ---
# 5. STATIONARY MODEL ANALYSIS
# ---
cat("\n--- STARTING STATIONARY MODEL ANALYSIS ---\n")
anchor_points_stat <- find_anchorpoints.lsm(coords_train, 1, TRUE)

# 5.1 Trace-Stationary Model
cat("Running Trace-Stationary model...\n")
variogram_trace_stat <- variogram.lsm(
  z = z_scores_trace_train, d = coords_train, a = anchor_points_stat$anchorpoints,
  n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "Identity", epsilon = 50
)
solution_trace_stat <- findsolutions.lsm(
  variogram_trace_stat, remove_not_convergent = TRUE, lower.delta = 0.5,
  upper.bound = c(50, 50, pi / 2, 40, 100), lower.bound = c(1, 1, 0, 1e-8, 1e-8),
  initial.position = c(1, 1, pi / 3, 30, 1), id = "exponentialnugget"
)
predictions_trace_stat <- predict.lsm(solution_trace_stat, coords_test, plot_output = FALSE)
mspe_trace_stat <- mean(rowSums((z_scores_trace_test - predictions_trace_stat$zpredicted)^2))

if (mspe_trace_stat > 300) {
  cat("-> Trace-Stationary MSPE: Diverged (>300)\n")
  mspe_trace_stat <- Inf
} else {
  cat(sprintf("-> Trace-Stationary MSPE: %.4f\n", mspe_trace_stat))
}

# 5.2 Op-Stationary Model
cat("Running Op-Stationary model...\n")
n_pcs_op <- 6
z_scores_op_train <- z_scores[train_idx, 1:n_pcs_op]
predictions_op_matrix_stat <- matrix(NA, nrow = n_test, ncol = n_pcs_op)

for (k in 1:n_pcs_op) {
  z_scalar_train <- z_scores_op_train[, k, drop = FALSE]
  variogram_scalar_stat <- variogram.lsm(
    z = z_scalar_train, d = coords_train, a = anchor_points_stat$anchorpoints,
    n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "Identity", epsilon = 50
  )
  solution_scalar_stat <- findsolutions.lsm(
    variogram_scalar_stat, remove_not_convergent = TRUE, lower.delta = 0.5,
    upper.bound = c(100, 100, pi / 2, 200, 500), lower.bound = c(2, 2, 0, 1e-8, 1e-8),
    initial.position = c(50, 50, pi / 3, 20, 100), id = "exponentialnugget"
  )
  predictions_scalar_stat <- predict.lsm(solution_scalar_stat, coords_test, plot_output = FALSE)
  predictions_op_matrix_stat[, k] <- predictions_scalar_stat$zpredicted
}

# Calculate MSPE for Op-Stat across all original components
full_predictions_op_stat <- matrix(0, nrow = n_test, ncol = n_pcs_trace)
full_predictions_op_stat[, 1:n_pcs_op] <- predictions_op_matrix_stat
error_matrix_op_stat <- z_scores_trace_test - full_predictions_op_stat
mspe_op_stat <- mean(rowSums(error_matrix_op_stat^2))

if (mspe_op_stat > 300) {
  cat("-> Op-Stationary MSPE: Diverged (>300)\n")
  mspe_op_stat <- Inf
} else {
  cat(sprintf("-> Op-Stationary MSPE: %.4f\n", mspe_op_stat))
}

# ---
# 6. NON-STATIONARY MODEL ANALYSIS (EPSILON TUNING)
# ---
cat("\n--- STARTING NON-STATIONARY MODEL ANALYSIS ---\n")
anchor_points_ns <- find_anchorpoints.lsm(coords_train, 12, TRUE)
epsilon_values <- c(5, 8, 10, 15, 20)
results_df <- data.frame()
solution_list <- list()

for (current_epsilon in epsilon_values) {
  cat(sprintf("\n--- Testing Epsilon = %d ---\n", current_epsilon))
  
  # 6.1 Trace-NS Model
  variogram_trace <- variogram.lsm(
    z = z_scores_trace_train, d = coords_train, a = anchor_points_ns$anchorpoints,
    epsilon = current_epsilon, n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "gaussian"
  )
  solution_trace <- findsolutions.lsm(
    variogram_trace, remove_not_convergent = TRUE, lower.delta = 0.5,
    upper.bound = c(current_epsilon, current_epsilon, pi / 2, 20, 50), lower.bound = c(2, 2, 0, 1e-8, 1e-8),
    initial.position = c(10, 10, pi / 3, 10, 1), id = "exponentialnugget"
  )
  solution_list[[as.character(current_epsilon)]] <- solution_trace
  
  predictions_trace <- predict.lsm(solution_trace, coords_test, plot_output = FALSE)
  mspe_trace <- mean(rowSums((z_scores_trace_test - predictions_trace$zpredicted)^2))
  
  if (mspe_trace > 300) {
    cat("-> Trace-NS MSPE: Diverged (>300)\n")
    mspe_trace <- Inf
  } else {
    cat(sprintf("-> Trace-NS MSPE: %.4f\n", mspe_trace))
  }
  
  # 6.2 Op-NS Model
  predictions_op_matrix <- matrix(NA, nrow = n_test, ncol = n_pcs_op)
  for (k in 1:n_pcs_op) {
    z_scalar_train <- z_scores_op_train[, k, drop = FALSE]
    variogram_scalar <- variogram.lsm(
      z = z_scalar_train, d = coords_train, a = anchor_points_ns$anchorpoints,
      epsilon = current_epsilon, n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "gaussian"
    )
    solution_scalar <- findsolutions.lsm(
      variogram_scalar, remove_not_convergent = TRUE, lower.delta = 0.5,
      upper.bound = c(current_epsilon, current_epsilon, pi / 2, 20, 50), lower.bound = c(2, 2, 0, 1e-8, 1e-8),
      initial.position = c(10, 10, pi / 3, 10, 1), id = "exponentialnugget"
    )
    predictions_scalar <- predict.lsm(solution_scalar, coords_test, plot_output = FALSE)
    predictions_op_matrix[, k] <- predictions_scalar$zpredicted
  }
  
  # Calculate MSPE for Op-NS across all original components
  full_predictions_op <- matrix(0, nrow = n_test, ncol = n_pcs_trace)
  full_predictions_op[, 1:n_pcs_op] <- predictions_op_matrix
  mspe_op <- mean(rowSums((z_scores_trace_test - full_predictions_op)^2))
  
  if (mspe_op > 300) {
    cat("-> Op-NS MSPE: Diverged (>300)\n")
    mspe_op <- Inf
  } else {
    cat(sprintf("-> Op-NS MSPE: %.4f\n", mspe_op))
  }
  
  results_df <- rbind(results_df, data.frame(
    epsilon = current_epsilon, trace_mspe = mspe_trace, op_mspe = mspe_op
  ))
}

# ---
# 7. FINAL COMPARISON & SAVING CHAMPION MODEL
# ---
cat("\n\n--- FINAL RESULTS SUMMARY ---\n")

# Create plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Save detailed results to a CSV file
write.csv(results_df, "plots/model_comparison_results.csv", row.names = FALSE)

# Find the best result for each non-stationary method
best_trace_ns_row <- results_df[which.min(results_df$trace_mspe), ]
best_op_ns_row <- results_df[which.min(results_df$op_mspe), ]

# Print the full summary
cat("--- Stationary Models ---\n")
cat(sprintf("Trace-Stat: MSPE = %.4f\n", mspe_trace_stat))
cat(sprintf("Op-Stat:    MSPE = %.4f\n", mspe_op_stat))
cat("\n--- Best Non-Stationary Models ---\n")
cat(sprintf("Best Trace-NS: Epsilon = %d, MSPE = %.4f\n", best_trace_ns_row$epsilon, best_trace_ns_row$trace_mspe))
cat(sprintf("Best Op-NS:    Epsilon = %d, MSPE = %.4f\n", best_op_ns_row$epsilon, best_op_ns_row$op_mspe))

# Determine the overall best model
all_results <- data.frame(
  model = c("Trace-Stat", "Op-Stat", "Trace-NS", "Op-NS"),
  mspe = c(mspe_trace_stat, mspe_op_stat, best_trace_ns_row$trace_mspe, best_op_ns_row$op_mspe)
)
best_model_info <- all_results[which.min(all_results$mspe), ]

cat(sprintf("\n--- Overall Winner: %s with an MSPE of %.4f ---\n", best_model_info$model, best_model_info$mspe))

# Save the champion model object for the plotting script
if (best_model_info$model == "Trace-NS") {
  champion_epsilon <- best_trace_ns_row$epsilon
  champion_solution <- solution_list[[as.character(champion_epsilon)]]
  cat(sprintf("\n--- SAVING CHAMPION MODEL (Trace-NS, epsilon = %d) ---\n", champion_epsilon))
  save(champion_solution, pca_result, file = "champion_model.RData")
} else {
  cat("\nChampion model was not Trace-NS. The plotting script requires a saved Trace-NS model.\n")
}

