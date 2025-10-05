# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SCRIPT 1: Data Generation for Simulation Study
#
# Purpose: This script generates the computationally intensive spatial data
# for all three simulation scenarios and saves (caches) the output to disk.
# Run this script once before starting the analysis.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 1. Setup ---
library(LocallyStationaryModels)

# -----------------------------------------------------------------------------
# Parameter Definitions for Each Scenario
# -----------------------------------------------------------------------------

# Scenario 1: Non-stationary Proportional Data
params_nsprop <- function(coords) {
  x <- coords[, 1]
  y <- coords[, 2]
  r <- sqrt(x^2 + y^2)
  
  lambda1 <- (18 + 6 * x) / 100
  lambda2 <- (8 + 3 * x) / 100
  phi <- (pi / 6) * (1 + x)
  alpha_s <- 8 - 3 * r
  
  # Coregionalization matrix A^K
  A_base <- matrix(c(1, -2, 0, 2), nrow = 2)
  A11 <- alpha_s * A_base[1, 1]
  A21 <- alpha_s * A_base[2, 1]
  A22 <- alpha_s * A_base[2, 2]
  
  # Mean m^K
  m1 <- 1 + x
  m2 <- 2 + y
  
  return(data.frame(lambda1, lambda2, phi, A11, A21, A22, m1, m2))
}

# Scenario 2: Non-stationary Non-proportional Data
params_nsnonprop <- function(coords) {
  x <- coords[, 1]
  y <- coords[, 2]
  
  lambda1 <- (18 + 6 * x) / 100
  lambda2 <- (8 + 3 * x) / 100
  phi <- (pi / 6) * (1 + x)
  
  # Spatially varying coregionalization matrix A^K_s
  A11 <- 1.5 + x
  A21 <- 6 + 2 * y
  A22 <- 1.5 - x
  
  # Mean m^K
  m1 <- 1 + x
  m2 <- 2 + y
  
  return(data.frame(lambda1, lambda2, phi, A11, A21, A22, m1, m2))
}

# Scenario 3: Stationary Data
params_stat <- function(coords) {
  n_points <- nrow(coords)
  
  lambda1 <- rep(18 / 100, n_points)
  lambda2 <- rep(8 / 100, n_points)
  phi <- rep(pi / 6, n_points)
  
  # Constant coregionalization matrix A^K
  A11 <- rep(5, n_points)
  A21 <- rep(-10, n_points)
  A22 <- rep(10, n_points)
  
  # Zero mean
  m1 <- rep(0, n_points)
  m2 <- rep(0, n_points)
  
  return(data.frame(lambda1, lambda2, phi, A11, A21, A22, m1, m2))
}

# --- 2. Data Generation and Caching Function ---
generate_and_cache_data <- function(scenario_name, params_func) {
  cat(sprintf("Generating data for scenario: %s\n", scenario_name))
  set.seed(123) # for reproducibility
  grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 50)))
  
  # Get true parameters for the entire grid
  true_params <- params_func(grid_points)
  
  # Generate the spatial random field
  sim_data_list <- sample.lsm(
    d = grid_points,
    variogram_id = "exponential",
    parameters = as.matrix(true_params),
    dim = 2,
    n_samples = 20
  )
  sim_data <- sim_data_list$simulated_processes
  
  # Save the generated data and parameters to a file
  file_name <- paste0("scenario_", gsub(" ", "_", tolower(scenario_name)), ".RData")
  save(grid_points, true_params, sim_data, file = file_name)
  cat(sprintf("-> Data for '%s' saved to '%s'\n", scenario_name, file_name))
}

# --- 3. Main Execution ---
# Generate and save data for all three scenarios
generate_and_cache_data("Non-Stationary Proportional", params_nsprop)
generate_and_cache_data("Non-Stationary Non-Proportional", params_nsnonprop)
generate_and_cache_data("Stationary", params_stat)

cat("\nAll data generation complete.\n")
