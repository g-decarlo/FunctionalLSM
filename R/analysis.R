# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SCRIPT: Full Simulation Study Analysis
#
# Purpose: This script analyzes pre-generated data for the simulation study
# presented in "Fully Non-Stationary Functional Spatial Random Fields".
# It implements all four comparison methods and uses the spatial average of
# the true parameters as starting points for optimization.
#
# Workflow:
# 1. Loads cached data for a given scenario.
# 2. If data is not found, it stops with an error.
# 3. Runs the analysis loop for all methods and repetitions.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 1. Setup ---
# Install required packages if they are not already installed
if (!require("LocallyStationaryModels")) {
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("elucasticus/FunctionalLSM-main")
}
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")


library(LocallyStationaryModels)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Parameter Definitions for Each Scenario
# (These functions define the true underlying models)
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
  
  A_base <- matrix(c(1, -2, 0, 2), nrow = 2)
  A11 <- alpha_s * A_base[1, 1]
  A21 <- alpha_s * A_base[2, 1]
  A22 <- alpha_s * A_base[2, 2]
  
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
  
  A11 <- 1.5 + x
  A21 <- 6 + 2 * y
  A22 <- 1.5 - x
  
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
  
  A11 <- rep(5, n_points)
  A21 <- rep(-10, n_points)
  A22 <- rep(10, n_points)
  
  m1 <- rep(0, n_points)
  m2 <- rep(0, n_points)
  
  return(data.frame(lambda1, lambda2, phi, A11, A21, A22, m1, m2))
}


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Calculates Mean Squared Error (MSE)
calculate_mse <- function(true_values, estimated_values) {
  if (is.vector(true_values)) {
    return(mean((true_values - estimated_values)^2, na.rm = TRUE))
  } else {
    return(mean(rowSums((true_values - estimated_values)^2), na.rm = TRUE))
  }
}

# Main function to run the simulation for a given scenario
run_simulation <- function(scenario_name, params_func, N_values, M_repetitions) {
  
  # --- Data Loading Logic ---
  file_name <- paste0("scenario_", gsub(" ", "_", tolower(scenario_name)), ".RData")
  if (file.exists(file_name)) {
    cat(sprintf("Loading cached data for scenario: %s\n", scenario_name))
    load(file_name)
  } else {
    stop(sprintf("Cached data file '%s' not found. Please run the data generation script first.", file_name))
  }
  
  N_tot <- nrow(grid_points)
  results_df <- data.frame()
  
  # --- Define Parameter Names and Optimization Bounds ---
  param_names_op <- c("lambda1", "lambda2", "phi", "A11", "A21", "A22")
  param_names_trace <- c("lambda1", "lambda2", "phi", "sigma")
  
  # Calculate bounds for Operator methods
  lower_bounds_op <- sapply(true_params[, param_names_op], min)
  upper_bounds_op <- sapply(true_params[, param_names_op], max)
  padding_op <- (upper_bounds_op - lower_bounds_op) * 0.1
  padding_op[padding_op == 0] <- 1
  lower_bounds_op <- lower_bounds_op - padding_op
  upper_bounds_op <- upper_bounds_op + padding_op
  
  # Calculate bounds for Trace methods
  true_sigma <- sqrt(true_params$A11^2 + true_params$A21^2 + true_params$A22^2)
  true_params_with_sigma <- cbind(true_params, sigma = true_sigma)
  lower_bounds_trace <- sapply(true_params_with_sigma[, param_names_trace], min)
  upper_bounds_trace <- sapply(true_params_with_sigma[, param_names_trace], max)
  padding_trace <- (upper_bounds_trace - lower_bounds_trace) * 0.1
  padding_trace[padding_trace == 0] <- 1
  lower_bounds_trace <- lower_bounds_trace - padding_trace
  upper_bounds_trace <- upper_bounds_trace + padding_trace
  
  # Ensure angle phi stays within its valid range [0, pi/2]
  lower_bounds_op["phi"] <- 0
  upper_bounds_op["phi"] <- pi / 2
  lower_bounds_trace["phi"] <- 0
  upper_bounds_trace["phi"] <- pi / 2
  
  # --- Analysis Loop ---
  for (N in N_values) {
    for (m in 1:M_repetitions) {
      
      # Select the correct 2 columns for the current repetition 'm'
      start_col <- (m - 1) * 2 + 1
      end_col <- m * 2
      current_realization <- sim_data[, start_col:end_col]
      
      set.seed(123 + N + m)
      train_indices <- sample(1:N_tot, N)
      train_coords <- grid_points[train_indices, ]
      train_data <- current_realization[train_indices, ]
      
      # Use the mean of the true parameters as the starting point
      initial_pos_vector_op <- colMeans(true_params[, param_names_op])
      initial_pos_vector_trace <- c(colMeans(true_params[, c("lambda1", "lambda2", "phi")]), mean(true_sigma))
      
      # Loop through each of the four methods
      for (method in c("Trace-NS", "Op-NS", "Trace-Stat", "Op-Stat")) {
        cat(sprintf("  -> Analyzing: %s, N=%d, Rep=%d, Method=%s\n", scenario_name, N, m, method))
        
        try({
          is_stationary <- grepl("-Stat", method)
          is_trace <- grepl("Trace-", method)
          
          # Set parameters based on method type
          model_dim <- if (is_trace) 1 else 2
          current_param_names <- if (is_trace) param_names_trace else param_names_op
          current_initial_pos <- if (is_trace) initial_pos_vector_trace else initial_pos_vector_op
          current_lower_bounds <- if (is_trace) lower_bounds_trace else lower_bounds_op
          current_upper_bounds <- if (is_trace) upper_bounds_trace else upper_bounds_op
          
          if (is_stationary) {
            # --- Stationary Method Logic ---
            single_anchor_point <- matrix(c(0, 0), nrow = 1)
            vario <- variogram.lsm(
              z = train_data, d = train_coords, anchorpoints = single_anchor_point,
              epsilon = 1e6, n_angles = 8, n_intervals = 8, kernel_id = "identity",
              dim = model_dim, print_output = FALSE
            )
            solu <- findsolutions.lsm(
              vario, id = "exponential", initial.position = current_initial_pos, 
              lower.bound = current_lower_bounds, upper.bound = current_upper_bounds, print_output = FALSE, remove_not_convergent = !is_stationary
            )
            stationary_params_vector <- as.numeric(solu$solutions)
            est_params <- t(replicate(N_tot, stationary_params_vector))
            colnames(est_params) <- current_param_names
          } else {
            # --- Non-Stationary Method Logic ---
            anchors <- find_anchorpoints.lsm(train_coords, n = 12, plot_output = FALSE)
            vario <- variogram.lsm(
              z = train_data, d = train_coords, anchorpoints = anchors$anchorpoints,
              epsilon = .8, n_angles = 8, n_intervals = 8, kernel_id = "gaussian",
              dim = model_dim, print_output = FALSE
            )
            solu <- findsolutions.lsm(
              vario, id = "exponential", initial.position = current_initial_pos,
              lower.bound = current_lower_bounds, upper.bound = current_upper_bounds, print_output = FALSE
            )
            smoothed_params_df <- smooth.lsm(solu, grid_points)
            est_params <- smoothed_params_df$parameters
            colnames(est_params) <- current_param_names
          }
          
          # Predict on the entire grid
          pred <- predict.lsm(solu, grid_points, print_output = FALSE, plot_output = FALSE)
          
          # Evaluate MSE
          if (is_trace) {
            est_sigma <- est_params[, "sigma"]
          } else {
            est_sigma <- sqrt(est_params[, "A11"]^2 + est_params[, "A21"]^2 + est_params[, "A22"]^2)
          }
          
          temp_results <- data.frame(
            Scenario = scenario_name, Method = method, N = N, Repetition = m,
            MSE_lambda1 = calculate_mse(true_params$lambda1, est_params[, "lambda1"]),
            MSE_lambda2 = calculate_mse(true_params$lambda2, est_params[, "lambda2"]),
            MSE_phi = calculate_mse(true_params$phi, est_params[, "phi"]),
            MSE_sigma = calculate_mse(true_sigma, est_sigma),
            MSE_mean = calculate_mse(true_params[, c("m1", "m2")], pred$smoothed_means),
            MSE_process = calculate_mse(current_realization, pred$zpredicted)
          )
          
          results_df <- rbind(results_df, temp_results)
          
        }, silent = FALSE)
      }
    }
  }
  return(results_df)
}

# -----------------------------------------------------------------------------
# Main Execution
# -----------------------------------------------------------------------------

# Simulation parameters from the paper
N_values <- c(128, 256, 512, 1024)
M_repetitions <- 20

# --- Simplified run for demonstration ---
# To run the full study, use the full N_values and M_repetitions
demo_N_values <- c(128)
demo_M_repetitions <- 20

all_results <- data.frame()

# Run all scenarios
scenarios <- list(
  "Non-Stationary Proportional" = params_nsprop,
  "Non-Stationary Non-Proportional" = params_nsnonprop,
  "Stationary" = params_stat
)

for (scen_name in names(scenarios)) {
  res <- run_simulation(scen_name, scenarios[[scen_name]], demo_N_values, demo_M_repetitions)
  all_results <- rbind(all_results, res)
}

# -----------------------------------------------------------------------------
# Results Aggregation and Plotting
# -----------------------------------------------------------------------------

# Function to create a standardized MSE plot for a given variable
create_mse_plot <- function(summary_data, var_name, N_values) {
  mean_col <- paste0("MSE_", var_name, "_mean")
  sd_col <- paste0("MSE_", var_name, "_sd")
  plot_title <- paste("Simulation Study Results:", tools::toTitleCase(var_name), "MSE")
  
  p <- ggplot(summary_data, aes_string(x = "N", y = mean_col, color = "Method", group = "Method")) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(
      aes_string(ymin = paste0(mean_col, " - ", sd_col), ymax = paste0(mean_col, " + ", sd_col)),
      width = 0.1 * min(diff(sort(unique(summary_data$N))), 1, na.rm = TRUE)
    ) +
    scale_y_log10() +
    scale_x_continuous(breaks = N_values) +
    facet_wrap(~Scenario, scales = "free_y") +
    labs(
      title = plot_title,
      x = "Number of Training Points (N)",
      y = "Mean Squared Error (log scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  
  return(p)
}


if (nrow(all_results) > 0) {
  summary_results <- all_results %>%
    group_by(Scenario, Method, N) %>%
    summarise(across(starts_with("MSE_"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), .groups = 'drop')
  
  print(summary_results)
  
  # List of variables to plot
  plot_vars <- c("lambda1", "lambda2", "phi", "sigma", "mean", "process")
  
  # Generate and print a plot for each variable
  for (var in plot_vars) {
    plot <- create_mse_plot(summary_results, var, N_values)
    print(plot)
  }
  
} else {
  cat("No results to plot. The simulation might have been skipped or encountered errors.\n")
}
