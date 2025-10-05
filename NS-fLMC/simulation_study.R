# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SCRIPT: Simulation for Trace-NS vs Op-NS Models
#
# DESCRIPTION:
# This script conducts a comprehensive oracle simulation study to compare the
# predictive performance of two multivariate non-stationary kriging models:
# 1. Op-NS Cokriging: Utilizes the full operator-based cross-covariance model.
# 2. Trace-NS Kriging: Employs a simplified model based on the trace-covariogram.
#
# The simulation investigates six coregionalization scenarios to assess performance
# under both correctly specified and misspecified model assumptions.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- SECTION 1: SCRIPT SETUP ---

# Install missing packages if necessary
packages <- c("dplyr", "ggplot2", "tidyr", "Matrix",
              "future.apply", "progressr", "gridExtra", "ellipse", "ggrepel")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(LocallyStationaryModels)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix)
library(future.apply)
library(progressr)
library(gridExtra) # For arranging plots
library(grid)      # For textGrob and gpar
library(ellipse)   # For anisotropy ellipses
library(ggrepel)   # For non-overlapping labels

# --- SECTION 2: CORE PARAMETER AND COVARIANCE FUNCTIONS ---

#' Generate All Parameters for a Given Scenario
#'
#' This function is the single source of truth for all model parameters. It defines
#' the spatially varying anisotropy and the six coregionalization scenarios.
#'
#' @param coords A matrix of spatial coordinates (x, y).
#' @param scenario A string specifying the coregionalization scenario.
#' @return A list containing parameters formatted for simulation and analysis.
generate_scenario_parameters <- function(coords,
                                         scenario = c("non-proportional", "prop_center",
                                                      "prop_v1", "prop_v2", "prop_v3", "prop_v4")) {
  scenario <- match.arg(scenario)
  
  x <- coords[, 1]; y <- coords[, 2]
  n_pts <- nrow(coords)
  
  # --- Common Spatially Varying Parameters ---
  lambda1 <- (18 + 6 * x) / 100
  lambda2 <- (8 + 3 * x) / 10
  phi <- (pi / 6) * (2 + 2 * x)
  m1 <- 0 * x; m2 <- 0 * y # Zero mean
  
  # Spatially varying scalar intensity r(s)
  r_s <- 0.3 * (sqrt(x^2 + y^2) + 1)
  
  # --- Scenario-Specific Structural Matrix K(s) ---
  if (scenario == "non-proportional") {
    # Fully non-proportional case: K(s) varies at every point
    K11 <- 2 - x;  K12 <- 0
    K21 <- -2;     K22 <- 2 + y
  } else {
    # Proportional cases: K is a fixed matrix, derived from a specific point
    fixed_coords <- switch(scenario,
                           "prop_center" = c(0, 0),    # Proportional @ (0, 0)
                           "prop_v1"     = c(-1, -1),  # Proportional @ (-1, -1)
                           "prop_v2"     = c(1, -1),   # Proportional @ (1, -1)
                           "prop_v3"     = c(-1, 1),   # Proportional @ (-1, 1)
                           "prop_v4"     = c(1, 1)     # Proportional @ (1, 1)
    )
    x_fixed <- fixed_coords[1]; y_fixed <- fixed_coords[2]
    
    # Calculate the fixed structural matrix K
    K_fixed <- matrix(c(2 - x_fixed, -2, 0, 2 + y_fixed), nrow = 2)
    
    # Replicate this fixed matrix for all locations
    K11 <- rep(K_fixed[1, 1], n_pts); K12 <- rep(K_fixed[1, 2], n_pts)
    K21 <- rep(K_fixed[2, 1], n_pts); K22 <- rep(K_fixed[2, 2], n_pts)
  }
  
  # Final Coregionalization Matrix A(s) = r(s) * K(s)
  A11 <- r_s * K11; A12 <- r_s * K12
  A21 <- r_s * K21; A22 <- r_s * K22
  
  # --- Assemble Output Formats ---
  params_for_sampling <- as.matrix(data.frame(
    lambda1 = lambda1, lambda2 = lambda2, phi = phi,
    A11 = A11, A21 = A21, A22 = A22
  ))
  
  params_list <- lapply(1:n_pts, function(i) {
    rotation_matrix <- matrix(c(cos(phi[i]), sin(phi[i]), -sin(phi[i]), cos(phi[i])), 2, 2)
    aniso_matrix <- rotation_matrix %*% diag(c(lambda1[i], lambda2[i])^2) %*% t(rotation_matrix)
    list(
      Sigma_aniso = aniso_matrix,
      A = matrix(c(A11[i], A21[i], A12[i], A22[i]), 2, 2),
      m = c(m1[i], m2[i])
    )
  })
  
  true_means <- as.matrix(data.frame(m1 = m1, m2 = m2))
  
  return(list(
    params_list = params_list,
    params_for_sampling = params_for_sampling,
    true_means = true_means
  ))
}

#' Compute the Non-Stationary Correlation R_NS(s_i, s_j)
#' @return A scalar correlation value.
compute_R_NS <- function(si, sj, params_i, params_j) {
  d <- 2
  Sigma_i <- params_i$Sigma_aniso; Sigma_j <- params_j$Sigma_aniso
  
  # Ensure matrices are positive definite to avoid numerical issues
  if(det(Sigma_i) <= 1e-12) Sigma_i <- as.matrix(nearPD(Sigma_i, base.matrix = TRUE)$mat)
  if(det(Sigma_j) <= 1e-12) Sigma_j <- as.matrix(nearPD(Sigma_j, base.matrix = TRUE)$mat)
  
  Sigma_sum <- Sigma_i + Sigma_j
  if(det(Sigma_sum) <= 1e-12) return(exp(-sqrt(sum((si - sj)^2)))) # Fallback
  
  term1 <- (2^(d / 2) * det(Sigma_i)^(1 / 4) * det(Sigma_j)^(1 / 4)) / det(Sigma_sum)^(1 / 2)
  Q_ij <- t(si - sj) %*% solve(Sigma_sum / 2) %*% (si - sj)
  R_S_val <- exp(-sqrt(Q_ij))
  
  return(as.numeric(term1 * R_S_val))
}

# --- SECTION 3: SIMULATION ENGINE ---

#' Run a Single Oracle Simulation Setting
#'
#' Executes M_repetitions for a given scenario and training size.
#'
#' @param scenario The coregionalization scenario name.
#' @param M_repetitions Number of Monte Carlo repetitions.
#' @param N_train Number of training points.
#' @param p Number of variables (dimensionality of the process).
#' @param nugget The nugget variance.
#' @return A data frame with performance metrics for each repetition.
run_oracle_simulation <- function(scenario, M_repetitions = 50, N_train = 100, p = 2, nugget = 1e-4) {
  
  # 1. Generate Parameters and Underlying Spatial Process
  grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 50)))
  all_params <- generate_scenario_parameters(grid_points, scenario)
  
  # Generate the smooth part of the spatial process (M_repetitions)
  sim_result_smooth <- LocallyStationaryModels:::samplelsm(
    d = grid_points, variogram_id = "exponential",
    parameters = all_params$params_for_sampling, dim = p, n_samples = M_repetitions
  )$simulated_processes
  
  # Add nugget effect to create the final observed data for all repetitions
  nugget_noise <- matrix(rnorm(nrow(grid_points) * p * M_repetitions, mean = 0, sd = sqrt(nugget)),
                         nrow = nrow(grid_points), ncol = p * M_repetitions)
  all_sim_residuals <- sim_result_smooth + nugget_noise
  
  # 2. Pre-computation (Invariant across repetitions)
  N_tot <- nrow(grid_points)
  train_indices <- sample(1:N_tot, N_train)
  test_indices <- setdiff(1:N_tot, train_indices)
  N_test <- length(test_indices)
  
  train_coords <- grid_points[train_indices, ]; test_coords <- grid_points[test_indices, ]
  train_params <- all_params$params_list[train_indices]; test_params <- all_params$params_list[test_indices]
  
  # Pre-compute covariance matrices and kriging weights
  # Op-NS Model Covariances
  C_block_train_Op <- matrix(0, nrow = N_train * p, ncol = N_train * p)
  for (i in 1:N_train) {
    for (j in i:N_train) { # Exploit symmetry
      R_ns_ij <- compute_R_NS(train_coords[i,], train_coords[j,], train_params[[i]], train_params[[j]])
      C_ij_op <- train_params[[i]]$A %*% t(train_params[[j]]$A) * R_ns_ij
      idx_i <- ((i - 1) * p + 1):(i * p); idx_j <- ((j - 1) * p + 1):(j * p)
      C_block_train_Op[idx_i, idx_j] <- C_ij_op
      if (i != j) C_block_train_Op[idx_j, idx_i] <- t(C_ij_op)
    }
  }
  inv_C_block_train_Op <- chol2inv(chol(C_block_train_Op + diag(nugget, N_train * p)))
  
  # Trace-NS Model Covariances
  C_trace_train <- matrix(0, nrow = N_train * p, ncol = N_train * p)
  for (i in 1:N_train) {
    for (j in i:N_train) {
      R_ns_ij <- compute_R_NS(train_coords[i,], train_coords[j,], train_params[[i]], train_params[[j]])
      C_trace_ij <- sum(diag(train_params[[i]]$A %*% t(train_params[[j]]$A))) * R_ns_ij
      for(k in 1:p){
        C_trace_train[(i-1)*p+k, (j-1)*p+k] <- C_trace_ij
        if(i!=j) C_trace_train[(j-1)*p+k, (i-1)*p+k] <- C_trace_ij
      }
    }
  }
  inv_C_trace_train <- chol2inv(chol(C_trace_train + diag(nugget, N_train * p)))
  
  
  # Cross-covariance matrices (Train <-> Test)
  C_block_traintest_Op <- matrix(0, nrow = N_train * p, ncol = N_test * p)
  C_trace_traintest <- matrix(0, nrow = N_train * p, ncol = N_test * p)
  
  for (i in 1:N_train) {
    for (j in 1:N_test) {
      R_ns_ij <- compute_R_NS(train_coords[i,], test_coords[j,], train_params[[i]], test_params[[j]])
      idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
      
      C_block_traintest_Op[idx_i, idx_j] <- train_params[[i]]$A %*% t(test_params[[j]]$A) * R_ns_ij
      C_trace_ij_test <- sum(diag(train_params[[i]]$A %*% t(test_params[[j]]$A))) * R_ns_ij
      for(k in 1:p){
        C_trace_traintest[(i-1)*p+k, (j-1)*p+k] <- C_trace_ij_test
      }
    }
  }
  
  weights_Op <- inv_C_block_train_Op %*% C_block_traintest_Op
  weights_Trace <- inv_C_trace_train %*% C_trace_traintest
  
  
  # 3. Main Repetition Loop (now faster)
  results_list <- lapply(1:M_repetitions, function(m) {
    start_col <- (m - 1) * p + 1; end_col <- m * p
    sim_data <- all_sim_residuals[, start_col:end_col] + all_params$true_means
    train_data <- sim_data[train_indices, ]; test_data <- sim_data[test_indices, ]
    residuals_train <- as.vector(t(train_data - all_params$true_means[train_indices,]))
    
    # Op-NS Prediction
    pred_op_ns <- matrix(t(weights_Op) %*% residuals_train, ncol = p, byrow = TRUE) + all_params$true_means[test_indices,]
    mspe_op_ns <- mean(rowSums((test_data - pred_op_ns)^2))
    
    # Trace-NS Prediction
    pred_trace_ns <- matrix(t(weights_Trace) %*% residuals_train, ncol = p, byrow = TRUE) + all_params$true_means[test_indices,]
    mspe_trace_ns <- mean(rowSums((test_data - pred_trace_ns)^2))
    
    data.frame(Repetition = m, MSPE_Op_NS = mspe_op_ns, MSPE_Trace_NS = mspe_trace_ns)
  })
  
  return(do.call(rbind, results_list))
}


# --- SECTION 4: VISUALIZATION OF SIMULATION SETUP ---

#' Create and save plots to visualize the simulation design.
#'
#' This function generates two separate plots: one for the K(s) matrix components
#' and one for the scalar intensity field r(s).
#'
create_and_save_setup_plots <- function() {
  # Create a dense grid for smooth plotting
  grid_dense <- as.matrix(expand.grid(x = seq(-1, 1, length.out=100),
                                      y = seq(-1, 1, length.out=100)))
  x <- grid_dense[, 1]; y <- grid_dense[, 2]
  
  # --- Data for K(s) components (from non-proportional scenario) ---
  K_components_data <- data.frame(
    x = x, y = y,
    K11 = 2 - x,
    K21 = -2,
    K22 = 2 + y
  ) %>%
    pivot_longer(cols=c(K11, K21, K22), names_to="component", values_to="value") %>%
    mutate(component = factor(component,
                              levels = c("K11", "K21", "K22"),
                              labels = c("K[11](s)", "K[21](s)", "K[22](s)")))
  
  # --- Data for proportional scenario points ---
  proportional_points <- data.frame(
    x = c(0, -1, 1, -1, 1),
    y = c(0, -1, -1, 1, 1),
    label = c("(0,0)", "(-1,-1)", "(1,-1)", "(-1,1)", "(1,1)")
  )
  
  # --- Plot 1: K(s) components with points ---
  p_K <- ggplot(K_components_data, aes(x, y, fill=value)) +
    geom_raster() +
    facet_wrap(~component, nrow=1, labeller = label_parsed) +
    geom_point(data = proportional_points, aes(x, y), inherit.aes = FALSE,
               color = "white", size = 4, shape = 18) +
    geom_text_repel(data = proportional_points, aes(x, y, label = label), inherit.aes = FALSE,
                    color = "white", size = 3, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.5,
                    segment.color = 'white', segment.size = 0.5) +
    scale_fill_viridis_c() +
    coord_fixed() +
    labs(title="Components of Structural Matrix K(s) in Non-Proportional Scenario",
         subtitle="White diamonds mark the locations used to define the fixed K for the Proportional scenarios.",
         x="Coordinate x", y="Coordinate y", fill="Value") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          strip.text = element_text(size=12, face="bold"),
          axis.title = element_text(face="bold"))
  
  # --- Data for r(s) scalar field ---
  r_s_data <- data.frame(
    x = x, y = y,
    r_s = 0.3 * (sqrt(x^2 + y^2) + 1)
  )
  
  # --- Plot 2: r(s) field ---
  p_r <- ggplot(r_s_data, aes(x, y, fill=r_s)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma") +
    coord_fixed() +
    labs(title="Scalar Intensity Field r(s)",
         subtitle="This field is common to all six simulation scenarios.",
         x="Coordinate x", y="Coordinate y", fill="Value") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=10),
          axis.title = element_text(face="bold"))
  
  # --- Save the plots ---
  ggsave("K_structure_plot.png", plot = p_K, width = 11, height = 4.5, dpi = 300)
  ggsave("r_intensity_plot.png", plot = p_r, width = 7, height = 6, dpi = 300)
  
  cat("--- Simulation setup plots saved to K_structure_plot.png and r_intensity_plot.png ---\n")
}

# Generate and save the setup plots before running the main simulation
create_and_save_setup_plots()


# --- SECTION 5: MAIN EXECUTION BLOCK ---
set.seed(0) # for reproducibility
M_rep <- 350       # Number of repetitions for stable estimates
N_values <- c(20, 50, 100, 120, 150, 200, 500) # Training sizes

scenarios_to_run <- c(
  "non-proportional",
  "prop_v1",          # Proportional based on vertex (-1, -1)
  "prop_v2",          # Proportional based on vertex (1, -1)
  "prop_center",
  "prop_v3",          # Proportional based on vertex (-1, 1)
  "prop_v4"           # Proportional based on vertex (1, 1)
)

# Set up parallel processing
plan(multisession, workers = availableCores() - 4)

# Create the full experimental grid
experiment_grid <- expand.grid(
  scenario = scenarios_to_run,
  n_train = N_values,
  stringsAsFactors = FALSE
)

# Enable progress bars for future_lapply
handlers(global = TRUE)
handlers("progress")

# Run the simulation experiment
cat("--- STARTING ORACLE SIMULATION ---\n")
with_progress({
  p_progress <- progressor(steps = nrow(experiment_grid))
  
  all_results_list <- future_lapply(1:nrow(experiment_grid), function(i) {
    setting <- experiment_grid[i, ]
    p_progress(sprintf("Scenario=%s, N_train=%g", setting$scenario, setting$n_train))
    
    results_for_setting <- run_oracle_simulation(
      scenario = setting$scenario,
      M_repetitions = M_rep,
      N_train = setting$n_train
    )
    
    results_for_setting$Scenario <- setting$scenario
    results_for_setting$N_train <- setting$n_train
    return(results_for_setting)
    
  }, future.seed = TRUE)
})
cat("\n--- ALL SIMULATIONS COMPLETE ---\n")

simulation_results <- do.call(rbind, all_results_list)


# --- SECTION 6: ANALYSIS AND PLOTTING RESULTS ---

# Create descriptive labels for facets based on model specification
simulation_results_labeled <- simulation_results %>%
  mutate(Scenario_Label = case_when(
    Scenario == "non-proportional" ~ "Non-Proportional\n",
    Scenario == "prop_center"      ~ "Proportional @ (0, 0)\n ",
    Scenario == "prop_v1"          ~ "Proportional @ (-1, -1)\n ",
    Scenario == "prop_v2"          ~ "Proportional @ (1, -1)\n ",
    Scenario == "prop_v3"          ~ "Proportional @ (-1, 1)\n ",
    Scenario == "prop_v4"          ~ "Proportional @ (1, 1)\n ",
    TRUE ~ Scenario
  )) %>%
  # Set the order for the facets
  mutate(Scenario_Label = factor(Scenario_Label, levels = c(
    "Non-Proportional\n",
    "Proportional @ (-1, -1)\n ",
    "Proportional @ (1, -1)\n ",
    "Proportional @ (0, 0)\n ",
    "Proportional @ (1, 1)\n ",
    "Proportional @ (-1, 1)\n "
  )))


# Calculate the paired difference for MSPE and the p-value from a paired t-test
summary_stats <- simulation_results_labeled %>%
  group_by(Scenario_Label, N_train) %>%
  summarise(
    Mean_Diff_MSPE = mean(MSPE_Trace_NS - MSPE_Op_NS, na.rm = TRUE),
    SE_Diff_MSPE = sd(MSPE_Trace_NS - MSPE_Op_NS, na.rm = TRUE) / sqrt(n()),
    P_Value = tryCatch({ t.test(MSPE_Trace_NS, MSPE_Op_NS, paired = TRUE)$p.value }, error = function(e) { NA_real_ }),
    .groups = 'drop'
  ) %>%
  mutate(
    CI_Lower_MSPE = Mean_Diff_MSPE - 1.96 * SE_Diff_MSPE,
    CI_Upper_MSPE = Mean_Diff_MSPE + 1.96 * SE_Diff_MSPE
  )

# Prepare data for plotting in a faceted grid
plot_data_diff <- summary_stats %>%
  select(Scenario_Label, N_train, Value = Mean_Diff_MSPE, CI_Lower = CI_Lower_MSPE, CI_Upper = CI_Upper_MSPE) %>%
  mutate(Metric = "MSPE Difference")

plot_data_pval <- summary_stats %>%
  select(Scenario_Label, N_train, Value = P_Value) %>%
  mutate(Metric = "Paired t-test p-value", CI_Lower = NA, CI_Upper = NA)

plot_data_final <- bind_rows(plot_data_diff, plot_data_pval) %>%
  mutate(Metric = factor(Metric, levels = c("MSPE Difference", "Paired t-test p-value")))

hline_data <- data.frame(
  Metric = factor(levels(plot_data_final$Metric), levels = levels(plot_data_final$Metric)),
  intercept = c(0, 0.05)
)

significance_plot <- ggplot(plot_data_final, aes(x = N_train, y = Value)) +
  geom_hline(data = hline_data, aes(yintercept = intercept), linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_ribbon(data = . %>% filter(Metric == "MSPE Difference"),
              aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.25, fill = "#0072B2") +
  geom_line(linewidth = 1, color = "#0072B2") +
  geom_point(size = 2.5, color = "#0072B2") +
  facet_grid(Metric ~ Scenario_Label, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(100, 300, 500)) +
  labs(
    title = "Oracle Kriging Performance: MSPE Difference and Statistical Significance",
    subtitle = "Positive MSPE difference favors Op-NS model. Lower panel shows p-values for paired t-test.",
    x = "Number of Training Points (N_train)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    strip.text.x = element_text(face = "bold", size = 9.5),
    strip.text.y = element_text(face = "bold", size = 12),
    strip.placement = "outside",
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Save the plot to a file
ggsave("mspe_significance_plot.png", plot = significance_plot, width = 14, height = 7, dpi = 300)

cat("\n--- PLOT SAVED to mspe_significance_plot.png ---\n")
print(significance_plot)

