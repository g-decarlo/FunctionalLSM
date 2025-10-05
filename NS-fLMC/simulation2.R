# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SCRIPT: Oracle Simulation for Trace-NS vs Op-NS (Resampling Version)
#
# ARCHITECTURE: This version has been modified to resample training points
# at each repetition, moving all matrix computations inside the main loop.
# This increases statistical independence at the cost of performance.
#
# Date: 2025-10-03
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Ensure required libraries are loaded
if (!require(LocallyStationaryModels)) install.packages("LocallyStationaryModels")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(Matrix)) install.packages("Matrix")
# For optional parallel processing
if (!require(future.apply)) install.packages("future.apply")
if (!require(progressr)) install.packages("progressr")


library(LocallyStationaryModels)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix) # For nearPD
library(future.apply) # For parallel processing
library(progressr) # To show progress bars with future


# -----------------------------------------------------------------------------
# SECTION 1: CORE LOGIC & PARAMETER GENERATION (Unchanged)
# -----------------------------------------------------------------------------

#' @title Generate All Parameters for a Given Scenario
#' @description This is the single source of truth for all model parameters.
generate_scenario_parameters <- function(coords, scenario = c("non-proportional", "proportional")) {
  scenario <- match.arg(scenario)
  
  x <- coords[, 1]; y <- coords[, 2]
  n_pts <- nrow(coords)
  
  # --- Common Parameters (using your latest version) ---
  lambda1 <- (18 + 6 * x) / 100 
  lambda2 <- (8 + 3 * x) / 10
  phi <- (pi / 6) * (2 + 2*x)
  m1 <-  0*x*y; m2 <- 0*y*x
  
  r <- sqrt(x^2 + y^2) + 1
  alpha_s <- .3*r
  
  # --- Scenario-Specific Coregionalization Matrix A(s) ---
  if (scenario == "non-proportional") {
    A11_base <-  2 - x; A12_base <- 0
    A21_base <- -2;     A22_base <- 2+y
    A11 <- alpha_s * A11_base; A12 <- alpha_s * A12_base
    A21 <- alpha_s * A21_base; A22 <- alpha_s * A22_base
  } else { # scenario == "proportional"
    A_fixed <- matrix(c(2, -2, 0, 2), nrow = 2)
    A11 <- alpha_s * A_fixed[1, 1]; A12 <- alpha_s * A_fixed[1, 2]
    A21 <- alpha_s * A_fixed[2, 1]; A22 <- alpha_s * A_fixed[2, 2]
  }
  
  # --- Assemble Output Formats ---
  params_for_sampling <- as.matrix(data.frame(
    lambda1 = lambda1,  
    lambda2 = lambda2,  
    phi = phi,
    A11 = A11, A21 = A21, A22 = A22
  ))
  
  params_list <- lapply(1:n_pts, function(i) {
    rotation_matrix <- matrix(c(cos(phi[i]), sin(phi[i]), -sin(phi[i]), cos(phi[i])), 2, 2)
    aniso_matrix <- rotation_matrix %*% diag(c(lambda1[i], lambda2[i])^2) %*% t(rotation_matrix)
    
    list(
      lambda1 = lambda1[i],
      lambda2 = lambda2[i],
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

# -----------------------------------------------------------------------------
# SECTION 2: HELPER FUNCTIONS & SIMULATION ORCHESTRATOR (MODIFIED)
# -----------------------------------------------------------------------------

sample.lsm <- function(d, variogram_id, parameters, dim = 1, n_samples = 1) {
  result <- LocallyStationaryModels:::samplelsm(d, variogram_id, parameters, dim, n_samples)
  return(result)
}

compute_R_NS <- function(si, sj, params_i, params_j) {
  d <- 2
  Sigma_i <- params_i$Sigma_aniso; Sigma_j <- params_j$Sigma_aniso
  
  if(det(Sigma_i) <= 1e-12) Sigma_i <- as.matrix(nearPD(Sigma_i, base.matrix = TRUE)$mat)
  if(det(Sigma_j) <= 1e-12) Sigma_j <- as.matrix(nearPD(Sigma_j, base.matrix = TRUE)$mat)
  
  Sigma_sum <- Sigma_i + Sigma_j
  if(det(Sigma_sum) <= 1e-12) return(exp(-sqrt(sum((si-sj)^2))))
  
  term1 <- (2^(d/2) * det(Sigma_i)^(1/4) * det(Sigma_j)^(1/4)) / det(Sigma_sum)^(1/2)
  Q_ij <- t(si - sj) %*% solve(Sigma_sum / 2) %*% (si - sj)
  R_S_val <- exp(-sqrt(Q_ij))
  
  return(as.numeric(term1 * R_S_val))
}


run_oracle_simulation <- function(scenario, M_repetitions = 50, N_train = 100, p = 2, nugget = 1e-4) {
  
  cat(paste("\n--- Running:", scenario, "scenario with N_train =", N_train, "---\n"))
  
  # 1. Generate ALL potential data ONCE
  grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 50)))
  all_params <- generate_scenario_parameters(grid_points, scenario)
  
  # Generate the smooth part of all spatial process realizations
  sim_result_smooth <- sample.lsm(d = grid_points, variogram_id = "exponential",
                                  parameters = all_params$params_for_sampling, dim = p, n_samples = M_repetitions)
  
  # Add the nugget effect to all realizations
  nugget_noise <- matrix(rnorm(nrow(grid_points) * p * M_repetitions, mean = 0, sd = sqrt(nugget)),
                         nrow = nrow(grid_points), ncol = p * M_repetitions)
  
  all_sim_residuals <- sim_result_smooth$simulated_processes + nugget_noise
  N_tot <- nrow(grid_points)
  
  # 2. Main Loop (Now includes all computations)
  results <- data.frame(Repetition = 1:M_repetitions, 
                        MSPE_Op_NS = NA, MSPE_Trace_NS = NA,
                        R2_Op_NS = NA, R2_Trace_NS = NA)
  
  for (m in 1:M_repetitions) {
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # START OF MOVED BLOCK: All computations are now inside the repetition loop
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 2A. Sample new training/test sets for THIS repetition
    train_indices <- sample(1:N_tot, N_train)
    test_indices <- setdiff(1:N_tot, train_indices)
    N_test <- length(test_indices)
    
    train_coords <- grid_points[train_indices, ]; test_coords <- grid_points[test_indices, ]
    train_params <- all_params$params_list[train_indices]; test_params <- all_params$params_list[test_indices]
    
    # 2B. Build covariance matrices for this specific training set
    C_block_train_Op <- matrix(0, nrow = N_train * p, ncol = N_train * p)
    C_trace_train_Naive <- matrix(0, nrow = N_train, ncol = N_train)
    
    for (i in 1:N_train) {
      sigma_i <- sqrt(sum(train_params[[i]]$A^2))
      for (j in i:N_train) { # Only compute upper triangle
        sigma_j <- sqrt(sum(train_params[[j]]$A^2))
        R_ns_ij <- compute_R_NS(train_coords[i,], train_coords[j,], train_params[[i]], train_params[[j]])
        
        C_ij_op <- train_params[[i]]$A %*% t(train_params[[j]]$A) * R_ns_ij
        idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
        C_block_train_Op[idx_i, idx_j] <- C_ij_op
        if (i != j) C_block_train_Op[idx_j, idx_i] <- t(C_ij_op)
        
        C_trace_train_Naive[i, j] <- sigma_i * sigma_j * R_ns_ij
        if (i != j) C_trace_train_Naive[j, i] <- C_trace_train_Naive[i, j]
      }
    }
    
    inv_C_block_train_Op <- chol2inv(chol(C_block_train_Op + diag(nugget, N_train * p)))
    inv_C_trace_train_Naive <- chol2inv(chol(C_trace_train_Naive + diag(p * nugget, N_train)))
    
    # 2C. Build cross-covariance matrices and calculate weights
    C_block_traintest_Op <- matrix(0, nrow = N_train * p, ncol = N_test * p)
    C_trace_traintest_Naive <- matrix(0, nrow = N_train, ncol = N_test)
    
    for (i in 1:N_train) {
      sigma_i <- sqrt(sum(train_params[[i]]$A^2))
      for (j in 1:N_test) {
        sigma_j <- sqrt(sum(test_params[[j]]$A^2))
        R_ns_ij <- compute_R_NS(train_coords[i,], test_coords[j,], train_params[[i]], test_params[[j]])
        
        idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
        C_block_traintest_Op[idx_i, idx_j] <- train_params[[i]]$A %*% t(test_params[[j]]$A) * R_ns_ij
        C_trace_traintest_Naive[i, j] <- sigma_i * sigma_j * R_ns_ij
      }
    }
    
    weights_Op <- inv_C_block_train_Op %*% C_block_traintest_Op
    weights_Naive <- inv_C_trace_train_Naive %*% C_trace_traintest_Naive
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # END OF MOVED BLOCK
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # 2D. Prediction and Evaluation for repetition 'm'
    start_col <- (m-1)*p + 1; end_col <- m*p
    sim_data <- all_sim_residuals[, start_col:end_col] + all_params$true_means
    train_data <- sim_data[train_indices, ]; test_data <- sim_data[test_indices, ]
    residuals_train <- train_data - all_params$true_means[train_indices,]
    
    # Op-NS Prediction
    pred_op_ns <- matrix(t(weights_Op) %*% as.vector(t(residuals_train)), ncol = p, byrow = TRUE) + all_params$true_means[test_indices,]
    results$MSPE_Op_NS[m] <- mean(rowSums((test_data - pred_op_ns)^2))
    
    # Naive Trace Prediction
    pred_naive <- t(weights_Naive) %*% residuals_train + all_params$true_means[test_indices,]
    results$MSPE_Trace_NS[m] <- mean(rowSums((test_data - pred_naive)^2))
    
    # R-squared calculations
    TSS <- sum(sweep(test_data, 2, colMeans(test_data), "-")^2)
    SSR_op_ns <- sum((test_data - pred_op_ns)^2)
    SSR_trace_ns <- sum((test_data - pred_naive)^2)
    
    if (TSS > 0) {
      results$R2_Op_NS[m] <- 1 - (SSR_op_ns / TSS)
      results$R2_Trace_NS[m] <- 1 - (SSR_trace_ns / TSS)
    } else {
      results$R2_Op_NS[m] <- NA
      results$R2_Trace_NS[m] <- NA
    }
  } # End of main repetition loop
  
  return(results)
}

# -----------------------------------------------------------------------------
# SECTION 4: MAIN EXECUTION BLOCK
# -----------------------------------------------------------------------------
set.seed(123)
N_values <- c(10, 20, 50, 100) 
# WARNING: Simulation is now much slower. Consider reducing M_rep for initial tests.
M_rep <- 1 
scenarios_to_run <- c("non-proportional", "proportional")

plan(multisession, workers = availableCores() - 1)

experiment_grid <- expand.grid(
  scenario = scenarios_to_run,
  n_train = N_values,
  stringsAsFactors = FALSE
)

handlers(global = TRUE)
handlers("progress")

with_progress({
  p <- progressor(steps = nrow(experiment_grid))
  
  all_results_list <- future_lapply(1:nrow(experiment_grid), function(i) {
    setting <- experiment_grid[i, ]
    p(sprintf("scenario=%s, n_train=%g", setting$scenario, setting$n_train))
    
    results_for_n <- run_oracle_simulation(
      scenario = setting$scenario,
      M_repetitions = M_rep,
      N_train = setting$n_train
    )
    
    results_for_n$Scenario <- setting$scenario
    results_for_n$N_train <- setting$n_train
    return(results_for_n)
    
  }, future.seed = TRUE)
})

simulation_results <- do.call(rbind, all_results_list)

cat("\n\n--- ALL SIMULATIONS COMPLETE ---\n")

# -----------------------------------------------------------------------------
# SECTION 5: AGGREGATION AND PLOTTING (Unchanged)
# -----------------------------------------------------------------------------

# --- PLOT 1: Direct Performance Comparison ---

# Reshape data into a long format suitable for ggplot faceting
summary_data_long <- simulation_results %>%
  pivot_longer(
    cols = c(MSPE_Op_NS, MSPE_Trace_NS, R2_Op_NS, R2_Trace_NS),
    names_to = "Metric_Method",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = ifelse(grepl("MSPE", Metric_Method), "MSPE", "R-squared"),
    Method = ifelse(grepl("Op_NS", Metric_Method), "Op-NS (Cokriging)", "Trace-NS (Simple Kriging)")
  ) %>%
  select(-Metric_Method)

# Calculate mean, SE, and CI for each group
summary_data <- summary_data_long %>%
  group_by(Scenario, N_train, Method, Metric) %>%
  summarise(
    Mean_Value = mean(Value, na.rm = TRUE),
    SE_Value = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(
    CI_Lower = Mean_Value - 1.96 * SE_Value,
    CI_Upper = Mean_Value + 1.96 * SE_Value,
    Scenario = ifelse(Scenario == "proportional", "Proportional (Trace-NS is Correct)", "Non-Proportional (Trace-NS is Misspecified)")
  )

# Generate the combined 2x2 plot
combined_plot <- ggplot(summary_data, aes(x = N_train, y = Mean_Value, color = Method, group = Method)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Method), alpha = 0.2, linetype = "dotted") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_grid(Metric ~ Scenario, scales = "free_y", switch = "y") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Oracle Kriging Performance vs. Training Sample Size",
    subtitle = "Comparison of MSPE and R-squared under different scenarios",
    x = "Number of Training Points (N_train)",
    y = NULL, 
    color = "Kriging Method",
    fill = "Kriging Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y = element_text(face = "bold", size = 12),
    strip.placement = "outside",
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(face="bold")
  )

# Print the final combined plot
cat("\n--- COMBINED PERFORMANCE PLOT (MSPE & R-SQUARED) ---\n")
print(combined_plot)


# --- PLOT 2: Visualizing the Paired Difference ---

# Calculate the paired difference for both metrics
summary_data_diff <- simulation_results %>%
  mutate(
    Difference_MSPE = MSPE_Trace_NS - MSPE_Op_NS,
    Difference_R2 = R2_Op_NS - R2_Trace_NS # Note the order is reversed for R2
  ) %>%
  pivot_longer(
    cols = c(Difference_MSPE, Difference_R2),
    names_to = "Metric",
    values_to = "Difference"
  ) %>%
  mutate(
    Metric = ifelse(grepl("MSPE", Metric), "MSPE Difference (Trace - Op)", "R-squared Difference (Op - Trace)")
  ) %>%
  group_by(Scenario, N_train, Metric) %>%
  summarise(
    Mean_Diff = mean(Difference, na.rm = TRUE),
    SE_Diff = sd(Difference, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  mutate(
    CI_Lower = Mean_Diff - 1.96 * SE_Diff,
    CI_Upper = Mean_Diff + 1.96 * SE_Diff,
    Scenario = ifelse(Scenario == "proportional", "Proportional (Trace-NS is Correct)", "Non-Proportional (Trace-NS is Misspecified)")
  )

# Generate the difference plot
difference_plot <- ggplot(summary_data_diff, aes(x = N_train, y = Mean_Diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, fill = "darkblue") +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  facet_grid(Metric ~ Scenario, scales = "free_y", switch = "y") +
  labs(
    title = "Paired Difference in Performance",
    subtitle = "Visualizing the statistical significance of the difference between models",
    x = "Number of Training Points (N_train)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    strip.placement = "outside",
    plot.title = element_text(face = "bold", size = 18),
    axis.title = element_text(face="bold")
  )

cat("\n--- PAIRED DIFFERENCE PLOT (MSPE & R-SQUARED) ---\n")
print(difference_plot)