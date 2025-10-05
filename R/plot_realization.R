# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# SECTION 7: VISUALIZING TRUE VS. PREDICTED FIELDS FOR p=3
#
# This script performs a single kriging prediction for a trivariate (p=3)
# process to visually compare the true spatial field against the predictions
# from all three models: Op-NS, Correct Trace-NS, and Misspecified Trace-NS.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# --- 0. Load Required Libraries ---
cat("--- Loading required libraries ---\n")
if (!require(LocallyStationaryModels)) install.packages("LocallyStationaryModels")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(Matrix)) install.packages("Matrix")
if (!require(viridis)) install.packages("viridis")

library(LocallyStationaryModels)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix)
library(viridis)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SECTION 1: CORE FUNCTIONS (ADAPTED FROM p=3 SIMULATION SCRIPT)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Generate All Parameters for a Given Scenario (p=3 version)
generate_scenario_parameters <- function(coords, scenario = c("non-proportional", "proportional")) {
  scenario <- match.arg(scenario)
  
  x <- coords[, 1]; y <- coords[, 2]
  n_pts <- nrow(coords)
  
  # --- Common Parameters (using your latest version) ---
  lambda1 <- (18 + 6 * x) / 100 
  lambda2 <- (8 + 3 * x) / 10
  phi <- (pi / 4) * (2 + 2*x)
  m1 <-  0*x*y; m2 <- 0*y*x; m3 <- 0*x*y # NEW: m3 added
  
  r <- sqrt(x^2 + y^2) + 2
  alpha_s <- .2*r
  
  # --- Scenario-Specific Coregionalization Matrix A(s) (NOW 3x3) ---
  if (scenario == "non-proportional") {
    # Define 9 spatially-varying components for the 3x3 matrix A(s)
    A11_base <-  1;      A12_base <- 0;         A13_base <- 0
    A21_base <- 1;  A22_base <- 1;     A23_base <- 0 
    A31_base <- 0;      A32_base <- 2;        A33_base <- 1
    
    A11 <- alpha_s * A11_base; A12 <- alpha_s * A12_base; A13 <- alpha_s * A13_base
    A21 <- alpha_s * A21_base; A22 <- alpha_s * A22_base; A23 <- alpha_s * A23_base
    A31 <- alpha_s * A31_base; A32 <- alpha_s * A32_base; A33 <- alpha_s * A33_base
    
  } else { # scenario == "proportional"
    # Define a fixed 3x3 matrix for A_0
    A_fixed <- matrix(c(2, 0, 0, -2, 2, 0, -2, -1, 3), nrow = 3)
    
    A11 <- alpha_s * A_fixed[1, 1]; A12 <- alpha_s * A_fixed[1, 2]; A13 <- alpha_s * A_fixed[1, 3]
    A21 <- alpha_s * A_fixed[2, 1]; A22 <- alpha_s * A_fixed[2, 2]; A23 <- alpha_s * A_fixed[2, 3]
    A31 <- alpha_s * A_fixed[3, 1]; A32 <- alpha_s * A_fixed[3, 2]; A33 <- alpha_s * A_fixed[3, 3]
  }
  
  # --- Assemble Output Formats ---
  params_for_sampling <- as.matrix(data.frame(
    lambda1 = lambda1,  lambda2 = lambda2,  phi = phi,
    A11 = A11, A21 = A21, A31 = A31, A22 = A22, A32 = A32, A33 = A33
  ))
  
  params_list <- lapply(1:n_pts, function(i) {
    rotation_matrix <- matrix(c(cos(phi[i]), sin(phi[i]), -sin(phi[i]), cos(phi[i])), 2, 2)
    aniso_matrix <- rotation_matrix %*% diag(c(lambda1[i], lambda2[i])^2) %*% t(rotation_matrix)
    
    list(
      lambda1 = lambda1[i], lambda2 = lambda2[i], Sigma_aniso = aniso_matrix,
      A = matrix(c(A11[i], A21[i], A31[i], A12[i], A22[i], A32[i], A13[i], A23[i], A33[i]), 3, 3),
      m = c(m1[i], m2[i], m3[i])
    )
  })
  
  true_means <- as.matrix(data.frame(m1 = m1, m2 = m2, m3 = m3))
  
  return(list(
    params_list = params_list,
    params_for_sampling = params_for_sampling,
    true_means = true_means
  ))
}

#' @title Non-stationary Correlation Function
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

#' @title LSM Sampling Function Wrapper
sample.lsm <- function(d, variogram_id, parameters, dim = 1, n_samples = 1) {
  result <- LocallyStationaryModels:::samplelsm(d, variogram_id, parameters, dim, n_samples)
  return(result)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SECTION 2: VISUALIZATION EXECUTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("\n--- GENERATING PLOT FOR TRUE VS. PREDICTED FIELDS (p=3) ---\n")

# --- 1. Set parameters for the visualization ---
set.seed(42) 
p <- 3 # CRITICAL: Set dimension to 3
nugget_vis <- 1e-4
grid_resolution <- 50  
N_train_vis <- 150 # Increased training points for better visual results with p=3

# --- 2. Generate the "True" Data ---
cat("...generating true spatial field (p=3)...\n")
plot_grid_points <- as.matrix(expand.grid(
  seq(-1, 1, length.out = grid_resolution),
  seq(-1, 1, length.out = grid_resolution)
))
N_tot_vis <- nrow(plot_grid_points)

plot_params <- generate_scenario_parameters(plot_grid_points, scenario = "non-proportional")

sim_smooth <- sample.lsm(d = plot_grid_points, variogram_id = "exponential",
                         parameters = plot_params$params_for_sampling, dim = p, n_samples = 1)
nugget_noise_vis <- matrix(rnorm(N_tot_vis * p, mean = 0, sd = sqrt(nugget_vis)),
                           nrow = N_tot_vis, ncol = p)
true_realization <- sim_smooth$simulated_processes + nugget_noise_vis

# --- 3. Sample Training Data and Setup Kriging ---
train_indices_vis <- sample(1:N_tot_vis, N_train_vis)
train_coords_vis <- plot_grid_points[train_indices_vis, ]
train_data_vis <- true_realization[train_indices_vis, ]
residuals_train_vis <- train_data_vis 

train_params_vis <- plot_params$params_list[train_indices_vis]
test_params_vis <- plot_params$params_list

# --- 4. Perform Kriging Prediction for ALL THREE Models ---
cat("...building covariance matrices for all 3 models (this may take a moment)...\n")
# A) Training-Training Covariance Matrices
C_block_train_Op <- matrix(0, nrow = N_train_vis * p, ncol = N_train_vis * p)
C_trace_train_Correct <- matrix(0, nrow = N_train_vis, ncol = N_train_vis)
C_trace_train_Misspecified <- matrix(0, nrow = N_train_vis, ncol = N_train_vis)

for (i in 1:N_train_vis) {
  sigma_i_misspecified <- sqrt(sum(train_params_vis[[i]]$A^2))
  for (j in i:N_train_vis) {
    sigma_j_misspecified <- sqrt(sum(train_params_vis[[j]]$A^2))
    R_ns_ij <- compute_R_NS(train_coords_vis[i,], train_coords_vis[j,], train_params_vis[[i]], train_params_vis[[j]])
    
    C_ij_op <- train_params_vis[[i]]$A %*% t(train_params_vis[[j]]$A) * R_ns_ij
    
    idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
    C_block_train_Op[idx_i, idx_j] <- C_ij_op
    if (i != j) C_block_train_Op[idx_j, idx_i] <- t(C_ij_op)
    
    C_trace_train_Correct[i, j] <- sum(diag(C_ij_op))
    if (i != j) C_trace_train_Correct[j, i] <- C_trace_train_Correct[i, j]
    
    C_trace_train_Misspecified[i, j] <- sigma_i_misspecified * sigma_j_misspecified * R_ns_ij
    if (i != j) C_trace_train_Misspecified[j, i] <- C_trace_train_Misspecified[i, j]
  }
}
inv_C_block_train_Op <- chol2inv(chol(C_block_train_Op + diag(nugget_vis, N_train_vis * p)))
inv_C_trace_train_Correct <- chol2inv(chol(C_trace_train_Correct + diag(p * nugget_vis, N_train_vis)))
inv_C_trace_train_Misspecified <- chol2inv(chol(C_trace_train_Misspecified + diag(p * nugget_vis, N_train_vis)))

# B) Training-Prediction Cross-Covariance Matrices
C_block_traintest_Op <- matrix(0, nrow = N_train_vis * p, ncol = N_tot_vis * p)
C_trace_traintest_Correct <- matrix(0, nrow = N_train_vis, ncol = N_tot_vis)
C_trace_traintest_Misspecified <- matrix(0, nrow = N_train_vis, ncol = N_tot_vis)

for (i in 1:N_train_vis) {
  sigma_i_misspecified <- sqrt(sum(train_params_vis[[i]]$A^2))
  for (j in 1:N_tot_vis) {
    sigma_j_misspecified <- sqrt(sum(test_params_vis[[j]]$A^2))
    R_ns_ij <- compute_R_NS(train_coords_vis[i,], plot_grid_points[j,], train_params_vis[[i]], test_params_vis[[j]])
    
    C_ij_op <- train_params_vis[[i]]$A %*% t(test_params_vis[[j]]$A) * R_ns_ij
    
    idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
    C_block_traintest_Op[idx_i, idx_j] <- C_ij_op
    
    C_trace_traintest_Correct[i, j] <- sum(diag(C_ij_op))
    
    C_trace_traintest_Misspecified[i, j] <- sigma_i_misspecified * sigma_j_misspecified * R_ns_ij
  }
}

# C) Compute Weights and Make Predictions for all models
cat("...computing weights and making predictions for all models...\n")
weights_Op <- inv_C_block_train_Op %*% C_block_traintest_Op
weights_Correct <- inv_C_trace_train_Correct %*% C_trace_traintest_Correct
weights_Misspecified <- inv_C_trace_train_Misspecified %*% C_trace_traintest_Misspecified

pred_op_ns <- matrix(t(weights_Op) %*% as.vector(t(residuals_train_vis)), ncol = p, byrow = TRUE)
pred_correct <- t(weights_Correct) %*% residuals_train_vis
pred_misspecified <- t(weights_Misspecified) %*% residuals_train_vis

# --- 5. Prepare Data for Plot (3 components, 4 sources) ---
plot_df_combined <- data.frame(
  x = plot_grid_points[, 1],
  y = plot_grid_points[, 2],
  True_Z1 = true_realization[, 1],
  OpNS_Z1 = pred_op_ns[, 1],
  Correct_Z1 = pred_correct[, 1],
  Misspec_Z1 = pred_misspecified[, 1],
  True_Z2 = true_realization[, 2],
  OpNS_Z2 = pred_op_ns[, 2],
  Correct_Z2 = pred_correct[, 2],
  Misspec_Z2 = pred_misspecified[, 2],
  True_Z3 = true_realization[, 3],
  OpNS_Z3 = pred_op_ns[, 3],
  Correct_Z3 = pred_correct[, 3],
  Misspec_Z3 = pred_misspecified[, 3]
)

plot_df_long <- plot_df_combined %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = c("Source", "Component"),
    names_sep = "_",
    values_to = "Value"
  )

# Set factor levels for correct plot ordering and labeling
plot_df_long$Source <- factor(plot_df_long$Source,
                              levels = c("True", "OpNS", "Correct", "Misspec"))
levels(plot_df_long$Source) <- c("True Realization", "Op-NS (Cokriging)", 
                                 "Trace-NS (Correct)", "Trace-NS (Misspecified)")

train_points_df <- as.data.frame(train_coords_vis)
names(train_points_df) <- c("x", "y")

# --- 6. Create the Combined 3x4 Plot ---
cat("...generating final plot...\n")
prediction_plot <- ggplot(plot_df_long, aes(x = x, y = y, fill = Value)) +
  geom_tile() +
  geom_point(data = train_points_df, aes(x = x, y = y), 
             color = "black", size = 0.5, alpha = 0.6, inherit.aes = FALSE) +
  facet_grid(Component ~ Source) +
  scale_fill_viridis_c(option = "plasma") +
  coord_fixed(ratio = 1) +
  labs(
    title = "Visual Comparison of True vs. Predicted Fields (p=3)",
    subtitle = paste0("Non-Proportional scenario with ", N_train_vis, " training points (black dots)"),
    x = NULL, y = NULL, fill = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(prediction_plot)