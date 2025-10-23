# --- SECTION 1: SCRIPT SETUP ---
library(devtools)
# install_github("g-decarlo/FunctionalLSM") # Assumes this is already installed
library(LocallyStationaryModels)

# Install missing packages if necessary
packages <- c("dplyr", "ggplot2", "tidyr", "Matrix",
              "future.apply", "progressr", "gridExtra", "ellipse", "ggrepel", "cowplot")
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
library(gridExtra)
library(grid)
library(ellipse)
library(ggrepel)
library(cowplot)

# --- SECTION 2: CORE PARAMETER AND COVARIANCE FUNCTIONS ---

generate_scenario_parameters <- function(coords,
                                         scenario = c("non-proportional", "prop_center",
                                                      "prop_v1", "prop_v2", "prop_v3", "prop_v4")) {
  scenario <- match.arg(scenario)
  
  x <- coords[, 1]; y <- coords[, 2]
  n_pts <- nrow(coords)
  
  # Common Spatially Varying Parameters
  lambda1 <- (18 + 6 * x) / 100
  lambda2 <- (8 + 3 * x) / 10
  phi <- (pi / 6) * (2 + 2 * x)
  m1 <- 0 * x; m2 <- 0 * y # Zero mean
  
  # Spatially varying scalar intensity r(s)
  r_s <- 0.3 * (sqrt(x^2 + y^2) + 1)
  
  # Scenario-Specific Structural Matrix K(s)
  if (scenario == "non-proportional") {
    K11 <- 2 - x;  K12 <- 0
    K21 <- -2;     K22 <- 2 + y
  } else {
    fixed_coords <- switch(scenario,
                           "prop_center" = c(0, 0),
                           "prop_v1"     = c(-1, -1),
                           "prop_v2"     = c(1, -1),
                           "prop_v3"     = c(-1, 1),
                           "prop_v4"     = c(1, 1)
    )
    x_fixed <- fixed_coords[1]; y_fixed <- fixed_coords[2]
    K_fixed <- matrix(c(2 - x_fixed, -2, 0, 2 + y_fixed), nrow = 2)
    
    K11 <- rep(K_fixed[1, 1], n_pts); K12 <- rep(K_fixed[1, 2], n_pts)
    K21 <- rep(K_fixed[2, 1], n_pts); K22 <- rep(K_fixed[2, 2], n_pts)
  }
  
  # Final Coregionalization Matrix A(s) = r(s) * K(s)
  A11 <- r_s * K11; A12 <- r_s * K12
  A21 <- r_s * K21; A22 <- r_s * K22
  
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

compute_R_NS <- function(si, sj, params_i, params_j) {
  d <- 2
  Sigma_i <- params_i$Sigma_aniso; Sigma_j <- params_j$Sigma_aniso
  
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

run_simulation <- function(scenario, M_repetitions = 50, N_train = 100, p = 2, nugget = 1e-4, fixed_seed = 0) {
  # 1. Generate Parameters and Underlying Spatial Process
  grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 50)))
  all_params <- generate_scenario_parameters(grid_points, scenario)
  
  sim_result_smooth <- sample.lsm(
    d = grid_points, variogram_id = "exponential",
    parameters = all_params$params_for_sampling, dim = p, n_samples = M_repetitions, seed = fixed_seed
  )$simulated_processes
  
  nugget_noise <- matrix(rnorm(nrow(grid_points) * p * M_repetitions, mean = 0, sd = sqrt(nugget)),
                         nrow = nrow(grid_points), ncol = p * M_repetitions)
  all_sim_residuals <- sim_result_smooth + nugget_noise
  
  # 2. Pre-computation
  N_tot <- nrow(grid_points)
  train_indices <- sample(1:N_tot, N_train)
  test_indices <- setdiff(1:N_tot, train_indices)
  N_test <- length(test_indices)
  
  train_coords <- grid_points[train_indices, ]; test_coords <- grid_points[test_indices, ]
  train_params <- all_params$params_list[train_indices]; test_params <- all_params$params_list[test_indices]
  
  C_block_train_Op <- matrix(0, nrow = N_train * p, ncol = N_train * p)
  for (i in 1:N_train) {
    for (j in i:N_train) {
      R_ns_ij <- compute_R_NS(train_coords[i,], train_coords[j,], train_params[[i]], train_params[[j]])
      C_ij_op <- train_params[[i]]$A %*% t(train_params[[j]]$A) * R_ns_ij
      idx_i <- ((i - 1) * p + 1):(i * p); idx_j <- ((j - 1) * p + 1):(j * p)
      C_block_train_Op[idx_i, idx_j] <- C_ij_op
      if (i != j) C_block_train_Op[idx_j, idx_i] <- t(C_ij_op)
    }
  }
  inv_C_block_train_Op <- chol2inv(chol(C_block_train_Op + diag(nugget, N_train * p)))
  
  sigma_train <- sapply(train_params, function(param) sqrt(sum(diag(param$A %*% t(param$A)))))
  sigma_test  <- sapply(test_params,  function(param) sqrt(sum(diag(param$A %*% t(param$A)))))
  
  C_trace_train <- matrix(0, nrow = N_train * p, ncol = N_train * p)
  for (i in 1:N_train) {
    for (j in i:N_train) {
      R_ns_ij <- compute_R_NS(train_coords[i,], train_coords[j,], train_params[[i]], train_params[[j]])
      C_trace_ij <- sigma_train[i] * sigma_train[j] * R_ns_ij
      for(k in 1:p){
        C_trace_train[(i-1)*p+k, (j-1)*p+k] <- C_trace_ij
        if(i!=j) C_trace_train[(j-1)*p+k, (i-1)*p+k] <- C_trace_ij
      }
    }
  }
  inv_C_trace_train <- chol2inv(chol(C_trace_train + diag(nugget, N_train * p)))
  
  C_block_traintest_Op <- matrix(0, nrow = N_train * p, ncol = N_test * p)
  C_trace_traintest <- matrix(0, nrow = N_train * p, ncol = N_test * p)
  
  for (i in 1:N_train) {
    for (j in 1:N_test) {
      R_ns_ij <- compute_R_NS(train_coords[i,], test_coords[j,], train_params[[i]], test_params[[j]])
      idx_i <- ((i-1)*p + 1):(i*p); idx_j <- ((j-1)*p + 1):(j*p)
      C_block_traintest_Op[idx_i, idx_j] <- train_params[[i]]$A %*% t(test_params[[j]]$A) * R_ns_ij
      C_trace_ij_test <- sigma_train[i] * sigma_test[j] * R_ns_ij
      for(k in 1:p){
        C_trace_traintest[(i-1)*p+k, (j-1)*p+k] <- C_trace_ij_test
      }
    }
  }
  
  weights_Op <- inv_C_block_train_Op %*% C_block_traintest_Op
  weights_Trace <- inv_C_trace_train %*% C_trace_traintest
  
  results_list <- lapply(1:M_repetitions, function(m) {
    start_col <- (m - 1) * p + 1; end_col <- m * p
    sim_data <- all_sim_residuals[, start_col:end_col] + all_params$true_means
    train_data <- sim_data[train_indices, ]; test_data <- sim_data[test_indices, ]
    residuals_train <- as.vector(t(train_data - all_params$true_means[train_indices,]))
    
    pred_op_ns <- matrix(t(weights_Op) %*% residuals_train, ncol = p, byrow = TRUE) + all_params$true_means[test_indices,]
    mspe_op_ns <- mean(rowSums((test_data - pred_op_ns)^2))
    
    pred_trace_ns <- matrix(t(weights_Trace) %*% residuals_train, ncol = p, byrow = TRUE) + all_params$true_means[test_indices,]
    mspe_trace_ns <- mean(rowSums((test_data - pred_trace_ns)^2))
    
    data.frame(Repetition = m, MSPE_Op_NS = mspe_op_ns, MSPE_Trace_NS = mspe_trace_ns)
  })
  
  return(do.call(rbind, results_list))
}


# --- SECTION 4: VISUALIZATION OF SIMULATION SETUP ---

create_and_save_K_structure_plot <- function() {
  grid_dense <- as.matrix(expand.grid(x = seq(-1, 1, length.out=100),
                                      y = seq(-1, 1, length.out=100)))
  x <- grid_dense[, 1]; y <- grid_dense[, 2]
  
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
  
  proportional_points <- data.frame(
    x = c(0, -1, 1, -1, 1),
    y = c(0, -1, -1, 1, 1),
    label = c("(0,0)", "(-1,-1)", "(1,-1)", "(-1,1)", "(1,1)")
  )
  
  p_K <- ggplot(K_components_data, aes(x, y, fill=value)) +
    geom_raster() +
    facet_wrap(~component, nrow=1, labeller = label_parsed) +
    geom_point(data = proportional_points, aes(x, y), inherit.aes = FALSE,
               color = "#E69F00", size = 5, shape = 18) +
    geom_text_repel(data = proportional_points, aes(x, y, label = label), inherit.aes = FALSE,
                    color = "#E69F00", size = 5, fontface = "bold",
                    box.padding = 0.6, point.padding = 0.6,
                    segment.color = '#E69F00', segment.size = 0.7) +
    scale_fill_viridis_c() +
    coord_fixed() +
    labs(title="Components of Structural Matrix K(s) in Non-Proportional Scenario",
         subtitle="Orange diamonds mark the locations used to define the fixed K for the Proportional scenarios.",
         x="Coordinate x", y="Coordinate y", fill="Value") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size=18, face="bold"),
          plot.subtitle = element_text(size=14),
          strip.text = element_text(size=16, face="bold"),
          axis.title = element_text(size=16, face="bold"),
          legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=12))
  
  ggsave("K_structure_plot.png", plot = p_K, width = 14, height = 6, dpi = 300)
  
  cat("--- Simulation setup plot saved to K_structure_plot.png ---\n")
}

create_and_save_merged_realization_plot <- function(n_curves_to_plot = 24, p = 2, nugget = 1e-4, fixed_seed = 0) {
  
  scenario_name <- "non-proportional"
  grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 50), seq(-1, 1, length.out = 50)))
  
  t_grid <- seq(-pi, pi, length.out = 201)
  e1_t <- (1/sqrt(pi)) * cos(t_grid)
  e2_t <- (1/sqrt(pi)) * sin(t_grid)
  
  all_params <- generate_scenario_parameters(grid_points, scenario = scenario_name)
  sim_result_smooth <- sample.lsm(
    d = grid_points, variogram_id = "exponential",
    parameters = all_params$params_for_sampling, dim = p, n_samples = 1, seed = fixed_seed
  )$simulated_processes
  nugget_noise <- matrix(rnorm(nrow(grid_points) * p, mean = 0, sd = sqrt(nugget)),
                         nrow = nrow(grid_points), ncol = p)
  sim_coeffs <- sim_result_smooth + nugget_noise + all_params$true_means
  
  sample_indices <- sample(1:nrow(grid_points), n_curves_to_plot)
  sampled_coeffs <- sim_coeffs[sample_indices, ]
  sampled_coords_df <- as.data.frame(grid_points[sample_indices, ])
  colnames(sampled_coords_df) <- c("x", "y")
  
  sampled_r_s <- 0.3 * (sqrt(sampled_coords_df$x^2 + sampled_coords_df$y^2) + 1)
  sampled_coords_df$r_s_value <- sampled_r_s
  
  curves_data <- lapply(1:n_curves_to_plot, function(j) {
    x_t <- sampled_coeffs[j, 1] * e1_t + sampled_coeffs[j, 2] * e2_t
    data.frame(
      t = t_grid,
      value = x_t,
      curve_id = factor(j),
      r_s_value = sampled_r_s[j]
    )
  })
  
  plot_df <- do.call(rbind, curves_data)
  
  p_funcs <- ggplot(plot_df, aes(x = t, y = value, group = curve_id, color = r_s_value)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", name = "r(s) Value") +
    labs(
      title = "Functional Realizations (Non-Proportional)",
      x = expression(italic(t)),
      y = expression(italic(X[s](t)))
    ) +
    scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi),
      labels = c(expression(-pi), expression(-pi/2), "0", expression(pi/2), expression(pi))
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(size=14, face="bold"),
      legend.text = element_text(size=12),
      legend.key.width = unit(2, "cm"),
      panel.grid.major = element_line(linetype = "dashed", color = "grey85"),
      panel.grid.minor = element_blank()
    )
  
  grid_dense <- as.matrix(expand.grid(x = seq(-1, 1, length.out=100),
                                      y = seq(-1, 1, length.out=100)))
  r_s_data <- data.frame(
    x = grid_dense[, 1], y = grid_dense[, 2],
    r_s = 0.3 * (sqrt(grid_dense[, 1]^2 + grid_dense[, 2]^2) + 1)
  )
  
  p_r <- ggplot(r_s_data, aes(x, y, fill=r_s)) +
    geom_raster() +
    scale_fill_viridis_c(option = "magma", name = "r(s) Value") +
    geom_point(data = sampled_coords_df, aes(x = x, y = y, fill = r_s_value), # Corrected aes to use x, y for scatter
               shape = 23, size = 5, color = "white", stroke = 1.2) +
    coord_fixed() +
    labs(title="Scalar Intensity Field r(s)",
         x="Coordinate x", y="Coordinate y") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size=18, face="bold"),
          plot.subtitle = element_text(size=14),
          axis.title = element_text(size=16, face="bold"),
          axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=12),
          legend.key.width = unit(2, "cm"))
  
  legend <- get_legend(p_r)
  
  p_r <- p_r + theme(legend.position = "none")
  p_funcs <- p_funcs + theme(legend.position = "none")
  
  plot_row <- plot_grid(p_r, p_funcs, ncol = 2, labels = "AUTO", label_size = 20, align = "h")
  
  p_combined <- plot_grid(plot_row, legend, ncol = 1, rel_heights = c(1, 0.15))
  
  ggsave("realizations_and_intensity_map.png", plot = p_combined, width = 18, height = 9, dpi = 300, bg = "white")
  cat("--- Merged realization plot saved to realizations_and_intensity_map.png ---\n")
}


# --- SECTION 5: MAIN EXECUTION BLOCK ---

# ---!!--- CONTROL FLAG ---!!---
# Set to TRUE to run the full, long simulation.
# Set to FALSE to only generate the setup plots.
RUN_FULL_SIMULATION <- FALSE
# ---!!--------------------!!---

set_seed = 0
set.seed(set_seed) # for reproducibility

# Always generate the setup plots
create_and_save_K_structure_plot()
create_and_save_merged_realization_plot(fixed_seed = set_seed)

if (RUN_FULL_SIMULATION) {
  
  cat("\n--- RUN_FULL_SIMULATION is TRUE. Starting main simulation. ---\n")
  
  M_rep <- 350
  N_values <- c(20, 30, 40, 60, 150, 200)
  
  scenarios_to_run <- c(
    "non-proportional",
    "prop_v1",
    "prop_v2",
    "prop_center",
    "prop_v3",
    "prop_v4"
  )
  
  plan(multisession, workers = availableCores() - 4)
  
  experiment_grid <- expand.grid(
    scenario = scenarios_to_run,
    n_train = N_values,
    stringsAsFactors = FALSE
  )
  
  handlers(global = TRUE)
  handlers("progress")
  
  with_progress({
    p_progress <- progressor(steps = nrow(experiment_grid))
    
    all_results_list <- future_lapply(1:nrow(experiment_grid), function(i) {
      setting <- experiment_grid[i, ]
      p_progress(sprintf("Scenario=%s, N_train=%g", setting$scenario, setting$n_train))
      
      results_for_setting <- run_simulation(
        scenario = setting$scenario,
        M_repetitions = M_rep,
        N_train = setting$n_train,
        fixed_seed = set_seed
      )
      
      results_for_setting$Scenario <- setting$scenario
      results_for_setting$N_train <- setting$n_train
      return(results_for_setting)
      
    }, future.seed = TRUE)
  })
  cat("\n--- ALL SIMULATIONS COMPLETE ---\n")
  
  simulation_results <- do.call(rbind, all_results_list)
  
  
  # --- SECTION 6: ANALYSIS AND PLOTTING RESULTS ---
  
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
    mutate(Scenario_Label = factor(Scenario_Label, levels = c(
      "Non-Proportional\n",
      "Proportional @ (-1, -1)\n ",
      "Proportional @ (1, -1)\n ",
      "Proportional @ (0, 0)\n ",
      "Proportional @ (1, 1)\n ",
      "Proportional @ (-1, 1)\n "
    )))
  
  
  threshold <- 4*2^(-23)
  summary_stats <- simulation_results_labeled %>%
    group_by(Scenario_Label, N_train) %>%
    summarise(
      Mean_Diff_MSPE = mean(MSPE_Trace_NS - MSPE_Op_NS, na.rm = TRUE),
      SE_Diff_MSPE = sd(MSPE_Trace_NS - MSPE_Op_NS, na.rm = TRUE) / sqrt(n()),
      P_Value = tryCatch({
        t.test(MSPE_Trace_NS, MSPE_Op_NS, paired = TRUE, alternative = "greater", mu = threshold)$p.value
      }, error = function(e) { NA_real_ }),
      .groups = 'drop'
    ) %>%
    mutate(
      CI_Lower_MSPE = Mean_Diff_MSPE - 1.96 * SE_Diff_MSPE,
      CI_Upper_MSPE = Mean_Diff_MSPE + 1.96 * SE_Diff_MSPE
    )
  
  plot_data_diff <- summary_stats %>%
    select(Scenario_Label, N_train, Value = Mean_Diff_MSPE, CI_Lower = CI_Lower_MSPE, CI_Upper = CI_Upper_MSPE) %>%
    mutate(Metric = "MSPE Difference")
  
  plot_data_pval <- summary_stats %>%
    select(Scenario_Label, N_train, Value = P_Value) %>%
    mutate(Metric = "One-tailed t-test p-value", CI_Lower = NA, CI_Upper = NA)
  
  plot_data_final <- bind_rows(plot_data_diff, plot_data_pval) %>%
    mutate(Metric = factor(Metric, levels = c("MSPE Difference", "One-tailed t-test p-value")))
  
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
    scale_x_continuous(breaks = N_values) +
    labs(
      title = "Kriging Performance with Known Parameters: MSPE Difference and Statistical Significance",
      subtitle = "Positive MSPE difference favors Op-NS. Lower panel shows p-values for one-tailed paired t-test.",
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
  
  ggsave("mspe_significance_plot.png", plot = significance_plot, width = 14, height = 7, dpi = 300, bg = "white")
  
  cat("\n--- PLOT SAVED to mspe_significance_plot.png ---\n")
  print(significance_plot)
  
} else {
  
  cat("\n--- RUN_FULL_SIMULATION is FALSE. Main simulation skipped. ---\n")
  
}

