# --- SECTION 1: SCRIPT SETUP ---
# Ensure required packages are installed
packages <- c("ggplot2", "dplyr", "tidyr", "Matrix", "cowplot", "LocallyStationaryModels")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(LocallyStationaryModels)
library(dplyr)
library(ggplot2)
library(tidyr)
library(Matrix)
library(cowplot)

# --- SECTION 2: CORE PARAMETER FUNCTION (from simulation script) ---

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
    K21 <- -0.5;     K22 <- 2 + y
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

# --- SECTION 3: GENERATE AND PLOT POINTWISE REALIZATION ---

set.seed(42)
p <- 2
nugget <- 1e-4

# 1. Define spatial grid and evaluation points
grid_points <- as.matrix(expand.grid(seq(-1, 1, length.out = 100), 
                                     seq(-1, 1, length.out = 100)))
t_eval <- c(-1, 1)

# 2. Define basis functions at t_eval
e1_at_t <- (1/sqrt(pi)) * cos(t_eval)
e2_at_t <- (1/sqrt(pi)) * sin(t_eval)

# 3. Generate parameters and one realization of coefficients
all_params <- generate_scenario_parameters(grid_points, scenario = "prop_v4")

sim_result_smooth <- sample.lsm(
  d = grid_points, variogram_id = "exponential",
  parameters = all_params$params_for_sampling, dim = p, n_samples = 1, seed = 42
)$simulated_processes

nugget_noise <- matrix(rnorm(nrow(grid_points) * p, mean = 0, sd = sqrt(nugget)),
                       nrow = nrow(grid_points), ncol = p)

sim_coeffs <- sim_result_smooth + nugget_noise + all_params$true_means

# 4. Calculate pointwise evaluations
# X_s(t) = C_1(s) * e_1(t) + C_2(s) * e_2(t)
X_s_t1 <- sim_coeffs[, 1] * e1_at_t[1] + sim_coeffs[, 2] * e2_at_t[1]
X_s_t2 <- sim_coeffs[, 1] * e1_at_t[2] + sim_coeffs[, 2] * e2_at_t[2]

# 5. Create plotting dataframe
plot_data <- data.frame(
  x = grid_points[, 1],
  y = grid_points[, 2],
  `-1` = X_s_t1,
  `1` = X_s_t2,
  check.names = FALSE
) %>%
  pivot_longer(
    cols = c("-1", "1"),
    names_to = "t_label",
    values_to = "value"
  ) %>%
  mutate(
    t_label = factor(t_label, levels = c("-1", "1"))
  )

# 6. Generate plot
# Use a diverging palette as the mean is 0
limit <- max(abs(plot_data$value)) * c(-1, 1)

p_pointwise <- ggplot(plot_data, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  facet_wrap(~t_label, labeller = label_bquote(italic(X[s])(t == .(as.character(t_label))))) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = limit, name = "Value") +
  coord_fixed() +
  labs(
    title = "Pointwise Evaluations of a Single Functional Realization",
    x = "Coordinate x",
    y = "Coordinate y"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "right"
  )

# 7. Save plot
ggsave("pointwise_realization_map.png", plot = p_pointwise, width = 14, height = 7, dpi = 300, bg = "white")

cat("--- Pointwise realization plot saved to pointwise_realization_map.png ---\n")

