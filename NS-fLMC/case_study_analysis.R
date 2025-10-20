# CASE STUDY ANALYSIS AND PLOTTING
#
# PURPOSE:
# This script loads the champion model selected by the model comparison script,
# refits it on the full (filtered) dataset, and generates all the
# publication-quality figures for the case study, including parameter maps,
# anisotropy ellipses, and the final hazard maps.
#
# OUTPUT:
# A series of high-resolution PNG images saved to the 'plots' subdirectory.

# ---
# 1. SETUP AND LOAD DATA
# ---
# Required packages
library(LocallyStationaryModels)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  compositions, data.table, robCompositions, ggplot2, dplyr, maps, 
  mapproj, cowplot, rworldmap, sp, RColorBrewer
)

# Ensure reproducibility
set.seed(42)

# Define key hyperparameters for the analysis
map_granularity <- 150          # Resolution of the prediction grid (e.g., 50x50)
station_filter_epsilon <- 0.0001  # Min distance between stations in degrees (~10m)
num_anchor_points <- 12        # Number of anchor points for the NS model

# Load the champion model configuration and PCA results
tryCatch({
  load("champion_model.RData")
}, error = function(e) {
  stop("Could not find 'champion_model.RData'. Please run '1_model_comparison.R' first.")
})


# Load the full, original dataset
tryCatch({
  coords_full <- as.matrix(read.csv("NS-fLMC/data/coordinatesrain.csv"))
  rain_obs_full <- read.csv("NS-fLMC/data/rainobservations.csv")
  density_mat_full <- as.matrix(read.table("NS-fLMC/data/density_matrix.prn", as.is = TRUE))
}, error = function(e) {
  stop("Data files not found. Please ensure they are in a 'data/' subdirectory.")
})


# ---
# 2. STATION FILTERING & DATA PREPROCESSING
# ---
cat("--- Preprocessing and filtering full dataset ---\n")
n_initial <- nrow(coords_full)
indices_to_process <- 1:n_initial
indices_to_keep <- c()

while(length(indices_to_process) > 0) {
  current_index <- indices_to_process[1]
  indices_to_keep <- c(indices_to_keep, current_index)
  if (length(indices_to_process) > 1) {
    current_coords <- coords_full[current_index, , drop = FALSE]
    remaining_indices <- indices_to_process[-1]
    remaining_coords <- coords_full[remaining_indices, , drop = FALSE]
    distances <- sqrt(rowSums(sweep(remaining_coords, 2, current_coords, "-")^2))
    close_indices <- remaining_indices[which(distances <= station_filter_epsilon)]
    indices_to_remove_from_pool <- c(current_index, close_indices)
  } else {
    indices_to_remove_from_pool <- current_index
  }
  indices_to_process <- setdiff(indices_to_process, indices_to_remove_from_pool)
}

cat(sprintf("Finished filtering. %d stations remaining.\n", length(indices_to_keep)))

coords_filtered <- coords_full[indices_to_keep, ]
rain_obs_filtered <- rain_obs_full[indices_to_keep, ]
density_mat_filtered <- density_mat_full[indices_to_keep, ]

# Preprocess functional data (densities)
density_mat_filtered <- matrix(unname(density_mat_filtered), nrow = dim(coords_filtered)[1], ncol = 60)
mode(density_mat_filtered) <- "numeric"
knots <- 10; xcp <- 1:60
knots_vector <- seq(from = range(xcp)[1] - 0.1, to = range(xcp)[2] + 0.1, length.out = knots)
smoothed_splines_filtered <- smoothSplines(
  k = 3, l = 2, alpha = 0.98, xcp = as.numeric(xcp), data = density_mat_filtered,
  knots = knots_vector, num_points = 256, prior = "sq"
)
z_scores_filtered <- predict(pca_result, newdata=smoothed_splines_filtered$Y_clr)

# Preprocess scalar data (rain probability)
prob_rain_filtered <- rowMeans(rain_obs_filtered > 0)
logit_rain_filtered <- qlogis(prob_rain_filtered)

# ---
# 3. MODEL FITTING ON FINAL DATASET
# ---
cat("--- Fitting final functional model on filtered data ---\n")
champion_epsilon <- champion_solution$epsilon
n_pcs <- ncol(champion_solution$initial_z)
anchor_points_filtered <- find_anchorpoints.lsm(coords_filtered, num_anchor_points, TRUE)

# Fit the champion functional model on the full filtered dataset
variogram_full <- variogram.lsm(
  z = z_scores_filtered[, 1:n_pcs], d = coords_filtered, a = anchor_points_filtered$anchorpoints,
  epsilon = champion_epsilon, n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "gaussian"
)
solution_full_functional <- findsolutions.lsm(
  variogram_full, remove_not_convergent = TRUE, lower.delta = 0.1,
  upper.bound = c(champion_epsilon, champion_epsilon, pi / 2, 20, 50), lower.bound = c(2, 2, 0, 1e-8,1e-8),
  initial.position = c(10, 10, pi / 3, 10, 1), id = "exponentialnugget"
)

cat("--- Fitting final scalar model for rain probability ---\n")
# Fit the scalar model using the same champion epsilon for consistency
variogram_prob <- variogram.lsm(
  z = logit_rain_filtered, d = coords_filtered, a = anchor_points_filtered$anchorpoints,
  epsilon = champion_epsilon, n_angles = 8, n_intervals = 16, dim = 1, kernel_id = "gaussian"
)
solution_full_prob <- findsolutions.lsm(
  variogram_prob, remove_not_convergent = TRUE, lower.delta = 0.1,
  upper.bound = c(2*champion_epsilon, 2*champion_epsilon, pi/2, 8, 200), lower.bound = c(1,1,0,.1,1e-8),
  initial.position = c(5, 5, pi/3, 0.8, 0), id = "exponentialnugget"
)

# ---
# 4. PREDICTION AND PLOTTING
# ---
cat("--- Setting up prediction grid and base map ---\n")
grid_pred <- expand.grid(
  lon = seq(from = -125, to = -67, length.out = map_granularity),
  lat = seq(from = 25, to = 50, length.out = map_granularity)
)

coords2country <- function(points_df) {
  countriesSP <- getMap(resolution = 'low')
  pointsSP <- SpatialPoints(points_df, proj4string = CRS(proj4string(countriesSP)))
  indices <- over(pointsSP, countriesSP)
  return(indices$ADMIN)
}
grid_pred$country <- coords2country(grid_pred)
grid_pred_us <- filter(grid_pred, country == "United States of America")
grid_pred_us_matrix <- as.matrix(grid_pred_us[, c("lon", "lat")])

us_map_full <- map_data("state")
us_map <- filter(us_map_full, !region %in% c("alaska", "hawaii"))

cat("--- Predicting on grid and back-transforming ---\n")
pred_functional <- predict.lsm(solution_full_functional, grid_pred_us_matrix, plot_output = FALSE)
pred_prob <- predict.lsm(solution_full_prob, grid_pred_us_matrix, plot_output = FALSE)

predicted_scores <- pred_functional$zpredicted
clr_predicted <- t(pca_result$center + t(predicted_scores %*% t(pca_result$rotation[,1:n_pcs])) * pca_result$scale)
predicted_densities_char <- clrInv(clr_predicted)
predicted_densities <- apply(predicted_densities_char, 2, as.numeric)

predicted_rain_prob <- plogis(pred_prob$zpredicted)

cat("--- Generating final figures ---\n")
if (!dir.exists("plots")) {
  dir.create("plots")
}

# --- Spatially Varying Parameter Maps ---
smoothed_params <- smooth.lsm(solution_full_functional, grid_pred_us_matrix)
plot_df <- cbind(grid_pred_us, smoothed_params$parameters)
colnames(plot_df) <- c("lon", "lat", "country", "lambda1", "lambda2", "phi", "sigma", "nugget")
plot_df$aniso_ratio <- plot_df$lambda1 / plot_df$lambda2

p_aniso <- ggplot() +
  geom_tile(data=plot_df, aes(x=lon, y=lat, fill=aniso_ratio)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="gray40") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  labs(title="Anisotropy Ratio (Functional Model)", x="Longitude", y="Latitude", fill="Ratio") +
  theme_bw(base_size = 14)

p_sigma <- ggplot() +
  geom_tile(data=plot_df, aes(x=lon, y=lat, fill=sigma)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="gray40") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(title=bquote("Spatially Varying" ~ sigma ~ "(Functional Model)"), x="Longitude", y="Latitude", fill=expression(sigma)) +
  theme_bw(base_size = 14)

ggsave("plots/Parameter_Anisotropy_Functional.png", p_aniso, width=10, height=7, dpi=300)
ggsave("plots/Parameter_Sigma_Functional.png", p_sigma, width=10, height=7, dpi=300)

# --- Anisotropy Ellipse Maps ---
plot_anisotropy_ellipses <- function(solution, anchor_points, map_data, title, filename, ellipse_scale) {
  anchor_coords <- anchor_points$anchorpoints
  smoothed_params_at_anchors <- smooth.lsm(solution, anchor_coords)$parameters
  
  anchor_points_df <- as.data.frame(anchor_coords)
  colnames(anchor_points_df) <- c("lon", "lat")
  
  ellipse_points <- function(center, lambda1, lambda2, phi, n=100) {
    angles <- seq(0, 2 * pi, length.out = n)
    ellipse_shape <- cbind(lambda1 * cos(angles), lambda2 * sin(angles))
    rotation_matrix <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow=2)
    rotated_ellipse <- t(rotation_matrix %*% t(ellipse_shape))
    return(data.frame(lon = center[1] + rotated_ellipse[,1], lat = center[2] + rotated_ellipse[,2]))
  }
  
  ellipses_df <- do.call(rbind, lapply(1:nrow(anchor_coords), function(i) {
    params <- smoothed_params_at_anchors[i, ]
    center <- anchor_coords[i, ]
    df <- ellipse_points(center, params[1]*ellipse_scale, params[2]*ellipse_scale, params[3]) # Use passed scale
    df$group <- i
    return(df)
  }))
  
  p <- ggplot() +
    geom_polygon(data=map_data, aes(x=long, y=lat, group=group), fill="gray90", color="white") +
    geom_polygon(data=ellipses_df, aes(x=lon, y=lat, group=group), color="firebrick", fill="firebrick", alpha=0.5) +
    geom_point(data=anchor_points_df, aes(x=lon, y=lat), color="black", size=1.5) +
    coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
    labs(title=title, x="Longitude", y="Latitude") +
    theme_bw(base_size = 14)
  
  ggsave(filename, p, width=10, height=7, dpi=300)
}

plot_anisotropy_ellipses(solution_full_functional, anchor_points_filtered, us_map,
                         "Anisotropy Ellipses (Functional Model)", "plots/Anisotropy_Ellipses_Functional.png",
                         ellipse_scale = 0.35)
plot_anisotropy_ellipses(solution_full_prob, anchor_points_filtered, us_map,
                         "Anisotropy Ellipses (Probability Model)", "plots/Anisotropy_Ellipses_Prob.png",
                         ellipse_scale = 0.15)


# --- Individual Hazard Maps ---
heavy_rain_threshold_idx <- floor((sqrt(130) / 60) * 256)
prob_heavy_rain_conditional <- rowSums(predicted_densities[, heavy_rain_threshold_idx:256])
prob_heavy_rain_total <- prob_heavy_rain_conditional * predicted_rain_prob

plot_df$prob_heavy_cond <- prob_heavy_rain_conditional
plot_df$prob_rain <- predicted_rain_prob
plot_df$prob_heavy_total <- prob_heavy_rain_total

p_heavy_cond <- ggplot() +
  geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_heavy_cond)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="gray40") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_distiller(palette = "Blues", direction = 1, labels=scales::percent) +
  labs(title="P(Rain > 130mm | Rain Occurred)", x="Longitude", y="Latitude", fill="Probability") +
  theme_bw(base_size = 14)

p_rain_occur <- ggplot() +
  geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_rain)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="gray40") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_distiller(palette = "Greens", direction = 1, labels=scales::percent) +
  labs(title="P(Rain Occurred)", x="Longitude", y="Latitude", fill="Probability") +
  theme_bw(base_size = 14)

p_heavy_total <- ggplot() +
  geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_heavy_total)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="gray40") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, labels=scales::percent) +
  labs(title="Total Probability of Heavy Rain (>130mm)", x="Longitude", y="Latitude", fill="Probability") +
  theme_bw(base_size = 14)

ggsave("plots/Hazard_Map_Conditional.png", p_heavy_cond, width=10, height=7, dpi=300)
ggsave("plots/Hazard_Map_Occurrence.png", p_rain_occur, width=10, height=7, dpi=300)
ggsave("plots/Hazard_Map_Total.png", p_heavy_total, width=10, height=7, dpi=300)

cat("--- Analysis and plotting complete. Check your directory for PNG files. ---\n")

