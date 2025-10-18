# ---
# 1. SETUP AND LOAD DATA
# ---
library(LocallyStationaryModels)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  compositions, data.table, robCompositions, ggplot2, dplyr, maps, mapproj, cowplot, rworldmap
)
set.seed(42)

# Load the saved champion model and PCA results
tryCatch({
  load("champion_model.RData")
}, error = function(e) {
  stop("Could not find 'champion_model.RData'. Please run the 'rainfall_simulation.R' script first.")
})

# Load the FULL, UNFILTERED dataset for the final analysis
tryCatch({
  coords_full <- as.matrix(read.csv("~/FunctionalLSM/NS-fLMC/data/coordinatesrain.csv"))
  rain_obs_full <- read.csv("~/FunctionalLSM/NS-fLMC/data/rainobservations.csv")
  density_mat_full <- as.matrix(read.table("~/FunctionalLSM/NS-fLMC/data/density_matrix.prn", as.is = TRUE))
}, error = function(e) {
  stop("Data files not found. Please ensure data files are in a 'data/' subdirectory or adjust the paths.")
})

# ---
# 2. FULL DATA PREPROCESSING & MODEL FITTING
# ---
cat("--- Preprocessing full dataset ---\n")
# --- Functional Data (Densities) ---
density_mat_full <- matrix(unname(density_mat_full), nrow = dim(coords_full)[1], ncol = 60)
mode(density_mat_full) <- "numeric"
knots <- 10; xcp <- 1:60; data_range <- range(xcp); knot_range <- c(data_range[1] - 0.1, data_range[2] + 0.1)
knots_vector <- seq(from = knot_range[1], to = knot_range[2], length.out = knots)
smoothed_splines_full <- smoothSplines(
  k = 3, l = 2, alpha = 0.98, xcp = as.numeric(xcp), data = density_mat_full,
  knots = knots_vector, num_points = 256, prior = "sq"
)
z_scores_full <- prcomp(smoothed_splines_full$Y_clr, scale. = TRUE, center = TRUE)$x

prob_rain_full <- rowMeans(rain_obs_full > 0)
prob_rain_comp_full <- cbind(prob_rain_full, 1 - prob_rain_full)
logit_rain_full <- qlogis(prob_rain_full)

cat("--- Fitting final functional model on all data ---\n")
champion_epsilon <- champion_solution$epsilon
n_pcs <- ncol(champion_solution$initial_z)
anchor_points_full <- find_anchorpoints.lsm(coords_full, 12, TRUE)

variogram_full <- variogram.lsm(
  z = z_scores_full[, 1:n_pcs], d = coords_full, a = anchor_points_full$anchorpoints,
  epsilon = champion_epsilon, n_angles = 6, n_intervals = 24, dim = 1, kernel_id = "gaussian"
)
solution_full_functional <- findsolutions.lsm(
  variogram_full, remove_not_convergent = TRUE, lower.delta = 0.5,
  upper.bound = c(champion_epsilon, champion_epsilon, pi/2, 20, 50), lower.bound = c(2,2,0,1e-8,1e-8),
  initial.position = c(10, 10, pi/3, 10, 1), id = "exponentialnugget"
)

# --- Fit the scalar model for rain probability ---
cat("--- Fitting final scalar model for rain probability ---\n")
variogram_prob <- variogram.lsm(
  z = logit_rain_full, d = coords_full, a = anchor_points_full$anchorpoints,
  epsilon = champion_epsilon, n_angles = 8, n_intervals = 16, dim = 1, kernel_id = "gaussian"
)
solution_full_prob <- findsolutions.lsm(
  variogram_prob, remove_not_convergent = TRUE, lower.delta = 1,
  upper.bound = c(200, 200, pi/2, 8, 200), lower.bound = c(2,2,0,0,1e-8),
  initial.position = c(1, 1, pi/3, 0.5, 0.1), id = "exponentialnugget"
)

# ---
# 3. PREDICTION GRID & PLOTTING SETUP
# ---
cat("--- Setting up prediction grid and base map ---\n")
# Create a prediction grid
grid_pred <- expand.grid(
  lon = seq(from = -125, to = -67, length.out = 150),
  lat = seq(from = 25, to = 50, length.out = 150)
)

# Function to filter grid points to be within the continental US
coords2country <- function(points_df) {
  countriesSP <- getMap(resolution = 'low')
  pointsSP <- SpatialPoints(points_df, proj4string = CRS(proj4string(countriesSP)))
  indices <- over(pointsSP, countriesSP)
  return(indices$ADMIN)
}
grid_pred$country <- coords2country(grid_pred)
grid_pred_us <- filter(grid_pred, country == "United States")

# Get US map data for ggplot
us_map <- map_data("state")

# ---
# 4. PREDICTION AND BACK-TRANSFORMATION
# ---
cat("--- Predicting on grid ---\n")
# Predict functional scores and scalar probabilities
pred_functional <- predict.lsm(solution_full_functional, as.matrix(grid_pred_us[,c("lon","lat")]), FALSE)
pred_prob <- predict.lsm(solution_full_prob, as.matrix(grid_pred_us[,c("lon","lat")]), FALSE)

# Back-transform predicted scores to densities
cat("--- Back-transforming predictions to densities ---\n")
predicted_scores <- pred_functional$zpredicted
clr_predicted <- t(pca_result$center + t(predicted_scores %*% t(pca_result$rotation[,1:n_pcs])) * pca_result$scale)
predicted_densities <- clrInv(clr_predicted)

# Back-transform predicted logits to probabilities
predicted_rain_prob <- plogis(pred_prob$zpredicted)

# ---
# 5. GENERATE PLOTS
# ---
cat("--- Generating plots ---\n")

# --- Plot 1: Spatially Varying Parameters ---
smoothed_params <- smooth.lsm(solution_full_functional, as.matrix(grid_pred_us[,c("lon","lat")]))
plot_df <- cbind(grid_pred_us, smoothed_params$parameters)
colnames(plot_df) <- c("lon", "lat", "country", "lambda1", "lambda2", "phi", "sigma", "nugget")
plot_df$aniso_ratio <- plot_df$lambda1 / plot_df$lambda2

p_aniso <- ggplot() + geom_tile(data=plot_df, aes(x=lon, y=lat, fill=aniso_ratio)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_viridis_c(option="C") + labs(title="Anisotropy Ratio", fill="Ratio") + theme_void()

p_sigma <- ggplot() + geom_tile(data=plot_df, aes(x=lon, y=lat, fill=sigma)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_viridis_c(option="D") + labs(title="Spatially Varying Sigma", fill="Sigma") + theme_void()

ggsave("Anisotropy_and_Sigma.png", plot_grid(p_aniso, p_sigma, ncol=2), width=12, height=5)

# --- Plot 2: Anisotropy Ellipses on Map ---
ellipse_points <- function(center, lambda1, lambda2, phi, n=100) {
  angles <- seq(0, 2 * pi, length.out = n)
  ellipse_shape <- cbind(lambda1 * cos(angles), lambda2 * sin(angles))
  rotation_matrix <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow=2)
  rotated_ellipse <- t(rotation_matrix %*% t(ellipse_shape))
  return(data.frame(lon = center[1] + rotated_ellipse[,1], lat = center[2] + rotated_ellipse[,2]))
}

ellipses_df <- do.call(rbind, lapply(1:nrow(anchor_points_full$anchorpoints), function(i) {
  params <- solution_full_functional$solutions[i, ]
  center <- anchor_points_full$anchorpoints[i,]
  df <- ellipse_points(center, params[1]*0.2, params[2]*0.2, params[3]) # scaled for visibility
  df$group <- i
  return(df)
}))

p_ellipses <- ggplot() +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill="gray90", color="white") +
  geom_polygon(data=ellipses_df, aes(x=lon, y=lat, group=group), color="firebrick", fill="firebrick", alpha=0.4) +
  geom_point(data=as.data.frame(anchor_points_full$anchorpoints), aes(x=V1, y=V2), color="black", size=1) +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  labs(title="Anisotropy Ellipses at Anchor Points") + theme_void()

ggsave("Anisotropy_Ellipses.png", p_ellipses, width=8, height=5)

# --- Plot 3: Hazard Maps ---
# Define heavy rain as sqrt(mm) > sqrt(130), which is ~11.4
# The domain of densities is 0 to 60 over 256 points
heavy_rain_threshold_idx <- floor((sqrt(130) / 60) * 256)
prob_heavy_rain_conditional <- rowSums(predicted_densities[, heavy_rain_threshold_idx:256])
prob_heavy_rain_total <- prob_heavy_rain_conditional * predicted_rain_prob

plot_df$prob_heavy_cond <- prob_heavy_rain_conditional
plot_df$prob_rain <- predicted_rain_prob
plot_df$prob_heavy_total <- prob_heavy_rain_total

p_heavy_cond <- ggplot() + geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_heavy_cond)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_viridis_c(option="A", labels=scales::percent) + 
  labs(title="P(Rain > 130mm | Rain Occurred)", fill="Probability") + theme_void()

p_rain_occur <- ggplot() + geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_rain)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_viridis_c(option="B", labels=scales::percent) +
  labs(title="P(Rain Occurred)", fill="Probability") + theme_void()

p_heavy_total <- ggplot() + geom_tile(data=plot_df, aes(x=lon, y=lat, fill=prob_heavy_total)) +
  geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=c(-125, -67), ylim=c(25,50)) +
  scale_fill_viridis_c(option="magma", labels=scales::percent) +
  labs(title="Total Probability of Heavy Rain (>130mm)", fill="Probability") + theme_void()

ggsave("Hazard_Maps.png", plot_grid(p_heavy_cond, p_rain_occur, p_heavy_total, ncol=3), width=18, height=4)

cat("--- Analysis and plotting complete. Check your directory for PNG files. ---\n")
