if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  LocallyStationaryModels, compositions, data.table, robCompositions,
  ggplot2, dplyr, maps, mapproj, cowplot, rworldmap, sp, RColorBrewer
)

set.seed(42)

# Define key hyperparameters for the analysis
map_granularity <- 150       # Resolution of the prediction grid
station_filter_epsilon <- 0.0001 # Min distance between stations (~10m)
num_anchor_points <- 12        # Number of anchor points for the NS model

# Load the champion model configuration and PCA results from the model selection step
tryCatch({
  load("champion_model.RData")
}, error = function(e) {
  stop("Could not find 'champion_model.RData'. Please run the model comparison script first.")
})

# Load the full, original dataset
tryCatch({
  coords_full <- as.matrix(read.csv("NS-fLMC/data/coordinatesrain.csv"))
  rain_obs_full <- read.csv("NS-fLMC/data/rainobservations.csv")
  density_mat_full <- as.matrix(read.table("NS-fLMC/data/density_matrix.prn", as.is = TRUE))
}, error = function(e) {
  stop("Data files not found. Please ensure they are in the 'NS-fLMC/data/' subdirectory.")
})


cat("--- Preprocessing and filtering full dataset ---\n")
n_initial <- nrow(coords_full)
indices_to_process <- 1:n_initial
indices_to_keep <- c()

# This loop efficiently filters stations that are too close to each other
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

# Load map data for the continental US, plus Canada and Mexico for context
us_map_full <- map_data("state")
us_map <- filter(us_map_full, !region %in% c("alaska", "hawaii"))
world_map <- map_data("world")
ca_mx_map <- filter(world_map, region %in% c("Canada", "Mexico"))

# Define consistent map coordinates for all plots to ensure proper alignment
map_xlim <- c(-127, -65)
map_ylim <- c(23, 51)


cache_file <- "prediction_cache.RData"
if (file.exists(cache_file)) {
  cat("--- Loading cached predictions ---\n")
  load(cache_file)
} else {
  cat("--- Predicting on grid (this may take a while) ---\n")
  pred_functional <- predict.lsm(solution_full_functional, grid_pred_us_matrix, plot_output = FALSE)
  pred_prob <- predict.lsm(solution_full_prob, grid_pred_us_matrix, plot_output = FALSE)
  
  cat("--- Back-transforming predictions ---\n")
  predicted_scores <- pred_functional$zpredicted
  clr_predicted <- t(pca_result$center + t(predicted_scores %*% t(pca_result$rotation[,1:n_pcs])) * pca_result$scale)
  predicted_densities_char <- clrInv(clr_predicted)
  predicted_densities <- apply(predicted_densities_char, 2, as.numeric)
  predicted_rain_prob <- plogis(pred_prob$zpredicted)
  
  cat(paste("--- Caching results to", cache_file, "---\n"))
  save(
    predicted_densities, predicted_rain_prob,
    file = cache_file
  )
}

cat("--- Generating final figures ---\n")
if (!dir.exists("plots")) {
  dir.create("plots")
}

larger_base_size <- 24

create_prediction_map <- function(plot_data, fill_var, palette, title, fill_label, base_font_size, is_prob = TRUE, legend_pos = "right") {
  p <- ggplot() +
    geom_tile(data = plot_data, aes(x = lon, y = lat, fill = {{fill_var}})) +
    geom_polygon(data = ca_mx_map, aes(x = long, y = lat, group = group), fill = "gray90", color = "white") +
    geom_polygon(data = us_map, aes(x = long, y = lat, group = group), fill = NA, color = "gray40") +
    coord_map(projection = "albers", lat0 = 39, lat1 = 45, xlim = map_xlim, ylim = map_ylim) +
    labs(title = title, x = "Longitude", y = "Latitude", fill = fill_label) +
    theme_bw(base_size = base_font_size)
  
  if (is_prob) {
    p <- p + scale_fill_distiller(palette = palette, direction = 1, labels = scales::percent, limits = c(0, NA))
  } else {
    p <- p + scale_fill_distiller(palette = palette, direction = 1)
  }
  
  if (legend_pos == "bottom") {
    p <- p + theme(legend.position = "bottom", 
                   legend.key.width = unit(2, "cm"), 
                   legend.title.position = "top", 
                   legend.title.hjust = 0.5)
  } else {
    p <- p + theme(legend.position = legend_pos)
  }
  
  return(p)
}


station_df <- data.frame(
  lon = coords_filtered[,1],
  lat = coords_filtered[,2],
  prob_rain = prob_rain_filtered
)
quantiles <- quantile(station_df$prob_rain, probs = c(0.2, 0.4, 0.6, 0.8, 0.95, 1))
selected_indices <- sapply(quantiles, function(q) which.min(abs(station_df$prob_rain - q)))
selected_stations_df <- station_df[selected_indices, ]

anchor_points_df <- as.data.frame(anchor_points_filtered$anchorpoints)
colnames(anchor_points_df) <- c("lon", "lat")

p_station_map <- ggplot() +
  geom_polygon(data = ca_mx_map, aes(x = long, y = lat, group = group), fill = "gray90", color = "white") +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group), fill = "gray80", color = "white") +
  geom_point(data = station_df, aes(x = lon, y = lat, color = prob_rain), size = 1.5, alpha = 0.8) +
  geom_point(data = anchor_points_df, aes(x = lon, y = lat), shape = 4, color = "black", size = 3, stroke = 1) +
  geom_point(data = selected_stations_df, aes(x = lon, y = lat), shape = 21, size = 4, stroke = 1.2, color = "red", fill = NA) +
  scale_color_distiller(palette = "Blues", direction = 1, labels = scales::percent) +
  coord_map(projection="albers", lat0=39, lat1=45, xlim=map_xlim, ylim=map_ylim) +
  labs(title="Station Rain Probability",
       x="Longitude", y="Latitude", color="P(Rain > 0)") +
  theme_bw(base_size = larger_base_size) +
  theme(legend.position = "none") 

x_sqrt_domain <- (1:256) * (60 / 256)
rainfall_mm_domain <- x_sqrt_domain^2

selected_densities <- smoothed_splines_filtered$Y[selected_indices, ]
plot_densities_df <- as.data.frame(t(selected_densities)) %>%
  mutate(rainfall_mm = rainfall_mm_domain) %>%
  tidyr::pivot_longer(cols = -rainfall_mm, names_to = "station_id", values_to = "density") %>%
  mutate(station_id = as.integer(gsub("V", "", station_id))) %>%
  left_join(
    selected_stations_df %>% mutate(station_id = 1:n()),
    by = "station_id"
  )

p_selected_densities <- ggplot(plot_densities_df, aes(x = rainfall_mm, y = density, group = station_id, color = prob_rain)) +
  geom_line(linewidth = 1.2) +
  scale_color_distiller(palette = "Blues", direction = 1, labels = scales::percent,
                        guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(xlim = c(0, 1000)) +
  labs(title="Selected Densities",
       x="Daily Rainfall (mm)", y="Density", color="P(Rain > 0)") +
  theme_bw(base_size = larger_base_size) +
  theme(legend.position = "none")

shared_legend_station_map <- get_legend(
  p_selected_densities + 
    theme(legend.position = "bottom", legend.key.width = unit(2.5, "cm"))
)

plots_row_station_map <- plot_grid(p_station_map, p_selected_densities, ncol = 2, align = 'h', axis = 'b')

p_station_map_combined <- plot_grid(plots_row_station_map, shared_legend_station_map, ncol = 1, rel_heights = c(1, 0.34))

ggsave("plots/Station_Map_and_Densities.png", p_station_map_combined, width = 20, height = 8.5, dpi = 300, bg = "white")


smoothed_params <- smooth.lsm(solution_full_functional, grid_pred_us_matrix)
plot_df <- cbind(grid_pred_us, smoothed_params$parameters)
colnames(plot_df) <- c("lon", "lat", "country", "lambda1", "lambda2", "phi", "sigma", "nugget")

p_sigma <- create_prediction_map(
  plot_data = plot_df,
  fill_var = sigma,
  palette = "YlOrRd",
  title = bquote("Spatially Varying Partial Sill (" * sigma * ")"),
  fill_label = expression(sigma),
  base_font_size = larger_base_size,
  is_prob = FALSE,
  legend_pos = "right"
)

plot_anisotropy_ellipses <- function(solution, anchor_points, map_data, border_data, title, ellipse_scale, base_font_size) {
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
    df <- ellipse_points(center, params[1]*ellipse_scale, params[2]*ellipse_scale, params[3])
    df$group <- i
    return(df)
  }))
  
  p <- ggplot() +
    geom_polygon(data=ca_mx_map, aes(x=long, y=lat, group=group), fill="gray90", color="white") +
    geom_polygon(data=us_map, aes(x=long, y=lat, group=group), fill="gray80", color="white") +
    geom_polygon(data=ellipses_df, aes(x=lon, y=lat, group=group), color="firebrick", fill="firebrick", alpha=0.5) +
    geom_point(data=anchor_points_df, aes(x=lon, y=lat), color="black", size=1.5) +
    coord_map(projection="albers", lat0=39, lat1=45, xlim=map_xlim, ylim=map_ylim) +
    labs(title=title, x="Longitude", y="Latitude") +
    theme_bw(base_size = base_font_size)
  
  return(p)
}

p_anisotropy <- plot_anisotropy_ellipses(
  solution_full_functional, anchor_points_filtered, us_map, ca_mx_map,
  title = "Anisotropy Ellipses",
  ellipse_scale = 0.35, base_font_size = larger_base_size
)

p_anisotropy_and_sill_combined <- plot_grid(p_anisotropy, p_sigma, ncol = 2, align = 'hv')

ggsave("plots/Anisotropy_and_Sill.png", p_anisotropy_and_sill_combined, width = 20, height = 7, dpi = 300, bg = "white")


heavy_rain_threshold_idx <- floor((sqrt(130) / 60) * 256)
prob_heavy_rain_conditional <- rowSums(predicted_densities[, heavy_rain_threshold_idx:256, drop = FALSE])
prob_heavy_rain_total <- prob_heavy_rain_conditional * predicted_rain_prob

plot_df$prob_heavy_cond <- prob_heavy_rain_conditional
plot_df$prob_rain <- predicted_rain_prob
plot_df$prob_heavy_total <- prob_heavy_rain_total

p_heavy_cond_no_legend <- create_prediction_map(
  plot_data = plot_df,
  fill_var = prob_heavy_cond,
  palette = "Blues",
  title = "P(Rain > 130mm | Rain Occurred)",
  fill_label = "Probability",
  base_font_size = larger_base_size,
  legend_pos = "none" # Remove legend
)

p_rain_occur_no_legend <- create_prediction_map(
  plot_data = plot_df,
  fill_var = prob_rain,
  palette = "Blues",
  title = "P(Rain Occurred)",
  fill_label = "Probability",
  base_font_size = larger_base_size,
  legend_pos = "none" # Remove legend
)

# Create shared legend for Hazard Maps
p_temp_legend <- create_prediction_map(
  plot_data = plot_df,
  fill_var = prob_rain,
  palette = "Blues",
  title = "P(Rain Occurred)",
  fill_label = "Probability",
  base_font_size = larger_base_size,
  legend_pos = "bottom"
)
shared_legend_hazard <- get_legend(p_temp_legend)

# Combine Hazard Maps
plots_row_hazard_maps <- plot_grid(
  p_heavy_cond_no_legend, 
  p_rain_occur_no_legend, 
  ncol = 2, 
  align = 'hv'
)

p_hazard_maps_combined <- plot_grid(
  plots_row_hazard_maps, 
  shared_legend_hazard, 
  ncol = 1, 
  rel_heights = c(1, 0.15) # Give legend 15% of vertical space
)

p_heavy_total <- create_prediction_map(
  plot_data = plot_df,
  fill_var = prob_heavy_total,
  palette = "Blues",
  title = "Total Probability of Heavy Rain (>130mm)",
  fill_label = "Probability",
  base_font_size = 14,
  legend_pos = "bottom"
)

# Save combined plot and the total hazard plot
ggsave("plots/Hazard_Maps_Combined.png", p_hazard_maps_combined, width = 20, height = 8.5, dpi = 300, bg = "white")
ggsave("plots/Total_Hazard_Map.png", p_heavy_total, width=10, height=7, dpi=300, bg = "white")

cat("--- Analysis and plotting complete. Check the 'plots' directory for PNG files. ---\n")

