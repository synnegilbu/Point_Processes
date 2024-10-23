# Load necessary libraries
library(geoR)      # For Gaussian Random Field simulation
library(scales)    # For color scales
library(spatstat)  # For spatial analysis
library(rgl)       # For 3D visualization

# Generalized function to simulate Log-Gaussian Cox Process in d dimensions
simulate_lgcp <- function(d, window, grid_size, mean = 0, variance = 1, scale = 1, mark_type = "continuous", mark_range = c(0, 1)) {
  # Generate grid for spatial domain
  grid <- do.call(expand.grid, lapply(window, function(w) seq(w[1], w[2], length.out = grid_size)))
  coords <- as.matrix(grid)  # Coordinates of the grid points
  
  # Simulate Gaussian Random Field (GRF) using geoR
  grf_sim <- grf(n = nrow(coords), cov.pars = c(variance, scale), cov.model = "gaussian")
  
  # Extract GRF values and compute intensity field
  grf_values <- grf_sim$data  # GRF values
  intensity <- exp(mean + grf_values)  # Intensity function
  
  # Simulate the Poisson point process based on the intensity field
  lambda_max <- max(intensity)  # Maximum intensity
  n_points <- rpois(1, lambda_max * prod(sapply(window, diff)))  # Number of points
  
  # Generate random points uniformly in the d-dimensional window
  points <- matrix(runif(n_points * d), ncol = d)
  
  # Transform points into the correct ranges based on the window
  for (i in 1:d) {
    points[, i] <- points[, i] * diff(window[[i]]) + window[[i]][1]
  }
  
  # Assign intensity values to the random points
  intensity_at_points <- apply(points, 1, function(pt) {
    # Find the closest grid point to the random point
    grid_index <- which.min(rowSums((coords - pt)^2))
    return(intensity[grid_index])
  })
  
  # Apply thinning to keep points with probability proportional to their intensity
  keep_points <- runif(n_points) <= intensity_at_points / lambda_max
  final_points <- points[keep_points, , drop = FALSE]
  
  # Assign marks to the valid points
  if (mark_type == "continuous") {
    marks <- runif(nrow(final_points), min = mark_range[1], max = mark_range[2])
  } else {
    stop("Invalid mark type. Choose 'continuous'.")
  }
  
  # Return results
  if (d == 2) {
    return(list(points = data.frame(x = final_points[, 1], y = final_points[, 2], mark = marks),
                intensity = intensity, coords = coords))
  } else if (d == 3) {
    return(list(points = data.frame(x = final_points[, 1], y = final_points[, 2], z = final_points[, 3], mark = marks),
                intensity = intensity, coords = coords))
  } else {
    stop("Invalid dimension. Choose either 2 or 3.")
  }
}

# Example usage for 3D continuous marks
d <- 3  # Set the number of dimensions (3 for 3D)
window <- list(c(0, 10), c(0, 10), c(0, 10))  # Define the observation region for 3D
grid_size <- 20  # Size of the grid for Gaussian Random Field
mean <- 0        # Mean of the GRF
variance <- 1    # Variance of the GRF
scale <- 1       # Scale parameter for the covariance function

# Run the simulation with continuous marks
simulated_lgcp <- simulate_lgcp(d, window, grid_size, mean, variance, scale, mark_type = "continuous", mark_range = c(0, 100))

# Extract the points and marks
points <- simulated_lgcp$points
marks <- points$mark
shot_locations <- simulated_lgcp$coords

# Plot the simulated point pattern for 3D
open3d()  # Open a 3D plotting window
# Define color gradient based on continuous marks
color_scale <- scales::col_numeric(palette = "Blues", domain = c(0, 100))

# Plot points with colors based on their continuous marks
plot3d(points$x, points$y, points$z, col = color_scale(marks), size = 3,
       xlab = "X", ylab = "Y", zlab = "Z", main = "3D Log-Gaussian Cox Process with Continuous Marks")

# Add the shot locations (cluster centers) in a distinct color
points3d(shot_locations[, 1], shot_locations[, 2], shot_locations[, 3], col = 'black', size = 8)

# Add a grading scale (color legend) for continuous marks
color_legend <- function(values, colors, title, position = "topright") {
  # Create a gradient of colors
  breaks <- seq(min(values), max(values), length.out = length(colors) + 1)
  
  # Plot the legend
  legend(position, legend = round(breaks, 1), fill = colors, title = title, border = "black")
}

# Add the grading scale (color legend) to the plot
color_legend(values = seq(0, 100, length.out = 10), colors = color_scale(seq(0, 100, length.out = 10)), title = "Continuous Marks")
