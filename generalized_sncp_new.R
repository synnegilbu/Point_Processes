# Load necessary libraries
library(scales)
library(spatstat)

# Generalized function to simulate Shot Noise Cox Process in d dimensions
simulate_sncp <- function(d, lambda, kappa, omega, region, n_shots, mark_type = "categorical", categories = NULL, mark_range = c(0, 1)) {
  # Check if region is appropriately defined
  if (length(region) != d || !all(sapply(region, length) == 2)) {
    stop("Region must specify min and max for each of the dimensions as a list of length 2.")
  }
  
  # Generate shot locations randomly within the specified region
  shot_locations <- matrix(runif(n_shots * d, 
                                 min = unlist(lapply(region, `[`, 1)), 
                                 max = unlist(lapply(region, `[`, 2))),
                           ncol = d)
  
  # Initialize a matrix to store the generated points and their marks
  points <- matrix(NA, nrow = 0, ncol = d)
  marks <- c()
  
  # For each shot location, generate a number of points around it
  for (i in 1:n_shots) {
    # Number of points to generate around this cluster center
    n_points <- rpois(1, lambda)
    
    # Generate points using a Gaussian dispersal kernel around the shot location
    new_points <- matrix(NA, nrow = n_points, ncol = d)
    for (j in 1:d) {
      new_points[, j] <- rnorm(n_points, mean = shot_locations[i, j], sd = omega)  # Gaussian spread
    }
    
    # Keep only points within the specified region
    valid_points <- new_points[apply(new_points, 1, function(pt) {
      all(pt >= unlist(lapply(region, `[`, 1)) & pt <= unlist(lapply(region, `[`, 2)))
    }), ]
    
    # Generate marks for the valid points
    if (mark_type == "categorical") {
      if (is.null(categories)) {
        stop("For categorical marks, please provide a vector of categories.")
      }
      new_marks <- sample(categories, nrow(valid_points), replace = TRUE)
    } else if (mark_type == "continuous") {
      new_marks <- runif(nrow(valid_points), min = mark_range[1], max = mark_range[2])
    } else {
      stop("Invalid mark type. Choose 'categorical' or 'continuous'.")
    }
    
    # Add the valid points and their marks to the list of all points
    points <- rbind(points, valid_points)
    marks <- c(marks, new_marks)
  }
  
  # Create the data frame for points, handling dimensions
  if (d == 2) {
    return(list(points = data.frame(x = points[, 1], y = points[, 2], mark = marks),
                shot_locations = shot_locations))
  } else if (d == 3) {
    return(list(points = data.frame(x = points[, 1], y = points[, 2], z = points[, 3], mark = marks),
                shot_locations = shot_locations))
  } else {
    stop("Invalid dimension. Choose either 2 or 3.")
  }
}

# Example usage for 2D continuous marks
d <- 2  # Set the number of dimensions (2 for 2D)
lambda <- 50    # Mean number of points per cluster
kappa <- 10     # Intensity of the Poisson process of cluster centers
omega <- 1.0    # Bandwidth of the cluster dispersal kernel
n_shots <- rpois(1, kappa)  # Number of shot locations determined by a Poisson process

# Define the observation region for 2D
region <- list(c(0, 10), c(0, 10))  # 2D region

# Run the simulation with continuous marks
simulated_sncp <- simulate_sncp(d, lambda, kappa, omega, region, n_shots, mark_type = "continuous", mark_range = c(0, 100))

# Extract the points and marks
points <- simulated_sncp$points
marks <- points$mark
shot_locations <- simulated_sncp$shot_locations

# Plot the simulated point pattern for 2D
plot(points$x, points$y, col = scales::col_numeric(palette = "Blues", domain = NULL)(marks), pch = 16,
     xlab = "X", ylab = "Y", main = "2D Shot Noise Cox Process with Continuous Marks")
points(shot_locations[, 1], shot_locations[, 2], col = "black", pch = 19, cex = 1.5)

# Add a grading scale (color legend) for continuous marks
color_legend <- function(values, colors, title, position = "topright") {
  # Create a gradient of colors
  breaks <- seq(min(values), max(values), length.out = length(colors) + 1)
  
  # Plot the legend
  legend(position, legend = round(breaks, 1), fill = colors, title = title, border = "black")
}

# Add the grading scale (color legend) to the plot
color_scale <- scales::col_numeric(palette = "Blues", domain = c(0, 100))
color_legend(values = seq(0, 100, length.out = 10), colors = color_scale(seq(0, 100, length.out = 10)), title = "Continuous Marks")



# Example usage for 3D continuous marks
d <- 3  # Set the number of dimensions (3 for 3D)
lambda <- 50    # Mean number of points per cluster
kappa <- 10     # Intensity of the Poisson process of cluster centers
omega <- 1.0    # Bandwidth of the cluster dispersal kernel
n_shots <- rpois(1, kappa)  # Number of shot locations determined by a Poisson process

# Define the observation region for 3D
region <- list(c(0, 10), c(0, 10), c(0, 10))  # 3D region

# Run the simulation with continuous marks
simulated_sncp <- simulate_sncp(d, lambda, kappa, omega, region, n_shots, mark_type = "continuous", mark_range = c(0, 100))

# Extract the points and marks
points <- simulated_sncp$points
marks <- points$mark
shot_locations <- simulated_sncp$shot_locations

# Plot the simulated point pattern for 3D
open3d()  # Open a 3D plotting window
# Define color gradient based on continuous marks
color_scale <- scales::col_numeric(palette = "Blues", domain = c(0, 100))

# Plot points with colors based on their continuous marks
plot3d(points$x, points$y, points$z, col = color_scale(marks), size = 3,
       xlab = "X", ylab = "Y", zlab = "Z", main = "3D Shot Noise Cox Process with Continuous Marks")

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

