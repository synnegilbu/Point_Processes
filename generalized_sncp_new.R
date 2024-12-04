library(rgl)

simulate_sncp <- function(d, lambda, omega, region, n_shots) {
  
  # Generate shot locations randomly within the specified region
  shot_locations <- matrix(runif(n_shots * d, 
                                 min = unlist(lapply(region, `[`, 1)), 
                                 max = unlist(lapply(region, `[`, 2))),
                           ncol = d)
  
  # Initialize a matrix to store the generated points
  points <- matrix(NA, nrow = 0, ncol = d)
  
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
    
    # Add the valid points to the list of all points
    points <- rbind(points, valid_points)
  }
  
  # Create a data frame for points
  point_columns <- paste0("dim", 1:d)
  points_df <- as.data.frame(points)
  colnames(points_df) <- point_columns
  
  # Return the results
  return(list(points = points_df, shot_locations = shot_locations))
}

# Example Usage
region <- list(c(0, 10), c(0, 10))  # Define region in 4D space
lambda <- 15           # Mean number of offspring per shot
omega <- 0.5           # Standard deviation of offspring dispersion
n_shots <- 5           # Number of shot centers
d <- 2                 # Number of dimensions

# Simulate process
result <- simulate_sncp(d, lambda, omega, region, n_shots)

# Check the structure of the result
str(result)

# If d = 2 or d = 3, you could visualize the result:

  plot(result$points[, 1:2], col = "blue", pch = 16, main = "Shot Noise Cox Process (2D)")
  points(result$shot_locations, col = "red", pch = 4)



region <- list(c(0, 10), c(0, 10), c(0, 10))  # Define region in 3D space
lambda <- 15           # Mean number of offspring per shot
omega <- 0.5           # Standard deviation of offspring dispersion
n_shots <- 5           # Number of shot centers
d <- 3                 # Number of dimensions

# Simulate process
result <- simulate_sncp(d, lambda, omega, region, n_shots)

# Visualization for 3D

open3d()
plot3d(result$points, col = "blue", size = 3, main = "Shot Noise Cox Process (3D)")
points3d(result$shot_locations, col = "red", size = 6)  # Shot locations in red


