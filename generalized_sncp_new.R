# Parameters:
#   d: dimensionality
#   lambda: cluster size intensity
#   omega: cluster dispersion
#   region: region bounds
#   n_shots: number of shots


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



