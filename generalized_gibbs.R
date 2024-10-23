library(rgl)

# Function to simulate a Gibbs process in d dimensions and plot the result
simulate_gibbs_process <- function(d = 2, n_points = 100, n_iter = 1000, 
                                   interaction_range = 0.1, beta = -0.5, 
                                   window_size = 1, plot_result = TRUE) {
  
  # Initialize the point coordinates randomly within the observation window
  points <- matrix(runif(n_points * d, min = 0, max = window_size), ncol = d)
  
  # Gibbs process parameters
  interaction_radius <- interaction_range
  beta_param <- beta
  
  # Function to calculate the potential energy of a configuration
  potential_function <- function(points) {
    total_potential <- 0
    n <- nrow(points)
    
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        distance <- sqrt(sum((points[i, ] - points[j, ])^2))
        if (distance < interaction_radius) {
          total_potential <- total_potential + beta_param
        }
      }
    }
    return(total_potential)
  }
  
  # Metropolis-Hastings MCMC loop
  for (iter in 1:n_iter) {
    # Propose a new point configuration
    idx <- sample(1:n_points, 1) # Choose a random point to move
    old_point <- points[idx, ]
    new_point <- runif(d, min = 0, max = window_size) # New location in the window
    
    # Temporary points matrix for the new configuration
    new_points <- points
    new_points[idx, ] <- new_point
    
    # Calculate potentials
    old_potential <- potential_function(points)
    new_potential <- potential_function(new_points)
    
    # Calculate acceptance probability
    accept_prob <- exp(old_potential - new_potential)
    
    # Accept or reject the new configuration based on the acceptance probability
    if (runif(1) < accept_prob) {
      points[idx, ] <- new_point
    }
  }
  
  # Plot the result if requested
  if (plot_result) {
    plot_points(points, d, window_size)
  }
  
  return(points)
}

# Function to plot the points based on the dimension
plot_points <- function(points, d, window_size) {
  if (d == 2) {
    plot(points[,1], points[,2], xlim = c(0, window_size), ylim = c(0, window_size), 
         xlab = "X", ylab = "Y", main = "2D Gibbs Process Simulation", pch = 19)
  } else if (d == 3) {
    plot3d(points[,1], points[,2], points[,3], xlim = c(0, window_size), 
           ylim = c(0, window_size), zlim = c(0, window_size), 
           xlab = "X", ylab = "Y", zlab = "Z", 
           main = "3D Gibbs Process Simulation", col = "blue", size = 5)
  } else {
    message("Plotting is not supported for dimensions higher than 3.")
  }
}

# Example usage
set.seed(42)
result_points <- simulate_gibbs_process(d = 3, n_points = 50, n_iter = 5000, 
                                        interaction_range = 3, beta = -1, 
                                        window_size = 1, plot_result = TRUE)


# Example usage for 2D Gibbs Process
set.seed(42)
result_points_2d <- simulate_gibbs_process(d = 2, n_points = 100, n_iter = 5000, 
                                           interaction_range = 1, beta = -1, 
                                           window_size = 1, plot_result = TRUE)

