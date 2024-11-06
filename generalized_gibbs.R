# Required Libraries
library(spatstat)
library(plotly)

# Function to simulate 2D Gibbs Point Process with optimization
gibbs_point_process_2d_optimized <- function(n, beta, r, iter = 500, window_size = 1) {
  # Create an initial random set of points within the unit square
  points <- matrix(runif(2 * n), ncol = 2) * window_size
  
  # Function to calculate pairwise distances efficiently using caching
  calculate_distances <- function(points) {
    dists <- as.matrix(dist(points))
    diag(dists) <- Inf  # Ignore self-interaction (distance = 0)
    return(dists)
  }
  
  # Metropolis-Hastings update step for Gibbs process with optimized distance calculation
  for (i in 1:iter) {
    dists <- calculate_distances(points)  # Precompute distances once per iteration
    
    for (j in 1:n) {
      # Select a point and propose a move within a small neighborhood
      current_point <- points[j, ]
      proposed_point <- current_point + runif(2, min = -0.05, max = 0.05)  # Small random move
      proposed_point <- pmin(pmax(proposed_point, 0), window_size)  # Keep point within bounds
      
      # Calculate distances between the proposed point and the others
      points[j, ] <- proposed_point
      new_dists <- calculate_distances(points)
      
      # Count the number of points within the repulsion distance `r`
      close_points_old <- sum(dists[j, ] < r)
      close_points_new <- sum(new_dists[j, ] < r)
      
      # Compute acceptance probability
      delta_U <- close_points_new - close_points_old  # Change in repulsion potential
      acceptance_prob <- min(1, exp(-beta * delta_U))
      
      # Accept or reject the proposed move
      if (runif(1) < acceptance_prob) {
        dists <- new_dists  # Update distances only if the move is accepted
      } else {
        points[j, ] <- current_point  # Reject move
      }
    }
  }
  
  return(points)
}

# Run optimized 2D Gibbs Point Process
result_2d_optimized <- gibbs_point_process_2d_optimized(n = 50, beta = 1, r = 1, iter = 500)

# Plot the result
plot(result_2d_optimized, xlab = "X", ylab = "Y", main = "Optimized 2D Gibbs Point Process", pch = 19, col = "blue")






# Function for Gibbs Process in 3D
gibbs_point_process_3d <- function(n, beta, r, domain_size = 1, iter = 5000) {
  # Initialize points randomly in the 3D domain
  points <- matrix(runif(3 * n), ncol = 3) * domain_size
  
  # Gibbs sampling loop
  for (i in 1:iter) {
    # Randomly select a point
    point_index <- sample(1:n, 1)
    current_point <- points[point_index, ]
    
    # Propose a new position within a small range around the point
    proposal <- current_point + runif(3, min = -0.05, max = 0.05)
    proposal <- pmin(pmax(proposal, 0), domain_size)  # Keep within bounds
    
    # Calculate the number of points within the distance `r` of the proposal
    distances <- sqrt(rowSums((points[-point_index, ] - proposal)^2))
    close_points <- sum(distances < r)
    
    # Calculate the probability of acceptance
    current_distances <- sqrt(rowSums((points[-point_index, ] - current_point)^2))
    current_close_points <- sum(current_distances < r)
    delta_U <- close_points - current_close_points  # Change in potential
    
    acceptance_prob <- min(1, exp(-beta * delta_U))
    
    # Accept or reject the new point
    if (runif(1) < acceptance_prob) {
      points[point_index, ] <- proposal
    }
  }
  
  # Return the points as a data frame for plotting
  return(data.frame(x = points[, 1], y = points[, 2], z = points[, 3]))
}

# Run 3D Gibbs point process
result_3d <- gibbs_point_process_3d(n = 100, beta = 1, r = 0.001, iter = 5000)

# Plot the 3D result
plot_ly(data = result_3d, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 2, color = 'blue', opacity = 0.5)) %>%
  layout(
    title = "3D Gibbs (Strauss) Point Process",
    scene = list(
      xaxis = list(title = 'X'),
      yaxis = list(title = 'Y'),
      zaxis = list(title = 'Z')
    )
  )
