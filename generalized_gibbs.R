
simulate_gibbs_with_covariates <- function(region, beta, interaction_radius, phi, n_iter, d = 3, covariate_field = NULL, covariate_coeff = 0) {
  
  region_volume <- prod(sapply(region, function(bounds) diff(bounds))) # Compute the volume of the region
  
  # Initialize the point configuration
  points <- matrix(nrow = 0, ncol = d)
  
  # Helper function to compute lambda(u, x)
  compute_lambda <- function(u, x, beta, interaction_radius, phi, covariate_field, covariate_coeff) {
    covariate_effect <- if (!is.null(covariate_field)) {
      exp(covariate_coeff * covariate_field(u))
    } else {
      1  # Default to no covariate effect
    }
    
    # Interaction penalty
    if (nrow(x) == 0) {
      return(beta * covariate_effect)
    }
    distances <- sqrt(rowSums((x - u)^2))
    interaction_penalty <- sum(phi(distances, interaction_radius))
    
    # Compute lambda
    lambda <- beta * covariate_effect * exp(-interaction_penalty)
    return(lambda)
  }
  
  # Helper function to compute r(u, x) for birth and death proposals
  compute_r <- function(proposed_point, current_points, add = TRUE) {
    if (add) {
      # Compute r(u, x) for adding a point
      lambda <- compute_lambda(proposed_point, current_points, beta, interaction_radius, phi, covariate_field, covariate_coeff)
      return((lambda * region_volume) / (nrow(current_points) + 1))
    } else {
      # Compute r(u, x) for removing a point
      remaining_points <- current_points[-proposed_point, , drop = FALSE]
      removed_point <- current_points[proposed_point, ]
      lambda <- compute_lambda(removed_point, remaining_points, beta, interaction_radius, phi, covariate_field, covariate_coeff)
      return(nrow(current_points) / (lambda * region_volume))
    }
  }
  
  accepted_points <- list()
  
  for (i in 1:n_iter) {
    proposal_type <- sample(c("add", "remove"), 1, prob = c(0.5, 0.5))
    acceptance_prob <- 0
    
    if (proposal_type == "add") {
      new_point <- sapply(region, function(bounds) runif(1, bounds[1], bounds[2]))
      r_ux <- compute_r(new_point, points, add = TRUE)
      acceptance_prob <- min(1, r_ux)
      
      if (runif(1) < acceptance_prob) {
        points <- rbind(points, new_point)  # Accept the new point
      }
      
    } else if (proposal_type == "remove" && nrow(points) > 0) {
      remove_idx <- sample(1:nrow(points), 1)
      r_ux <- compute_r(remove_idx, points, add = FALSE)
      acceptance_prob <- min(1, r_ux)
      
      if (runif(1) < acceptance_prob) {
        points <- points[-remove_idx, , drop = FALSE]  # Accept the removal
      }
    }
    
    accepted_points[[i]] <- points
  }
  
  # Return the final configuration and all accepted configurations
  list(final_points = points, all_points = accepted_points)
}

# Define a 2D region
region_2d <- list(c(0, 10), c(0, 10))  # Region is [0, 10] x [0, 10]

# Define a covariate field
covariate_field_2d <- function(point) {
  x <- point[1]
  y <- point[2]
  return(sin(x / 2) + cos(y / 3))  # Example: Oscillatory covariate field
}

# Define the potential function (hard-core process with soft repulsion)
phi <- function(r, interaction_radius) {
  ifelse(r < interaction_radius, 1 - r / interaction_radius, 0)  # Linear penalty
}

# Parameters for the Gibbs process
beta <- 50          # Baseline intensity
interaction_radius <- 10.0  # Interaction radius
n_iter <- 1000      # Number of iterations
covariate_coeff <- 0.5  # Weight of the covariate effect

# Run the simulation
result_2d <- simulate_gibbs_with_covariates(region_2d, beta, interaction_radius, phi, n_iter, d = 2,
                                            covariate_field = covariate_field_2d, covariate_coeff = covariate_coeff)

# Combine all points from accepted configurations
all_points_2d <- do.call(rbind, result_2d$all_points)

# Plot the final configuration
if (!is.null(all_points_2d) && nrow(all_points_2d) > 0) {
  plot(all_points_2d[, 1], all_points_2d[, 2], 
       pch = 16, col = "blue", 
       xlab = "X", ylab = "Y", 
       main = "Simulated Gibbs Process (2D) with Covariates",
       xlim = region_2d[[1]], ylim = region_2d[[2]])
} else {
  message("No points generated in the simulation.")
}


# Define a 3D region
region_3d <- list(c(0, 1), c(0, 1), c(0, 1))  # 3D region

# Define a covariate field
covariate_field_3d <- function(point) {
  x <- point[1]
  y <- point[2]
  z <- point[3]
  return(2 * x + y - z)  # Linear gradient covariate
}

# Define the potential function (hard-core process with soft repulsion)
phi <- function(r, interaction_radius) {
  ifelse(r < interaction_radius, 1 - r / interaction_radius, 0)  # Linear penalty
}

# Parameters for the Gibbs process
beta <- 1        # Baseline intensity
interaction_radius <- 5 # Interaction radius
n_iter <- 1000      # Number of iterations
covariate_coeff <- 2.0  # Strong covariate effect

# Run the simulation
result_3d <- simulate_gibbs_with_covariates(region_3d, beta, interaction_radius, phi, n_iter, d = 3,
                                            covariate_field = covariate_field_3d, covariate_coeff = covariate_coeff)

# Combine all points from accepted configurations
all_points_3d <- do.call(rbind, result_3d$all_points)

# Plot the final configuration
if (!is.null(all_points_3d) && nrow(all_points_3d) > 0) {
  library(rgl)
  plot3d(all_points_3d[, 1], all_points_3d[, 2], all_points_3d[, 3],
         col = "blue", size = 5,
         xlab = "X", ylab = "Y", zlab = "Z",
         main = "Simulated Gibbs Process (3D) with Covariates")
} else {
  message("No points generated in the simulation.")
}

