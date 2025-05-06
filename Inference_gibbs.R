library(rstan)

simulate_gibbs <- function(region, beta, interaction_radius, phi, n_iter, d = 3, covariate_field = NULL, covariate_coeff = 0, buffer_size = NULL) {
  
  # Add a buffer region to account for edge effects
  buffer_region <- lapply(region, function(bounds) c(bounds[1] - buffer_size, bounds[2] + buffer_size))
  region_volume <- prod(sapply(buffer_region, function(bounds) diff(bounds))) # Compute the volume of the extended region
  
  # Initialize the point configuration
  points <- matrix(nrow = 0, ncol = d)
  
  # Function to compute lambda(u, x)
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
  
  # Function to compute r(u, x) for birth and death proposals
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
  
  # Initialize accepted points
  accepted_points <- list()
  
  for (i in 1:n_iter) {
    proposal_type <- sample(c("add", "remove"), 1, prob = c(0.5, 0.5))
    acceptance_prob <- 0
    
    if (proposal_type == "add") {
      new_point <- sapply(buffer_region, function(bounds) runif(1, bounds[1], bounds[2]))
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
  
  # Filter points to retain only those within the original region
  final_points <- points[apply(points, 1, function(pt) {
    all(sapply(1:d, function(dim) region[[dim]][1] <= pt[dim] && pt[dim] <= region[[dim]][2]))
  }), ]
  
  # Return the final configuration and all accepted configurations
  list(final_points = final_points, all_points = accepted_points)
}

# Define a 3D region
region_3d <- list(c(0, 10), c(0, 10), c(0, 10))

# Define a 3D covariate field (optional)
covariate_field_3d <- function(point) {
  x <- point[1]
  y <- point[2]
  z <- point[3]
  return(2 * x + y - z) 
}

# Define interaction function (repulsion)
phi <- function(r, interaction_radius) {
  ifelse(r < interaction_radius, -1 + (r / interaction_radius), 0)  # Rescaled for better estimation
}


# Parameters
beta_true <- 10  # True intensity
phi_coeff_true <- -1  # Repulsion strength
interaction_radius <- 1  # Influence range
n_iter <- 5000  # Number of iterations
covariate_coeff <- 2.0  # Covariate influence
buffer_size <- 2  # Extra space to reduce edge effects

# Run the Gibbs process simulation
result_3d <- simulate_gibbs(
  region_3d, beta_true, interaction_radius, phi, n_iter,
  d = 3, covariate_field = covariate_field_3d, covariate_coeff = covariate_coeff,
  buffer_size = buffer_size
)

# Prepare data list for Stan
stan_data <- list(
  d = 3,  # 3D space
  N = nrow(all_points_3d),  # Number of simulated points
  sampled_points = all_points_3d,  # Points as input
  interaction_radius = 1  # Same as used in simulation
)

# Run the Stan model
fit <- stan(
  file = "~/Documents/Skole/H2024/master/code/Point_Processes/inference_gibbs.stan",  # Stan file
  data = stan_data,
  iter = 5000,
  chains = 4
)

posterior <- extract(fit)
beta_est <- posterior$beta
phi_est <- posterior$phi_coeff
