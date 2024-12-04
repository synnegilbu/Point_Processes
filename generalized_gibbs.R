simulate_gibbs_with_covariates <- function(region, beta, gamma, r, n_iter, burn_in, 
                                           covariate_field, covariate_coeff, points = NULL, d = 2) {
  
  # Function to calculate pairwise interaction energy
  pairwise_interaction <- function(points, r) {
    if (nrow(points) < 2) return(0)  
    dist_matrix <- as.matrix(dist(points))
    sum(dist_matrix[upper.tri(dist_matrix)] < r)  # Count pairs within radius r
  }
  
  # Compute the covariate effect for the current configuration of points
  compute_covariate_effect <- function(points, covariate_field, covariate_coeff) {
    if (nrow(points) == 0) return(0)
    cov_values <- apply(points, 1, covariate_field)  
    sum(covariate_coeff * cov_values)  
  }
  
  # Initialize points matrix if not provided
  if (is.null(points)) {
    points <- matrix(nrow = 0, ncol = d)
  }
  
  # Initialize accepted configurations
  accepted_points <- list()
  
  for (i in 1:n_iter) {
    # Initialize acceptance probability
    acceptance_prob <- 0
    
    # Proposal: Add or Remove a point
    proposal_type <- sample(c("add", "remove"), 1, prob = c(0.5, 0.5))
    
    if (proposal_type == "add") {
      new_point <- sapply(region, function(bounds) runif(1, bounds[1], bounds[2]))
      proposed_points <- rbind(points, new_point)
      
      # Compute energies and covariate effects
      current_energy <- pairwise_interaction(points, r)
      proposed_energy <- pairwise_interaction(proposed_points, r)
      current_covariate <- compute_covariate_effect(points, covariate_field, covariate_coeff)
      proposed_covariate <- compute_covariate_effect(proposed_points, covariate_field, covariate_coeff)
      
      acceptance_prob <- min(1, beta * gamma^proposed_energy * exp(proposed_covariate - current_covariate))
      
    } else if (proposal_type == "remove" && nrow(points) > 0) {
      # Remove a randomly chosen point
      remove_idx <- sample(1:nrow(points), 1)
      proposed_points <- points[-remove_idx, , drop = FALSE]
      
      # Compute energies and covariate effects
      current_energy <- pairwise_interaction(points, r)
      proposed_energy <- pairwise_interaction(proposed_points, r)
      current_covariate <- compute_covariate_effect(points, covariate_field, covariate_coeff)
      proposed_covariate <- compute_covariate_effect(proposed_points, covariate_field, covariate_coeff)
      
      acceptance_prob <- min(1, nrow(points) / (beta * gamma^current_energy * exp(current_covariate - proposed_covariate)))
    } else {
      proposed_points <- points  # If no points exist to remove
    }
    
    # Accept or reject the proposal based on the acceptance probability
    if (runif(1) < acceptance_prob) {
      points <- proposed_points
    }
    
    # Save accepted configuration after burn-in
    if (i > burn_in) {
      accepted_points[[i - burn_in]] <- points
    }
  }
  
  # Return the final configuration and all saved configurations
  list(final_points = points, all_points = accepted_points)
}


set.seed(123)  # For reproducibility

# Define the 2D region
region <- list(c(0, 10), c(0, 10))  # Region is [0, 10] x [0, 10]

# Define a covariate field
covariate_field <- function(point) {
  # A simple example: linear gradient covariate
  x <- point[1]
  y <- point[2]
  return(x + y)
}

# Parameters for the Gibbs process
beta <- 60          # Baseline intensity
gamma <- 0.8       # Interaction strength
r <- 1.5           # Interaction radius
n_iter <- 1000     # Number of iterations
burn_in <- 500     # Burn-in period
covariate_coeff <- 0.1  # Weight of the covariate effect

# Run the simulation
result <- simulate_gibbs_with_covariates(region, beta, gamma, r, n_iter, burn_in, 
                                         covariate_field, covariate_coeff, d = 2)

# Extract final points
final_points <- result$final_points

# Plot the points
plot(final_points[, 1], final_points[, 2], 
     pch = 16, col = "blue", xlab = "X", ylab = "Y", 
     main = "Simulated Gibbs Process (2D)",
     xlim = region[[1]], ylim = region[[2]])



set.seed(42)

# Define 3D region
region_3d <- list(c(0, 1), c(0, 1), c(0, 1))

# Covariate function: simple linear combination of coordinates
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(x + y + z)  # Simple covariate field
}

# Simulation parameters
beta <- 150                  # High intensity
gamma <- 0.2                 # Weak repulsion (gamma < 1)
r <- 0.05                     # Small interaction radius
n_iter <- 5000               # Moderate iterations
burn_in <- 500               # Burn-in iterations
covariate_coeff <- 0.5       # Moderate covariate effect

# Run simulation
result_1 <- simulate_gibbs_with_covariates(region = region_3d, beta = beta, gamma = gamma, r = r, 
                                           n_iter = n_iter, burn_in = burn_in, 
                                           covariate_field = covariate_field, covariate_coeff = covariate_coeff,
                                           points = matrix(nrow = 0, ncol = 3), d = 3)

# Plot the final configuration
final_points_1 <- result_1$final_points
library(rgl)
plot3d(final_points_1[, 1], final_points_1[, 2], final_points_1[, 3], col = "red", size = 5,
       main = "High Intensity, Weak Repulsion")




set.seed(42)
# Define 3D region
region_3d <- list(c(0, 1), c(0, 1), c(0, 1))

# Covariate function: a simple spatial gradient
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(2 * x + y)  # Linear gradient covariate
}

# Simulation parameters
beta <- 50                  # Low intensity
gamma <- 0.3                 # Strong repulsion (gamma < 1)
r <- 0.1                     # Larger interaction radius
n_iter <- 5000               # Moderate iterations
burn_in <- 500               # Burn-in iterations
covariate_coeff <- 0.2       # Low covariate effect

# Run simulation
result_2 <- simulate_gibbs_with_covariates(region = region_3d, beta = beta, gamma = gamma, r = r, 
                                           n_iter = n_iter, burn_in = burn_in, 
                                           covariate_field = covariate_field, covariate_coeff = covariate_coeff,
                                           points = matrix(nrow = 0, ncol = 3), d = 3)

# Plot the final configuration
final_points_2 <- result_2$final_points
plot3d(final_points_2[, 1], final_points_2[, 2], final_points_2[, 3], col = "green", size = 5,
       main = "Low Intensity, Strong Repulsion")

set.seed(42)

# Define 3D region
region_3d <- list(c(0, 1), c(0, 1), c(0, 1))

# Covariate function: quadratic covariate for added complexity
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(x^2 + y^2 + z^2)  # Quadratic covariate
}

# Simulation parameters
beta <- 100                  # High intensity
gamma <- 1.5                 # Attraction (gamma > 1)
r <- 0.05                    # Small interaction radius
n_iter <- 5000               # Moderate iterations
burn_in <- 500               # Burn-in iterations
covariate_coeff <- 0.3       # Moderate covariate effect

# Run simulation
result_3 <- simulate_gibbs_with_covariates(region = region_3d, beta = beta, gamma = gamma, r = r, 
                                           n_iter = n_iter, burn_in = burn_in, 
                                           covariate_field = covariate_field, covariate_coeff = covariate_coeff,
                                           points = matrix(nrow = 0, ncol = 3), d = 3)

# Plot the final configuration
final_points_3 <- result_3$final_points
plot3d(final_points_3[, 1], final_points_3[, 2], final_points_3[, 3], col = "blue", size = 5,
       main = "High Intensity, Attraction (gamma > 1)")


set.seed(42)

# Define 3D region
region_3d <- list(c(0, 1), c(0, 1), c(0, 1))

# Covariate function: quadratic covariate to enhance spatial variability
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(2 * x^2 + 3 * y^2 + z)  # Complex covariate
}

# Simulation parameters
beta <- 150                  # High intensity
gamma <- 2.0                 # Strong attraction (gamma > 1)
r <- 0.1                     # Moderate interaction radius
n_iter <- 5000               # Moderate iterations
burn_in <- 500               # Burn-in iterations
covariate_coeff <- 1.0       # Strong covariate effect

# Run simulation
result_4 <- simulate_gibbs_with_covariates(region = region_3d, beta = beta, gamma = gamma, r = r, 
                                           n_iter = n_iter, burn_in = burn_in, 
                                           covariate_field = covariate_field, covariate_coeff = covariate_coeff,
                                           points = matrix(nrow = 0, ncol = 3), d = 3)

# Plot the final configuration
final_points_4 <- result_4$final_points
plot3d(final_points_4[, 1], final_points_4[, 2], final_points_4[, 3], col = "purple", size = 5,
       main = "High Intensity, Strong Attraction, Strong Covariate Effect")
