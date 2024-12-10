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

