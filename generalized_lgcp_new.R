library(ggplot2)
library(reshape2)
library(tidyverse)
library(plotly)
library(MASS)
library(viridis)


generate_poisson_gp_with_covariates <- function(d, bounds, m, length_scale, covariate_field, covariate_coeff, seed = 42) {
  # Parameters
  set.seed(seed)
  xlims <- bounds  # List of bounds for each dimension: e.g., list(c(0, 10), c(0, 10))
  grid_sides <- lapply(xlims, function(lim) seq(lim[1], lim[2], length.out = m))
  
  # Create grid coordinates
  grid_coords <- expand.grid(grid_sides)
  X <- as.matrix(grid_coords)
  grid_step <- sapply(xlims, function(lim) diff(lim) / m)
  X <- sweep(X, 2, grid_step / 2, "+")  # Center on grid square
  
  # RBF Kernel Function
  rbf <- function(X1, X2, length_scale) {
    dists <- as.matrix(dist(rbind(X1, X2)))
    exp(-dists[1:nrow(X1), (nrow(X1)+1):(nrow(X1)+nrow(X2))] ^ 2 / (2 * length_scale ^ 2))
  }
  
  # Draw sample from GP
  K <- rbf(X, X, length_scale)
  Y <- MASS::mvrnorm(mu = rep(0, nrow(X)), Sigma = K)
  
  # Covariate effect: apply covariate field
  covariate_values <- apply(X, 1, covariate_field)  # Get covariate value at each point
  covariate_effect <- exp(covariate_coeff * covariate_values)  # Apply covariate coefficient
  
  # Convert to rates (incorporate covariate)
  Z <- exp(Y) * covariate_effect  # Modify the rate with the covariate effect
  zmax <- max(Z)
  
  # Draw Poisson (adjust N based on zmax and volume)
  volume <- prod(sapply(xlims, diff))  # Volume of the d-dimensional space
  N <- min(1e6, rpois(1, zmax * volume))  # Limit N to avoid excessive memory usage
  
  # Generate random points in the d-dimensional space, but limit the total points
  points <- matrix(runif(N * d), ncol = d)
  for (i in seq_len(d)) {
    points[, i] <- xlims[[i]][1] + points[, i] * diff(xlims[[i]])
  }
  
  # Find grid square indices for events
  indices <- lapply(seq_len(d), function(i) findInterval(points[, i], grid_sides[[i]]))
  flat_indices <- Reduce(function(a, b) (a - 1) * m + b, indices)
  
  # Perform thinning
  kp <- which(runif(N) <= (Z[flat_indices] / zmax))
  
  # Return accepted points and rates
  list(accepted_points = points[kp, , drop = FALSE], rates = Z, grid_coords = X)
}

