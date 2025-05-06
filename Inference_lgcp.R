# Load necessary libraries
library(fmesher)
library(geometry)
library(ggplot2)
library(INLA)
library(pracma)
library(rstan)
library(FNN)
library(bayesplot)
# --- STEP 1: SIMULATE LGCP DATA ---

simulate_lgcp <- function(d, bounds, m, length_scale, covariate_field, covariate_coeff, seed = 123) {
  
  set.seed(seed)
  xlims <- bounds  # List of bounds for each dimension
  grid_sides <- lapply(xlims, function(lim) seq(lim[1], lim[2], length.out = m)) # Calculates the size of each grid cell along each dimension
  
  # Create grid coordinates
  grid_coords <- expand.grid(grid_sides)
  X <- as.matrix(grid_coords)
  grid_step <- sapply(xlims, function(lim) diff(lim) / m)  # Size of each grid cell
  
  # Powered Exponential Covariance function
  rbf <- function(X1, X2, length_scale) {
    dists <- as.matrix(dist(rbind(X1, X2)))
    exp(-dists[1:nrow(X1), (nrow(X1)+1):(nrow(X1)+nrow(X2))] ^ 2 / (2 * length_scale ^ 2))
  }
  
  # Draw sample from Gaussian Process
  K <- rbf(X, X, length_scale)
  Y <- MASS::mvrnorm(mu = rep(0, nrow(X)), Sigma = K) #Draws samples from a multivariate normal distribution
  # Y represents the Gaussian Random Field
  
  
  # Covariate effect: apply covariate field
  covariate_values <- apply(X, 1, covariate_field)
  covariate_effect <- exp(covariate_coeff * covariate_values) 
  
  # Incorporate covariates
  Z <- exp(Y) * covariate_effect  # Modify the rate with the covariate effect
  
  # Generate points proportional to GRF rates
  volume <- prod(sapply(xlims, diff))  # Volume of the d-dimensional space
  total_rate <- sum(Z) * (volume / (m^d))  # Integral approximation over grid
  N <- rpois(1, total_rate)  # Total number of points
  
  # Sample grid cells proportional to GRF rates
  probabilities <- Z / sum(Z)  # Normalize rates to probabilities
  sampled_indices <- sample(1:length(Z), size = N, replace = TRUE, prob = probabilities)
  
  # Generate random points within the sampled grid cells
  sampled_points <- matrix(0, nrow = N, ncol = d)
  for (i in seq_len(d)) {
    cell_centers <- grid_coords[sampled_indices, i]
    offsets <- runif(N, min = -grid_step[i] / 2, max = grid_step[i] / 2)  # Random offset within cell
    sampled_points[, i] <- cell_centers + offsets
  }
  
  # Return all points and rates
  list(
    sampled_points = sampled_points,  # Points sampled from the LGCP
    rates = Z,                        # GRF rates
    grid_coords = X                   # Grid coordinates
  )
}

# Define the 2D domain and parameters
d <- 2  
bounds <- list(c(0, 10), c(0, 10))  # 2D space
m <- 20  # Resolution
length_scale <- 1.0  
covariate_coeff <- 0.5  

# Define Covariate Function (to be replaced with SPDE-based field later)
covariate_field <- function(x) sin(x[1] / 2) + cos(x[2] / 3)  

# Simulate LGCP in 2D
lgcp_result <- simulate_lgcp(d, bounds, m, length_scale, covariate_field, covariate_coeff, seed = 123)

# Extract results
lgcp_points <- lgcp_result$sampled_points
rates <- lgcp_result$rates
X <- lgcp_result$grid_coords


mesh <- inla.mesh.2d(
  loc = rbind(lgcp_points, X),  # Use observed points and grid
  max.edge = c(1, 2),  # Increase to coarsen the mesh
  cutoff = 0.2  # Avoid tiny triangles
  )

# Plot mesh
plot(mesh)
points(lgcp_points, col = "red", pch = 20)


nu <- 1           # Smoothness parameter (commonly 1 for SPDE)
sigma <- 1.0      # Marginal standard deviation of the GRF

# Compute kappa
kappa <- sqrt(8 * nu) / length_scale

# Compute tau using the correct formula
tau <- sigma * kappa^nu * sqrt(gamma(nu + d/2) * (4 * pi)^(d/2) / gamma(nu))


# Define SPDE model with computed parameters
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
theta <- c(log(tau), log(kappa))  # SPDE hyperparameters
Q <- inla.spde.precision(spde, theta = theta)  # Compute precision matrix

# Convert precision matrix to standard format for Stan
Q_matrix <- as.matrix(Q)
write.table(Q_matrix, "Q_matrix.txt", row.names = FALSE, col.names = FALSE)


lgcp_points <- as.matrix(lgcp_points)
mesh_index <- apply(lgcp_points, 1, function(pt) {
  which.min(rowSums((mesh$loc - matrix(pt, nrow = nrow(mesh$loc), ncol = ncol(mesh$loc), byrow = TRUE))^2))
})



# Compute prior means from SPDE parameters
log_tau_prior_mean <- log(tau)
log_kappa_prior_mean <- log(kappa)

# Prepare Stan data
stan_data <- list(
  N = nrow(Q_matrix),  # Number of mesh nodes
  M = nrow(lgcp_points),  # Number of observed points
  Q = Q_matrix,  # Precision matrix
  Y = rep(1, nrow(lgcp_points)),  # Poisson count data
  mesh_index = mesh_index,  # Mapped indices
  log_tau_prior_mean = log_tau_prior_mean,
  log_kappa_prior_mean = log_kappa_prior_mean
)

# Run Stan
fit <- stan(
  file = "~/Documents/Skole/H2024/master/code/Point_Processes/lgcp_model.stan",
  data = stan_data,
  iter = 500, cores = 4
)

log_tau_post <- extract(fit, pars = "log_tau")$log_tau
log_kappa_post <- extract(fit, pars = "log_kappa")$log_kappa

mcmc_rhat(rhat(fit))
mcmc_neff(neff_ratio(fit))
mcmc_trace(fit, pars = c("log_tau", "log_kappa"))
