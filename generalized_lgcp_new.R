# Load necessary libraries
library(ggplot2)
library(rgl)
# Parameters:
#   d: dimensions
#   bounds: space bounds
#   m: grid resolution
#   length_scale: correlation length scale
#   covariate_field: user-defined function
#   covariate_coeff: covariate coefficient
#   seed: random seed

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


#2D simulation with GRF underneath
bounds <- list(c(0, 10), c(0, 10))  
m <- 60                           
length_scale <- 1                  
covariate_field <- function(point) sin(point[1] / 2) + cos(point[2] / 3)
covariate_coeff <- 0.5           

# Simulate LGCP
result <- simulate_lgcp(
  d = 2,
  bounds = bounds,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Extract results
sampled_points <- result$sampled_points
rates <- result$rates
grid_coords <- result$grid_coords

grf_data <- data.frame(
  x = grid_coords[, 1],
  y = grid_coords[, 2],
  z = rates
)

# Visualize Sampled Points + GRF
points_data <- data.frame(
  x = sampled_points[, 1],
  y = sampled_points[, 2]
)

ggplot() +
  geom_tile(data = grf_data, aes(x = x, y = y, fill = z)) +
  scale_fill_viridis_c() +
  geom_point(data = points_data, aes(x = x, y = y), color = "red", size = 0.5) +
  labs(title = "Gaussian Random Field with Thinned Poisson Points",
       x = "X", y = "Y", fill = "Rate") +
  theme_minimal()


# Define parameters for 3D simulation
result <- simulate_lgcp(
  d = 3,
  bounds = list(c(0, 10), c(0, 10), c(0, 10)),
  m = 10,
  length_scale = 2,
  covariate_field = function(x) { sin(x[1]) + cos(x[2]) + 0.5 * x[3] },
  covariate_coeff = 0.5,
  seed = 123
)


# Extract results
sampled_points <- result$sampled_points
rates <- result$rates
grid_coords <- result$grid_coords


# Plot sampled points
plot3d(sampled_points[, 1], sampled_points[, 2], sampled_points[, 3],
       col = "red", type = "p",
       xlab = "X", ylab = "Y", zlab = "Z", main = "3D LGCP Simulated Points")

# To check if this gives a valid model, you could run the K-function in the
# spatstat-package. This will show that we have clustering
