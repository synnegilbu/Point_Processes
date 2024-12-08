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

# Define a covariate field (simple example: linear gradient in 3D)
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(2 * x + 3 * y - z)  # Example: a linear combination of coordinates
}

# Set parameters for the Poisson process
d <- 3  # Dimensionality
bounds <- list(c(0, 10), c(0, 10), c(0, 10))  # Bounds for the 3D region
m <- 10  # Grid resolution
length_scale <- 2  # Length scale for the Gaussian process kernel
covariate_coeff <- 0.5  # Covariate coefficient (controls the effect of the covariate)

# Generate the Poisson point process with covariates
result <- generate_poisson_gp_with_covariates(
  d = d,
  bounds = bounds,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Extract the accepted points
accepted_points <- result$accepted_points

# 3D Plot
plot_ly(
  x = accepted_points[, 1],
  y = accepted_points[, 2],
  z = accepted_points[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "blue", opacity = 0.6)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z")
    ),
    title = "Non-Homogeneous Poisson Process with Covariates in 3D"
  )

set.seed(42)

# Define the 3D region
region_3d <- list(c(0, 10), c(0, 10), c(0, 10))

# Covariate function: exponential decay
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(exp(-(x^2 + y^2 + z^2) / 10))  # Exponential decay based on distance from origin
}

# Simulation parameters
beta <- 100
gamma <- 0.5
r <- 0.1
m <- 10  # Grid resolution
length_scale <- 1
covariate_coeff <- 2  # Strong covariate effect

# Run simulation
result_2 <- generate_poisson_gp_with_covariates(
  d = 3,
  bounds = region_3d,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Plot the result
accepted_points_2 <- result_2$accepted_points
plot_ly(
  x = accepted_points_2[, 1],
  y = accepted_points_2[, 2],
  z = accepted_points_2[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "red", opacity = 0.6)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z")
    ),
    title = "Exponential Covariate Field with High Covariate Coefficient"
  )
set.seed(42)

# Define the 3D region
region_3d <- list(c(0, 10), c(0, 10), c(0, 10))

# Covariate function: sine-wave oscillation
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]; z <- u[3]
  return(sin(x) * cos(y) * sin(z))  # Sine-wave oscillation in 3D space
}

# Simulation parameters
beta <- 100
gamma <- 0.5
r <- 0.1
m <- 10  # Grid resolution
length_scale <- 1
covariate_coeff <- 0.5  # Moderate covariate effect

# Run simulation
result_3 <- generate_poisson_gp_with_covariates(
  d = 3,
  bounds = region_3d,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Plot the result
accepted_points_3 <- result_3$accepted_points
plot_ly(
  x = accepted_points_3[, 1],
  y = accepted_points_3[, 2],
  z = accepted_points_3[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "green", opacity = 0.6)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z")
    ),
    title = "Sine-Wave Covariate Field"
  )


set.seed(42)

# Define the 3D region
region_3d <- list(c(0, 10), c(0, 10), c(0, 10))

# Covariate function: uniform covariate field
covariate_field <- function(u) {
  return(1)  # Constant value across all space (uniform covariate field)
}

# Simulation parameters
beta <- 100
gamma <- 0.5
r <- 0.1
m <- 10  # Grid resolution
length_scale <- 1
covariate_coeff <- 0 # Strong covariate effect

# Run simulation
result_4 <- generate_poisson_gp_with_covariates(
  d = 3,
  bounds = region_3d,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Plot the result
accepted_points_4 <- result_4$accepted_points
plot_ly(
  x = accepted_points_4[, 1],
  y = accepted_points_4[, 2],
  z = accepted_points_4[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 3, color = "purple", opacity = 0.6)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z")
    ),
    title = "Uniform Covariate Field with Strong Effect"
  )

# Define a covariate field (simple example: linear gradient in 2D)
covariate_field <- function(u) {
  x <- u[1]; y <- u[2]
  return(2 * x + 3 * y)  # Example: linear covariate, sum of coordinates
}

# Set parameters for the Poisson process
d <- 2  # Dimensionality (2D)
bounds <- list(c(0, 10), c(0, 10))  # Bounds for the 2D region
m <- 60  # Grid resolution (density of points in the grid)
length_scale <- 2  # Length scale for the Gaussian process kernel
covariate_coeff <- 0.009  # Moderate covariate effect

# Generate the Poisson point process with covariates
result_2d <- generate_poisson_gp_with_covariates(
  d = d,
  bounds = bounds,
  m = m,
  length_scale = length_scale,
  covariate_field = covariate_field,
  covariate_coeff = covariate_coeff,
  seed = 42
)

# Check if the result is empty
if (length(result_2d$accepted_points) == 0) {
  stop("No points were accepted. Adjust the rate and volume settings.")
}

# Extract the accepted points, grid coordinates, and rates
accepted_points_2d <- result_2d$accepted_points
grid_coords_2d <- result_2d$grid_coords
rates_2d <- result_2d$rates

# Debugging: Check the range of accepted points and grid
cat("Range of accepted points:\n")
print(range(accepted_points_2d))

cat("\nRange of grid coordinates:\n")
print(range(grid_coords_2d))

# Debugging: Check the rate (Z) and zmax
cat("\nMaximum rate (zmax):", max(rates_2d), "\n")
cat("Rate values (Z):", summary(rates_2d), "\n")

# Convert grid coordinates and rates into a data frame for plotting
grid_df <- data.frame(x = grid_coords_2d[, 1], y = grid_coords_2d[, 2], rate = rates_2d)

# Convert accepted points into a data frame for plotting
accepted_points_df <- data.frame(x = accepted_points_2d[, 1], y = accepted_points_2d[, 2])

# 1. Plot the Gaussian Random Field (GRF) as a heatmap
ggplot(grid_df, aes(x = x, y = y, fill = rate)) +
  geom_tile() +  # Creates the heatmap for the GRF
  scale_fill_viridis() +  # Color scale for the GRF
  labs(title = "Gaussian Random Field (GRF) with Covariate Effect", x = "X", y = "Y") +
  theme_minimal()

# 2. 2D Plot of the Poisson points on top of the GRF heatmap
ggplot() +
  geom_tile() +  # Creates the heatmap for the GRF
  scale_fill_viridis() +  # Color scale for the GRF
  geom_point(data = accepted_points_df, aes(x = x, y = y), color = "red", size = 2) +  # Poisson points in red
  labs(x = "X", y = "Y") +
  theme_minimal()
