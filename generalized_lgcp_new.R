library(ggplot2)
library(reshape2)
library(tidyverse)
library(plotly)

generate_poisson_gp <- function(d, bounds, m, length_scale, seed = 42) {
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
  
  # Convert to rates
  Z <- exp(Y)
  zmax <- max(Z)
  
  # Draw Poisson
  volume <- prod(sapply(xlims, diff))  # Volume of the d-dimensional space
  N <- rpois(1, zmax * volume)
  
  # Generate random points in the d-dimensional space
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
  list(accepted_points = points[kp, , drop = FALSE], rates = Z)
}


result <- generate_poisson_gp(
  d = 3,
  bounds = list(c(0, 10), c(0, 10), c(0, 10)),
  m = 10,  
  length_scale = 2,
  seed = 42
)

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
    title = "Non-Homogeneous Poisson Process in 3D"
  )

