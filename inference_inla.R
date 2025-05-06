# Required Libraries
library(INLA)
library(sp)
library(rgl)
library(fmesher)

# General LGCP Inference Framework
lgcp_inference <- function(d, sampled_points, bounds, m, covariate_field, seed = 123) {
  set.seed(seed)
  
  ### 1. Create the Mesh
  # Define non-convex hull boundary
  boundary <- fm_nonconvex_hull(as.matrix(sampled_points), convex = -0.1)
  
  # Create the mesh using fmesher
  mesh <- fm_mesh_2d_inla(
    loc = sampled_points,   # Input points
    boundary = boundary,    # Boundary of the domain
    max.edge = c(1, 2),     # Maximum edge length of triangles
    cutoff = 0.5            # Minimum distance between mesh nodes
  )
  
  # Visualize the Mesh
  if (d == 3) {
    library(rgl)
    plot3d(mesh$loc, col = "blue", size = 2, main = "Mesh in 3D")
    points3d(sampled_points[, 1], sampled_points[, 2], sampled_points[, 3], col = "red", size = 3)
  }
  
  ### 2. Define Covariate Field
  covariate_values <- apply(sampled_points, 1, covariate_field)
  
  ### 3. Define SPDE Model
  spde <- inla.spde2.pcmatern(
    mesh = mesh,
    alpha = 2,                 # Smoothness parameter
    prior.range = c(1, 0.5),   # Prior for range (P(range > 1) = 0.5)
    prior.sigma = c(1, 0.5)    # Prior for variance (P(sigma > 1) = 0.5)
  )
  
  ### 4. Aggregate Points into Grid (Counts)
  grid_sides <- lapply(bounds, function(lim) seq(lim[1], lim[2], length.out = m))
  grid_coords <- expand.grid(grid_sides)
  
  # Count points in each grid cell
  counts <- numeric(nrow(grid_coords))
  grid_step <- sapply(bounds, function(lim) diff(lim) / m)
  
  for (i in seq_along(counts)) {
    lower_bounds <- grid_coords[i, ] - grid_step / 2
    upper_bounds <- grid_coords[i, ] + grid_step / 2
    counts[i] <- sum(apply(sampled_points, 1, function(pt) all(pt >= lower_bounds & pt < upper_bounds)))
  }
  
  ### 5. Define INLA Stack
  A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(grid_coords))
  
  stack <- inla.stack(
    data = list(y = counts),  # Observed counts
    A = list(A, 1),          # Projector matrix and covariate effects
    effects = list(
      spatial = 1:spde$n.spde,  # Spatial effects (SPDE random field)
      covariate = covariate_values
    )
  )
  
  ### 6. Model Formulation in INLA
  formula <- y ~ 1 + covariate + f(spatial, model = spde)
  
  ### 7. Run INLA for Inference
  result <- inla(
    formula,
    family = "poisson",  # Poisson likelihood for counts
    data = inla.stack.data(stack),
    control.predictor = list(A = inla.stack.A(stack)),
    control.inla = list(int.strategy = "eb"),  # Use Laplace approximation
    control.compute = list(dic = TRUE, waic = TRUE)  # Compute model fit metrics
  )
  
  ### 8. Output Results
  list(
    result = result,  # INLA result object
    mesh = mesh,      # Mesh for visualization
    stack = stack,    # INLA stack
    grid_coords = grid_coords,  # Grid coordinates
    counts = counts   # Observed counts
  )
}


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


# Simulate 3D LGCP Data
data <- simulate_lgcp(
  d = 3,  # 3 Dimensions
  bounds = list(c(0, 5), c(0, 5), c(0, 5)),  # 3D bounds
  m = 5,  # Grid resolution (number of grid cells per dimension)
  length_scale = 1,  # Length scale for Gaussian random field
  covariate_field = function(x) { x[1] + x[2] + 0.5 * x[3] },  # Covariate function
  covariate_coeff = 0.5,  # Effect size of the covariate
  seed = 123
)

# Access simulated data
sampled_points <- data$sampled_points  # Observed points in 3D space
grid_coords <- data$grid_coords        # Grid coordinates
rates <- data$rates                    # Intensity values at grid points
sampled_points <- sampled_points[1:500, ]  # First 500 points

# Perform Inference on Simulated 3D Data
system.time({
  results <- lgcp_inference(
    d = 3,
    sampled_points = data$sampled_points,
    bounds = list(c(0, 5), c(0, 5), c(0, 5)),
    m = 5,
    covariate_field = function(x) { x[1] + x[2] + 0.5 * x[3] }
  )
})


