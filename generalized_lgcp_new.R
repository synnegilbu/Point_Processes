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


# Load necessary packages
library(spatstat)

# Define parameters for the LGCP simulation
bounds <- list(c(0, 1), c(0, 1))  # 2D space bounds (x and y)
m <- 50  # Grid size
length_scale <- 0.3  # Gaussian process length scale
covariate_field <- function(x) x[1]  # Covariate: x-coordinate
covariate_coeff <- 0.5  # Covariate effect

# Simulate LGCP points using your function
sim_result <- simulate_lgcp(d = 2, bounds = bounds, m = m,
                            length_scale = length_scale,
                            covariate_field = covariate_field,
                            covariate_coeff = covariate_coeff)

# Convert sampled points to a `ppp` object (spatstat point pattern format)
sampled_points <- sim_result$sampled_points
ppp_obj <- ppp(x = sampled_points[, 1], y = sampled_points[, 2],
               window = owin(c(0, 1), c(0, 1)))

# Perform likelihood-based inference using kppm (for LGCP)
lgcp_fit <- kppm(ppp_obj ~ 1, clusters = "LGCP")

# View results
summary(lgcp_fit)
plot(lgcp_fit)  # Plot fitted intensity


# INFERENCE USING INLA
# Install necessary packages
if (!require("INLA")) install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")

library(INLA)
library(sp)
library(MASS)  # For multivariate normal sampling

# Increase covariate effect and grid resolution
bounds <- list(c(0, 1), c(0, 1))  # Domain
m <- 100  # Finer grid
covariate_coeff <- 2  # Stronger covariate effect

# Simulate LGCP
sim_result <- simulate_lgcp(d = 2, bounds = bounds, m = m,
                            length_scale = 0.1,  # Shorter length scale for more variation
                            covariate_field = function(x) 5 * x[1] + 2 * x[2],
                            covariate_coeff = covariate_coeff)

# Check if points are generated
if (is.null(nrow(sim_result$sampled_points)) || nrow(sim_result$sampled_points) == 0) {
  stop("No points generated. Try increasing covariate effect or grid resolution.")
} else {
  print(paste("Generated", nrow(sim_result$sampled_points), "points."))
}

sampled_points <- sim_result$sampled_points  # Extract sampled points
points_data <- data.frame(
  x = sampled_points[, 1],
  y = sampled_points[, 2]
)

# Step 2: Convert sampled points to SpatialPoints object
coords <- as.data.frame(sampled_points)
colnames(coords) <- c("x", "y")
coords$x <- as.numeric(coords$x)
coords$y <- as.numeric(coords$y)
coordinates(coords) <- ~ x + y

# Step 3: Create spatial mesh for INLA
mesh <- inla.mesh.2d(loc = coords@coords, 
                     max.edge = c(0.05, 0.1),  # Smaller values for finer triangles
                     cutoff = 0.001)  # Lower value to include more vertices close to each other
plot(mesh, main = "Adjusted 2D Spatial Mesh")


# Step 4: Define SPDE (Stochastic Partial Differential Equation) model for Gaussian process
spde <- inla.spde2.pcmatern(mesh = mesh, alpha = 2, prior.range = c(0.1, 0.5), prior.sigma = c(1, 0.01))

# Step 5: Prepare data for INLA
n_points <- nrow(sampled_points)
response <- rep(1, n_points)  # Dummy response for Poisson point pattern
A <- inla.spde.make.A(mesh, loc = coords@coords)  # Projector matrix

# Stack the data (INLA format)
stack <- inla.stack(data = list(y = response),
                    A = list(A, 1),
                    effects = list(list(spatial_field = 1:mesh$n),
                                   data.frame(intercept = 1)))

# Step 6: Fit the LGCP model using INLA
formula <- y ~ 0 + intercept + f(spatial_field, model = spde)

result <- inla(formula,
               family = "poisson",
               data = inla.stack.data(stack),
               control.predictor = list(A = inla.stack.A(stack)),
               control.inla = list(strategy = "adaptive"))

# Step 7: View posterior results
summary(result)

# Step 8: Plot posterior mean intensity
inla.plot.mesh(mesh)
inla.plot.field(result$summary.random$spatial_field$mean, mesh)
