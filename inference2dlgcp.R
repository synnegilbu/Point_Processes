library(geometry)
library(rgl)
library(Matrix)
library(FNN) 
library(rstan)


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


result <- simulate_lgcp(
  d = 3,
  bounds = list(c(0, 10), c(0, 10), c(0, 10)),
  m = 10,
  length_scale = 2,
  covariate_field = function(x) { x[1]+ x[2] + 0.5 * x[3] },
  covariate_coeff = 0.05,
  seed = 123
)


# Extract results
sampled_points <- result$sampled_points
rates <- result$rates
grid_coords <- result$grid_coords


compute_fem_matrices <- function(mesh, points, d) {
  num_nodes <- max(mesh)  # Total number of mesh nodes
  
  # Initialize sparse matrices for mass (C) and stiffness (G)
  C <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(num_nodes, num_nodes))
  G <- sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(num_nodes, num_nodes))
  
  for (i in 1:nrow(mesh)) {
    simplex_nodes <- mesh[i, ]
    coords <- points[simplex_nodes, ]
    
    # Compute volume of the simplex 
    B <- cbind(1, coords)
    volume <- abs(det(B)) / factorial(d)
    
    # Compute gradient of basis functions
    B_inv <- solve(B)
    grad_phi <- t(B_inv[-1, ])  
    
    # Initialize local matrices
    local_G <- matrix(0, d + 1, d + 1)
    local_C <- matrix(0, d + 1, d + 1)
    
    local_G <- volume * (grad_phi %*% t(grad_phi))
    local_C <- (volume / (d + 1)) * matrix(1, d + 1, d + 1)
    
    # Assemble into global matrices
    for (j in 1:(d + 1)) {
      for (k in 1:(d + 1)) {
        C[simplex_nodes[j], simplex_nodes[k]] <- C[simplex_nodes[j], simplex_nodes[k]] + local_C[j, k]
        G[simplex_nodes[j], simplex_nodes[k]] <- G[simplex_nodes[j], simplex_nodes[k]] + local_G[j, k]
      }
    }
  }
  
  return(list(C = C, G = G))
}

compute_precision_matrix <- function(C, G, kappa) {
  C_inv <- solve(C)  # Compute inverse of C
  Q <- C_inv %*% (kappa^2 * C + G) %*% C_inv  # Compute precision matrix
  return(as(Q, "sparseMatrix"))  # Convert to sparse format
}


# Generate 3D sample points
set.seed(123)



# Generate Delaunay triangulation mesh
mesh <- delaunayn(sampled_points)
num_points <- max(mesh)

# Compute FEM matrices
fem_matrices <- compute_fem_matrices(mesh, sampled_points, d = 3)
C <- fem_matrices$C
G <- fem_matrices$G

# Define Matérn SPDE parameter
kappa <- sqrt(8*(7/2))/2  

# Compute the precision matrix
Q <- compute_precision_matrix(C, G, kappa)

# Convert Q to CSR format for Stan
Q_sparse <- as(Q, "dgCMatrix")

# Assuming `Q` is a sparse matrix
Q_dense <- as.matrix(Q_sparse)  # Convert sparse matrix Q to dense format
mesh_index <- knnx.index(sampled_points, grid_coords, k = 1)

# Generate synthetic observation data (Poisson-distributed counts)
Y <- rpois(length(mesh_index), lambda = 1)



mesh_index <- as.integer(c(mesh_index))  


# Convert precision matrix Q to sparse format (Compressed Sparse Row - CSR)
Q_sparse <- as(Q, "dgCMatrix")  

# Extract CSR components properly
Q_w <- Q_sparse@x              # Non-zero values
Q_v <- Q_sparse@i + 1          # Row indices (1-based for Stan)
Q_u <- Q_sparse@p + 1          # Column start indices (1-based for Stan)


stan_data <- list(
  N = nrow(Q_sparse),   # Number of mesh nodes
  M = length(mesh_index), 
  Q_w = Q_w,            # Non-zero values
  Q_v = Q_v,            # Row indices (1-based)
  Q_u = Q_u,            # Column start indices (1-based)
  Y = Y,                # Observed data
  mesh_index = mesh_index,  
  log_kappa_prior_mean = log(kappa),
  nnz_Q = length(Q_w)   # Number of non-zero elements
)



fit <- stan(
  file = "~/Documents/Skole/H2024/master/code/Point_Processes/try_inference.stan",
  data = stan_data,
  iter = 500, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)


