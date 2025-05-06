library(geometry)
library(Matrix)
library(MASS)
library(INLA)
library(FNN)

# General mesh in d dimensions
build_mesh <- function(d, m, bounds) {
  grid_axes <- lapply(bounds, function(lim) seq(lim[1], lim[2], length.out = m))
  grid_points <- expand.grid(grid_axes)
  pts <- as.matrix(grid_points)
  
  simplices <- delaunayn(pts, options = "QJ")  
  
  list(points = pts, simplices = simplices)
}

# FEM matrix builder
compute_fem_matrices <- function(points, simplices) {
  n <- nrow(points)
  d <- ncol(points)
  C <- Matrix(0, n, n, sparse = TRUE)
  G <- Matrix(0, n, n, sparse = TRUE)
  
  for (i in 1:nrow(simplices)) {
    idx <- simplices[i, ]
    verts <- points[idx, , drop = FALSE]
    
    # Build matrix B for the affine transformation
    B <- t(verts[-1, , drop = FALSE] - matrix(verts[1, ], d, d, byrow = TRUE))  # d x d
    detB <- det(B)
    
    # Skip degenerate simplex
    if (abs(detB) < 1e-12) {
      warning(paste("Skipping degenerate simplex at row", i, "with det =", detB))
      next
    }
    
    vol <- abs(detB) / factorial(d)
    
    # Gradients of basis functions in physical space
    Ghat <- cbind(-1, diag(d))  # Gradients on reference element
    grads <- solve(t(B), Ghat)  # Each column is grad(phi_i)
    
    
    # Stiffness matrix (G)
    GK <- matrix(0, d+1, d+1)
    for (j in 1:(d+1)) {
      for (k in 1:(d+1)) {
        GK[j, k] <- vol * sum(grads[, j] * grads[, k])
      }
    }
    
    # Mass matrix (C) using exact integration
    MK <- matrix(1, d+1, d+1)
    diag(MK) <- 2
    MK <- MK * vol / ((d + 1) * (d + 2))
    
    # Assemble global matrices
    C[idx, idx] <- C[idx, idx] + MK
    G[idx, idx] <- G[idx, idx] + GK
  }
  
  list(C = C, G = G)
}


# Precision matrix
assemble_precision_matrix <- function(C, G, tau = 1, kappa = 1, jitter = 1e-5) {
  Q <- tau^2 * (kappa^2 * C + G)
  Q + Diagonal(nrow(Q), jitter)
}

# Latent field
simulate_latent_field <- function(Q) {
  as.vector(inla.qsample(n = 1, Q = Q))
}

# Intensity function
build_intensity <- function(Y, coords, covariate_fn, beta, scale_intensity = 1000) {
  cov_vals <- apply(coords, 1, covariate_fn)
  eta <- Y + beta * cov_vals
  eta <- pmin(pmax(eta, -6), 6)  # clamp for stability
  lambda <- scale_intensity * exp(eta)
  list(lambda = lambda, covariate = cov_vals)
}

# LGCP point simulation
simulate_lgcp_points_continuous <- function(Y, coords, covariate_fn, beta, bounds,
                                            scale_intensity = 1000, seed = 123) {
  set.seed(seed)
  d <- ncol(coords)
  volume <- prod(sapply(bounds, function(b) diff(b)))
  
  # Step 1: Sample a homogeneous Poisson process at high enough intensity
  eta_vals <- Y + beta * apply(coords, 1, covariate_fn)
  eta_vals <- pmin(pmax(eta_vals, -10), 10)
  eta_max <- max(eta_vals)
  lambda_max <- scale_intensity * exp(eta_max)
  N_max <- rpois(1, lambda_max * volume)
  
  # Step 2: Uniform samples in the space
  points <- matrix(runif(N_max * d), ncol = d)
  for (i in 1:d) {
    points[, i] <- bounds[[i]][1] + points[, i] * diff(bounds[[i]])
  }
  
  # Step 3: Interpolate eta using nearest neighbor
  nn <- get.knnx(coords, points, k = 1)
  eta_interp <- Y[nn$nn.index] + beta * apply(points, 1, covariate_fn)
  eta_interp <- pmin(pmax(eta_interp, -10), 10)
  lambda_interp <- scale_intensity * exp(eta_interp)
  
  # Step 4: Thinning
  keep <- runif(N_max) < (lambda_interp / lambda_max)
  retained <- points[keep, , drop = FALSE]
  
  cat("Simulated", nrow(retained), "points from continuous LGCP\n")
  return(retained)
}




d <- 3
m <- 100
covariate_fn = function(x) x[1] + 2 * x[2]
beta <- 0.2
estimate_beta = TRUE
scale_intensity <- 1000

bounds <- replicate(d, c(0, 1), simplify = FALSE)
mesh <- build_mesh(d, m, bounds)
fem <- compute_fem_matrices(mesh$points, mesh$simplices)
Q <- assemble_precision_matrix(fem$C, fem$G)
Y <- simulate_latent_field(Q) 
  cat("Latent field variance:", var(Y), "\n")
  
intensity_out <- build_intensity(Y, mesh$points, covariate_fn, beta, scale_intensity)
lgcp_points <- simulate_lgcp_points_continuous(
    Y = Y,
    coords = mesh$points,
    covariate_fn = covariate_fn,
    beta = beta,
    bounds = bounds,
    scale_intensity = scale_intensity
  )
  
library(rgl)
plot3d(lgcp_points[, 1], lgcp_points[, 2], lgcp_points[, 3],
       col = "red", type = "p",
       xlab = "X", ylab = "Y", zlab = "Z", main = "3D LGCP Simulated Points")



d <- 2
m <- 150  # Higher resolution for smoother heatmap
covariate_fn <- function(x) x[1] + 2 * x[2]
beta <- 0.2
estimate_beta <- TRUE
scale_intensity <- 1000

bounds <- replicate(d, c(0, 1), simplify = FALSE)
mesh <- build_mesh(d, m, bounds)
fem <- compute_fem_matrices(mesh$points, mesh$simplices)
Q <- assemble_precision_matrix(fem$C, fem$G)
Y <- simulate_latent_field(Q)
cat("Latent field variance:", var(Y), "\n")

intensity_out <- build_intensity(Y, mesh$points, covariate_fn, beta, scale_intensity)
lgcp_points <- simulate_lgcp_points_continuous(
  Y = Y,
  coords = mesh$points,
  covariate_fn = covariate_fn,
  beta = beta,
  bounds = bounds,
  scale_intensity = scale_intensity
)

# --- Plot: 2D GRF + LGCP points ---
# Reshape GRF to grid
x_vals <- unique(mesh$points[,1])
y_vals <- unique(mesh$points[,2])
z_mat <- matrix(Y, nrow = length(x_vals), ncol = length(y_vals))

image(x_vals, y_vals, z_mat,
      xlab = "X", ylab = "Y", main = "2D GRF and LGCP Points", col = terrain.colors(100))
points(lgcp_points[,1], lgcp_points[,2], pch = 20, col = "red")

