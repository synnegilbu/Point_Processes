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


# INLA inference
run_inla_inference <- function(counts, covariate, Q, estimate_beta = TRUE) {
  df <- data.frame(y = counts, idx = 1:length(counts), cov = covariate)
  
  fmla <- if (estimate_beta) {
    y ~ cov + f(idx, model = "generic0", Cmatrix = Q,
                hyper = list(prec = list(initial = 2, fixed = FALSE)))
  } else {
    y ~ f(idx, model = "generic0", Cmatrix = Q,
          hyper = list(prec = list(initial = 2, fixed = FALSE)))
  }
  
  inla(
    fmla,
    family = "nbinomial",
    data = df,
    control.predictor = list(compute = TRUE),
    control.inla = list(strategy = "laplace")
  )
}

# Pipeline
run_spde_lgcp_pipeline <- function(
    d = 3, m = 10, covariate_fn = function(x) sum(x), beta = 0.1,
    estimate_beta = TRUE, scale_intensity = 1000
) {
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
  
  counts <- tabulate(get.knnx(mesh$points, lgcp_points, k = 1)$nn.index, nbins = nrow(mesh$points))
  
  
  cat("Sampled points:", nrow(lgcp_points), "\n")
  cat("Non-zero bins:", sum(counts > 0), "/", length(counts), "\n")
  
  result <- run_inla_inference(counts, intensity_out$covariate, Q, estimate_beta)
  
  list(
    mesh_points = mesh$points,
    latent_field = Y,
    estimated_field = result$summary.random$idx$mean,
    beta_estimate = if (estimate_beta) result$summary.fixed else NULL,
    counts = counts,
    result = result
  )
}

res <- run_spde_lgcp_pipeline(
  d = 3,
  m = 10,
  covariate_fn = function(x) x[1] + 2 * x[2],  
  beta = 0.5,
  estimate_beta = TRUE
)

plot(res$latent_field, res$estimated_field,
     xlab = "True Latent Field",
     ylab = "INLA Estimated Field",
     main = "SPDE LGCP Recovery",
     col = rgb(0, 0, 1, 0.3), pch = 16)
abline(0, 1, col = "red")

print(res$beta_estimate)
cor(res$latent_field, res$estimated_field)      # Correlation
mean((res$latent_field - res$estimated_field)^2) # MSE
plot(res$latent_field, res$latent_field - res$estimated_field,
     xlab = "True Field", ylab = "Residual", main = "Residual Plot")
abline(h = 0, col = "red")








# --- Run 4D Example ---
res4d <- run_spde_lgcp_pipeline(
  d = 4,
  m = 4,  # 4^4 = 256 nodes — a good balance of resolution vs compute
  covariate_fn = function(x) x[1] + 0.5 * x[2] - x[3] + 0.2 * x[4],
  beta = 0.7,
  estimate_beta = TRUE,
  scale_intensity = 1500
)

# 1. Extract posterior mean of latent field
Y_hat <- res4d$estimated_field

# 2. Recompute covariate values (must match what was used in simulation)
covariate_vals <- apply(res4d$mesh_points, 1, function(x) x[1] + 0.5 * x[2] - x[3] + 0.2 * x[4])

# 3. Get estimated beta from INLA output
beta_hat <- res4d$beta_estimate["cov", "mean"]

# 4. Compute estimated log-intensity and intensity
eta_hat <- Y_hat + beta_hat * covariate_vals
lambda_hat <- exp(eta_hat)

# Correlation between estimated and true log-intensity
cat("Corr(η):", cor(eta_hat, eta_true), "\n")

# Correlation between estimated and true intensity
cat("Corr(λ):", cor(lambda_hat, lambda_true), "\n")

# True latent field from simulation
Y_true <- res4d$latent_field

# True beta used during simulation
beta_true <- 0.7

# Compute true log-intensity and intensity
eta_true <- Y_true + beta_true * covariate_vals
lambda_true <- exp(eta_true)


# Compare visually
plot(lambda_true, lambda_hat,
     xlab = "True λ", ylab = "Estimated λ",
     main = "Estimated vs True Intensity (4D)",
     pch = 16, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red")




