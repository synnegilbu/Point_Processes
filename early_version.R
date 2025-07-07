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
    T <- t(verts[-1, , drop = FALSE] - matrix(verts[1, ], d, d, byrow = TRUE))  # d x d
    detT <- det(T)
    if (abs(detT) < 1e-12) {
      warning(paste("Skipping degenerate simplex at row", i, "with det =", detT))
      next
    }
    
    vol <- abs(detT) / factorial(d)
    
    # Gradients of basis functions in physical space
    Ghat <- cbind(-1, diag(d))  # Gradients on reference element
    grads <- solve(t(T), Ghat)  # Each column is grad(phi_i)
    
    
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
assemble_precision_matrix <- function(C, G, tau = 0.1, kappa = 2, jitter = 1e-5) {
  Q <- tau^2 * (kappa^2 * C + G)
  Q + Diagonal(nrow(Q), jitter)
}

# Latent field
simulate_latent_field <- function(Q) {
  as.vector(inla.qsample(n = 1, Q = Q))
}

# Intensity function
build_intensity <- function(Y, coords, covariate_fn, beta, scale_intensity = 2000) {
  cov_vals <- apply(coords, 1, covariate_fn)
  eta <- Y + beta * cov_vals
  eta <- pmin(pmax(eta, -6), 6)  # clamp for stability
  lambda <- scale_intensity * exp(eta)
  list(lambda = lambda, covariate = cov_vals)
}

# LGCP point simulation
simulate_lgcp_points_continuous <- function(Y, coords, covariate_fn, beta, bounds,
                                            scale_intensity = 2000, seed = 123) {
  set.seed(seed)
  d <- ncol(coords)
  volume <- prod(sapply(bounds, function(b) diff(b)))
  
  # Step 1: Sample a homogeneous Poisson process at high enough intensity
  eta_vals <- Y + beta * apply(coords, 1, covariate_fn)
  eta_vals <- pmin(pmax(eta_vals, -10), 10)
  eta_max <- max(eta_vals)
  lambda_max <- min(scale_intensity * exp(6), 1e5)  # clamp max
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
    estimate_beta = TRUE, scale_intensity = 2000
) {
  bounds <- replicate(d, c(0, 1), simplify = FALSE)
  mesh <- build_mesh(d, m, bounds)
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  Q <- assemble_precision_matrix(fem$C, fem$G)
  Y <- simulate_latent_field(Q) 
  Y <- scale(Y)  # zero mean, unit variance
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

# res <- run_spde_lgcp_pipeline(
#   d = 3,
#   m = 10,
#   covariate_fn = function(x) x[1] + 2 * x[2],  
#   beta = 0.5,
#   estimate_beta = TRUE
# )
# 
# hist(res$latent_field, breaks = 50, main = "Latent Field Distribution", xlab = "Y", col = "skyblue")
# 
# plot(res$latent_field, res$estimated_field, pch = 16, cex = 0.5,
#      xlab = "True Y", ylab = "Estimated Y", main = "Latent vs Estimated Field")
# abline(0, 1, col = "red")
# 
# residuals <- res$latent_field - res$estimated_field
# hist(residuals, breaks = 50, col = "lightcoral", main = "Residuals", xlab = "Y - Y_hat")
# hist(res$counts, breaks = 30, main = "Histogram of Poisson Counts", xlab = "Counts", col = "grey")


benchmark_mse_dim1 <- function(reps = 3, m = 10, beta = FALSE) {
  mse_values <- numeric(reps)
  
  for (rep in 1:reps) {
    cat("1D Repetition", rep, "\n")
    
    res <- run_spde_lgcp_pipeline(
      d = 1,
      m = m,
      covariate_fn = function(x) x[1],
      beta = beta,
      estimate_beta = TRUE
    )
    
    mse_values[rep] <- mean((res$latent_field - res$estimated_field)^2)
    print(cor(res$latent_field, res$estimated_field))
    
  }
  
  data.frame(dimension = 1, repetition = 1:reps, latent_mse = mse_values)
}
benchmark_mse_dim2 <- function(reps = 3, m = 10, beta = FALSE) {
  mse_values <- numeric(reps)
  
  for (rep in 1:reps) {
    cat("2D Repetition", rep, "\n")
    
    res <- run_spde_lgcp_pipeline(
      d = 2,
      m = m,
      covariate_fn = function(x) x[1],
      beta = beta,
      estimate_beta = TRUE
    )
    
    mse_values[rep] <- mean((res$latent_field - res$estimated_field)^2)
    print(cor(res$latent_field, res$estimated_field))

  }
  
  data.frame(dimension = 2, repetition = 1:reps, latent_mse = mse_values)
}
benchmark_mse_dim3 <- function(reps = 3, m = 10, beta = FALSE) {
  mse_values <- numeric(reps)
  
  for (rep in 1:reps) {
    cat("3D Repetition", rep, "\n")
    
    res <- run_spde_lgcp_pipeline(
      d = 3,
      m = m,
      covariate_fn = function(x) x[1],
      beta = beta,
      estimate_beta = TRUE
    )
    
    mse_values[rep] <- mean((res$latent_field - res$estimated_field)^2)
  }
  
  data.frame(dimension = 3, repetition = 1:reps, latent_mse = mse_values)
}
benchmark_mse_dim4 <- function(reps = 3, m = 10, beta = FALSE) {
  mse_values <- numeric(reps)
  
  for (rep in 1:reps) {
    cat("4D Repetition", rep, "\n")
    
    res <- run_spde_lgcp_pipeline(
      d = 4,
      m = m,
      covariate_fn = function(x) x[1],
      beta = beta,
      estimate_beta = TRUE
    )
    
    mse_values[rep] <- mean((res$latent_field - res$estimated_field)^2)
  }
  
  data.frame(dimension = 3, repetition = 1:reps, latent_mse = mse_values)
}

set.seed(42)

mse1 <- benchmark_mse_dim1()
mse2 <- benchmark_mse_dim2()
print(rbind(mean(mse1$latent_mse), mean(mse2$latent_mse)))

mse3 <- benchmark_mse_dim3()
mse4 <- benchmark_mse_dim4()

rbind(mean(mse1$latent_mse), mean(mse2$latent_mse), mean(mse3$latent_mse), mean(mse4$latent_mse))





res2d <- run_spde_lgcp_pipeline(
  d = 2,
  m = 100,  # use higher resolution for nicer plots
  covariate_fn = function(x) x[1] + 2 * x[2],
  beta = 0.5,
  estimate_beta = TRUE,
  scale_intensity = 7000
)


library(fields)
# Prepare data
x <- unique(res2d$mesh_points[, 1])
y <- unique(res2d$mesh_points[, 2])
z_true <- matrix(res2d$latent_field, nrow = length(x), byrow = FALSE)
z_est <- matrix(res2d$estimated_field, nrow = length(x), byrow = FALSE)

# Shared color limits
zlim <- range(c(z_true, z_est), na.rm = TRUE)

# Set layout
layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
par(mar = c(4, 4, 3, 1))  # Margins: bottom, left, top, right

# Plot true field
image.plot(x, y, z_true,
           main = "True Latent Field", xlab = "x", ylab = "y", asp = 1,
           zlim = zlim)

# Plot estimated field
image.plot(x, y, z_est,
           main = "Estimated Latent Field (INLA)", xlab = "x", ylab = "y", asp = 1,
           zlim = zlim)

pts <- res2d$mesh_points
Y_true <- res2d$latent_field
Y_est  <- res2d$estimated_field
residuals <- Y_true - Y_est

cat("RMSE:", sqrt(mean(residuals^2)), "\n")
cat("Mean residual:", mean(residuals), "\n")
cat("SD of residuals:", sd(residuals), "\n")


