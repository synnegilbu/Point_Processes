# Required libraries
library(geometry)
library(Matrix)
library(INLA)
library(FNN)

# Jitter points to avoid degenerate simplices
jitter_points <- function(pts, eps = 1e-4) {
  pts + matrix(runif(length(pts), -eps, eps), ncol = ncol(pts))
}

# Build mesh in d dimensions
build_mesh <- function(d, m, bounds) {
  grid_axes <- lapply(bounds, function(lim) seq(lim[1], lim[2], length.out = m))
  grid_points <- expand.grid(grid_axes)
  pts <- jitter_points(as.matrix(grid_points))
  simplices <- delaunayn(pts, options = "QJ")
  list(points = pts, simplices = simplices)
}

# Compute FEM matrices
compute_fem_matrices <- function(points, simplices) {
  n <- nrow(points)
  d <- ncol(points)
  C <- Matrix(0, n, n, sparse = TRUE)
  G <- Matrix(0, n, n, sparse = TRUE)
  for (i in 1:nrow(simplices)) {
    idx <- simplices[i, ]
    verts <- points[idx, , drop = FALSE]
    T <- t(verts[-1, , drop = FALSE] - matrix(verts[1, ], d, d, byrow = TRUE))
    detT <- det(T)
    if (abs(detT) < 1e-12) next
    vol <- abs(detT) / factorial(d)
    Ghat <- cbind(-1, diag(d))
    grads <- solve(t(T), Ghat)
    GK <- vol * t(grads) %*% grads
    MK <- matrix(1, d+1, d+1); diag(MK) <- 2
    MK <- MK * vol / ((d + 1) * (d + 2))
    C[idx, idx] <- C[idx, idx] + MK
    G[idx, idx] <- G[idx, idx] + GK
  }
  list(C = C, G = G)
}

# Assemble precision matrix
assemble_precision_matrix <- function(C, G, tau = 0.1, kappa = 5, jitter = 1e-4) {
  Q <- tau^2 * (kappa^2 * C + G)
  Q <- forceSymmetric(Q)
  diag_vals <- diag(Q)
  diag(Q)[diag_vals < jitter] <- jitter
  return(Q)
}

# Sanity check for Q
check_matrix_sanity <- function(Q) {
  dvals <- diag(Q)
  cat("Min diag(Q):", min(dvals), "\n")
  cat("Any NA?:", any(is.na(Q)), "\n")
  cat("Any Inf?:", any(is.infinite(Q)), "\n")
  if (any(dvals < 0)) stop("Negative diagonal entry in Q — aborting.")
}

# Simulate latent field
simulate_latent_field <- function(Q) {
  as.vector(inla.qsample(n = 1, Q = Q))
}

# Simulate LGCP points
simulate_lgcp_points_continuous <- function(Y, coords, covariate_fn, beta, bounds, scale_intensity = 2000, seed = 123) {
  set.seed(seed)
  d <- ncol(coords)
  volume <- prod(sapply(bounds, function(b) diff(b)))
  eta_vals <- Y + beta * apply(coords, 1, covariate_fn)
  eta_vals <- pmin(pmax(eta_vals, -10), 10)
  lambda_max <- min(scale_intensity * exp(6), 1e5)
  N_max <- rpois(1, lambda_max * volume)
  points <- matrix(runif(N_max * d), ncol = d)
  for (i in 1:d) points[, i] <- bounds[[i]][1] + points[, i] * diff(bounds[[i]])
  nn <- get.knnx(coords, points, k = 1)
  eta_interp <- Y[nn$nn.index] + beta * apply(points, 1, covariate_fn)
  eta_interp <- pmin(pmax(eta_interp, -10), 10)
  lambda_interp <- scale_intensity * exp(eta_interp)
  keep <- runif(N_max) < (lambda_interp / lambda_max)
  points[keep, , drop = FALSE]
}

# Projector matrix A
build_projector_matrix <- function(mesh_points, simplices, eval_points) {
  d <- ncol(mesh_points)
  A <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                    dims = c(nrow(eval_points), nrow(mesh_points)))
  ts <- tsearchn(mesh_points, simplices, eval_points, bary = TRUE)
  valid <- which(!is.na(ts$idx))
  for (i in valid) {
    simplex_id <- ts$idx[i]
    bary <- ts$p[i, ]
    if (any(is.na(bary))) next
    verts <- simplices[simplex_id, ]
    A[i, verts] <- bary
  }
  A
}

# Construct likelihood data
construct_likelihood_data <- function(mesh, observed_points, covariate_fn, alpha_weights) {
  mesh_points <- mesh$points
  n_mesh <- nrow(mesh_points)
  n_obs <- nrow(observed_points)
  locations <- rbind(mesh_points, observed_points)
  A <- build_projector_matrix(mesh_points, mesh$simplices, locations)
  cov_values <- apply(locations, 1, covariate_fn)
  y <- c(rep(0, n_mesh), rep(1, n_obs))
  weight <- c(alpha_weights, rep(0, n_obs))
  idx <- c(1:n_mesh, rep(NA, n_obs))
  list(y = y, weight = weight, covariate = cov_values, idx = idx, A = A)
}

# Run INLA model
run_inla_continuous <- function(likelihood_data, Q, estimate_beta = TRUE, beta = 0.0) {
  y <- as.numeric(likelihood_data$y)
  weight <- as.numeric(likelihood_data$weight)
  covariate <- as.numeric(likelihood_data$covariate)
  mesh_idx <- likelihood_data$idx
  A <- likelihood_data$A
  
  offset <- beta * covariate
  n_mesh <- max(mesh_idx, na.rm = TRUE)
  
  # Index of latent field values: mesh nodes only
  idx_latent <- 1:n_mesh
  
  # Build INLA stack
  stk <- inla.stack(
    data = list(y = y, E = weight + 1e-10, offset = offset),
    A = list(A, 1),
    effects = list(
      idx = 1:ncol(A),
      covariate = covariate
    ),
    tag = "spatial"
  )
  
  # Try running INLA with error handling
  result <- tryCatch({
    inla(
      formula = y ~ covariate + f(idx, model = "generic0", Cmatrix = Q),
      family = "poisson",
      data = inla.stack.data(stk),
      E = inla.stack.data(stk)$E,
      offset = inla.stack.data(stk)$offset,
      control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
      control.inla = list(strategy = "laplace"),
      verbose = TRUE
    )
  }, error = function(e) {
    cat("INLA failed with error:\n", conditionMessage(e), "\n")
    return(NULL)
  })
  
  return(result)
}


# Pipeline
run_spde_lgcp_pipeline_continuous <- function(
    d = 3,
    m = 10,  # finer mesh
    covariate_fn = function(x) x[1],
    beta = 1.0,
    estimate_beta = TRUE,
    scale_intensity = 3000
) {
  bounds <- replicate(d, c(0, 1), simplify = FALSE)
  
  # Step 1: Build mesh
  mesh <- build_mesh(d, m, bounds)
  
  # Step 2: Compute FEM matrices
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  
  # Step 3: Precision matrix
  Q <- assemble_precision_matrix(fem$C, fem$G)
  check_matrix_sanity(Q)
  
  # Step 4: Simulate unscaled latent field
  Y <- simulate_latent_field(Q)
  
  # Step 5: Simulate LGCP points using true field
  lgcp_points <- simulate_lgcp_points_continuous(
    Y, mesh$points, covariate_fn, beta, bounds, scale_intensity
  )
  
  # Step 6: Construct likelihood
  alpha_weights <- rep(prod(sapply(bounds, function(b) diff(b))) / nrow(mesh$points),
                       nrow(mesh$points))
  likelihood_data <- construct_likelihood_data(mesh, lgcp_points, covariate_fn, alpha_weights)
  
  # Step 7: Run INLA model
  result <- run_inla_continuous(likelihood_data, Q, estimate_beta, beta = beta)
  
  list(
    mesh_points = mesh$points,
    latent_field = Y,
    estimated_field = result$summary.random$idx$mean,
    beta_estimate = if (estimate_beta) result$summary.fixed else NULL,
    observed_points = lgcp_points,
    result = result
  )
}


# Example run
set.seed(123)
result3d <- run_spde_lgcp_pipeline_continuous(
  d = 3,
  m = 10,
  covariate_fn = function(x) x[1],
  beta = 1.5,
  estimate_beta = TRUE,
  scale_intensity = 3000
)



cat("Estimated beta (3D):\n")
print(result3d$beta_estimate)

cat("\nLatent field RMSE (3D):\n")
rmse <- sqrt(mean((result3d$latent_field - result3d$estimated_field)^2))
print(rmse)

cat("\nCorrelation between true and estimated field (3D):\n")
print(cor(result3d$latent_field, result3d$estimated_field))

plot_field_slice <- function(result, slice_dim = 3, slice_val = 0.5, tol = 0.05, title_prefix = "Latent Field") {
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  
  df <- as.data.frame(result$mesh_points)
  colnames(df) <- c("x", "y", "z")
  df$true <- result$latent_field
  df$estimated <- log(result$estimated_field)
  
  dims <- c("x", "y", "z")
  fixed_dim <- dims[slice_dim]
  other_dims <- setdiff(dims, fixed_dim)
  
  df_slice <- df %>% filter(abs(.data[[fixed_dim]] - slice_val) < tol)
  
  cat("Number of points in slice:", nrow(df_slice), "\n")
  if (nrow(df_slice) == 0) {
    warning("No points found near the specified slice value. Try increasing `tol`.")
    return(NULL)
  }
  
  fill_limits <- range(c(df_slice$true, df_slice$estimated), na.rm = TRUE)
  
  p1 <- ggplot(df_slice, aes(x = .data[[other_dims[1]]], y = .data[[other_dims[2]]], fill = true)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(limits = fill_limits) +
    labs(title = paste(title_prefix, "- True"), x = other_dims[1], y = other_dims[2])
  
  p2 <- ggplot(df_slice, aes(x = .data[[other_dims[1]]], y = .data[[other_dims[2]]], fill = estimated)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(limits = fill_limits) +
    labs(title = paste(title_prefix, "- Estimated"), x = other_dims[1], y = other_dims[2])
  
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}


plot_field_slice(result3d, slice_dim = 3, slice_val = 0.5, tol = 0.1)


df <- as.data.frame(result3d$mesh_points)
colnames(df) <- c("x", "y", "z")
df$true <- result3d$latent_field
df$true <- result$latent_field           # Already log-intensity
df$estimated <- log(result3d$estimated_field)  # Convert back to log-scale


df_slice <- df %>% filter(abs(z - 0.5) < 0.2)
nrow(df_slice)
head(df_slice)





# Updated 2D simulation run
test_result_2d <- run_spde_lgcp_pipeline_continuous(
  d = 2,
  m = 40,
  covariate_fn = function(x) x[1],
  beta = 1.0,
  estimate_beta = TRUE,
  scale_intensity = 10000
)


plot_latent_vs_estimated <- function(result, title_prefix = "Field") {
  library(ggplot2)
  library(akima)
  library(gridExtra)
  
  df <- as.data.frame(result$mesh_points)
  colnames(df) <- c("x", "y")
  df$true <- result$latent_field
  df$estimated <- result$estimated_field
  
  # Interpolate onto a regular grid
  interp_true <- with(df, akima::interp(x, y, true, duplicate = "mean"))
  interp_est <- with(df, akima::interp(x, y, estimated, duplicate = "mean"))
  
  df_true <- expand.grid(x = interp_true$x, y = interp_true$y)
  df_true$z <- as.vector(interp_true$z)
  
  df_est <- expand.grid(x = interp_est$x, y = interp_est$y)
  df_est$z <- as.vector(interp_est$z)
  
  p1 <- ggplot(df_true, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    coord_equal() +
    scale_fill_viridis_c() +
    ggtitle(paste0(title_prefix, ": True Latent Field"))
  
  p2 <- ggplot(df_est, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = TRUE) +
    coord_equal() +
    scale_fill_viridis_c() +
    ggtitle(paste0(title_prefix, ": Estimated Latent Field"))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}


# Visual comparison
plot_latent_vs_estimated(test_result_2d, title_prefix = "2D LGCP")


cat("Estimated beta (2D):\n")
print(test_result_2d$beta_estimate)

cat("\nLatent field RMSE (2D):\n")
rmse_2d <- sqrt(mean((test_result_2d$latent_field - test_result_2d$estimated_field)^2))
print(rmse_2d)

cat("\nCorrelation between true and estimated field (2D):\n")
print(cor(test_result_2d$latent_field, test_result_2d$estimated_field))

