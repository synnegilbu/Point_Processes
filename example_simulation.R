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
assemble_precision_matrix <- function(C, G, tau = 1, kappa = 1, jitter = 1e-2) {
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
construct_likelihood_data <- function(mesh, observed_points, covariate_fn, bounds) {
  mesh_points <- mesh$points
  n_mesh <- nrow(mesh_points)
  n_obs <- nrow(observed_points)
  total_volume <- prod(sapply(bounds, function(b) diff(b)))
  alpha_weights <- rep(total_volume / n_mesh, n_mesh)  
  locations <- rbind(mesh_points, observed_points)
  A <- build_projector_matrix(mesh_points, mesh$simplices, locations)
  cov_values <- apply(locations, 1, covariate_fn)
  y <- c(rep(0, n_mesh), rep(1, n_obs))
  weight <- c(alpha_weights, rep(0, n_obs))
  idx <- c(1:n_mesh, rep(n_mesh + 1L, n_obs))  
  list(
    y = y,
    weight = weight,
    covariate = cov_values,
    idx = idx,
    A = A
  )
}


# Run INLA model
run_inla_continuous <- function(likelihood_data, Q, estimate_beta = TRUE, beta = 0.0) {
  y <- as.numeric(likelihood_data$y)
  weight <- as.numeric(likelihood_data$weight)
  covariate <- as.numeric(likelihood_data$covariate)
  mesh_idx <- likelihood_data$idx
  A <- likelihood_data$A
  
  offset <- if (!estimate_beta) beta * covariate else 0
  n_mesh <- max(mesh_idx, na.rm = TRUE)
  idx_latent <- 1:n_mesh
  
  stk <- inla.stack(
    data = list(y = y, E = weight + 1e-10, offset = offset),
    A = list(A, 1),
    effects = list(
      idx = 1:ncol(A),
      covariate = covariate
    ),
    tag = "spatial"
  )
  cat("Minimum eigenvalue of Q: ")
  print(min(eigen(as.matrix(Q), only.values = TRUE)$values))
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



run_spde_lgcp_pipeline_continuous <- function(
    d = 3,
    m = 10,
    covariate_fn = function(x) x[1],
    beta = 1.0,
    estimate_beta = TRUE,
    scale_intensity = 3000
) {
  bounds <- replicate(d, c(0, 1), simplify = FALSE)
  
  # ---- First pass: initial mesh and LGCP simulation ----
  mesh_initial <- build_mesh(d, m, bounds)
  fem_initial <- compute_fem_matrices(mesh_initial$points, mesh_initial$simplices)
  Q_initial <- assemble_precision_matrix(fem_initial$C, fem_initial$G)
  check_matrix_sanity(Q_initial)
  Y_initial <- simulate_latent_field(Q_initial)
  
  lgcp_points_initial <- simulate_lgcp_points_continuous(
    Y_initial, mesh_initial$points, covariate_fn, beta, bounds, scale_intensity
  )
  
  # ---- Second pass: enforce spacing and build new mesh ----
  enforce_min_spacing <- function(points, min_dist = 0.03) {
    if (!is.matrix(points)) stop("Points must be a matrix.")
    if (ncol(points) < 2) stop("Points must be at least 2D.")
    keep_indices <- c(1)
    for (i in 2:nrow(points)) {
      pt <- points[i, , drop = FALSE]
      kept <- points[keep_indices, , drop = FALSE]
      pt_rep <- matrix(rep(pt, each = nrow(kept)), nrow = nrow(kept))
      dists <- sqrt(rowSums((pt_rep - kept)^2))
      if (all(dists > min_dist)) keep_indices <- c(keep_indices, i)
    }
    points[keep_indices, , drop = FALSE]
  }
  
  regular_grid <- build_mesh(d = d, m = 15, bounds = bounds)$points
  combined <- rbind(lgcp_points_initial, regular_grid)
  filtered_points <- enforce_min_spacing(combined, min_dist = 0.025)
  new_mesh_points <- jitter_points(filtered_points, eps = 1e-4)
  new_simplices <- delaunayn(new_mesh_points, options = "QJ")
  mesh <- list(points = new_mesh_points, simplices = new_simplices)
  trimesh(mesh$simplices, mesh$points)
  
  # ---- Volume diagnostics (optional) ----
  volumes <- apply(mesh$simplices, 1, function(i) {
    verts <- mesh$points[i, , drop = FALSE]
    if (d == 2) {
      abs(det(cbind(verts[2, ] - verts[1, ], verts[3, ] - verts[1, ]))) / 2
    } else if (d == 3) {
      abs(det(verts[2:4, ] - matrix(verts[1, ], 3, 3, byrow = TRUE))) / 6
    } else {
      NA
    }
  })
  print(summary(volumes))
  
  # ---- Final mesh, latent field, and data ----
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  Q <- assemble_precision_matrix(fem$C, fem$G,, tau=1)
  check_matrix_sanity(Q)
  Y <- simulate_latent_field(Q)
  Y <- Y - mean(Y)
  
  
  lgcp_points <- simulate_lgcp_points_continuous(
    Y, mesh$points, covariate_fn, beta, bounds, scale_intensity
  )
  points(lgcp_points)
  
  likelihood_data <- construct_likelihood_data(mesh, lgcp_points, covariate_fn, bounds)
  
  result <- run_inla_continuous(likelihood_data, Q, estimate_beta = FALSE, beta = 1.0)

  
  list(
    mesh_points = mesh$points,
    latent_field = Y,
    estimated_field = result$summary.random$idx$mean,
    beta_estimate = if (estimate_beta) result$summary.fixed else NULL,
    observed_points = lgcp_points,
    likelihood_data = likelihood_data,
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
  scale_intensity = 500 
)



test_result_2d <- run_spde_lgcp_pipeline_continuous(
  d = 2,
  m = 40,
  covariate_fn = function(x) x[1],
  beta = 1.0,
  estimate_beta = TRUE,
  scale_intensity = 500
)


