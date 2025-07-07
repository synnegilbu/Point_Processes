# Required libraries
library(readxl)
library(geometry)
library(Matrix)
library(INLA)
library(FNN)


avonet <- read_excel("Documents/Skole/H2024/master/code/Point_Processes/datasets/16586228/AVONET Supplementary dataset 1.xlsx", 
                                             sheet = "AVONET3_BirdTree")
# Filter complete cases
avonet_filtered <- avonet[complete.cases(avonet[, c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width", "Trophic.Level")]), ]

# Extract 3D coordinates (trait space) and covariate
coords_beak_raw <- as.matrix(avonet_filtered[, c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width")])
trophic_covariate <- as.numeric(factor(avonet_filtered$Trophic.Level))


# Scale coordinates to [0,1] for unit cube mesh
coords_beak <- scale(coords_beak_raw, center = apply(coords_beak_raw, 2, min), scale = apply(coords_beak_raw, 2, max) - apply(coords_beak_raw, 2, min))

# Define covariate function using nearest neighbor
covariate_fn <- function(x_row) {
  idx <- get.knnx(coords_beak, matrix(x_row, nrow = 1), k = 1)$nn.index
  trophic_covariate[idx]
}

# Helper functions
jitter_points <- function(pts, eps = 1e-4) {
  pts + matrix(runif(length(pts), -eps, eps), ncol = ncol(pts))
}

compute_fem_matrices <- function(points, simplices) {
  n <- nrow(points); d <- ncol(points)
  C <- Matrix(0, n, n, sparse = TRUE)
  G <- Matrix(0, n, n, sparse = TRUE)
  for (i in 1:nrow(simplices)) {
    idx <- simplices[i, ]
    verts <- points[idx, , drop = FALSE]
    T <- t(verts[-1, , drop = FALSE] - matrix(verts[1, ], d, d, byrow = TRUE))
    detT <- det(T); if (abs(detT) < 1e-12) next
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

assemble_precision_matrix <- function(C, G, tau = 1, kappa = 1, jitter = 1e-2) {
  Q <- tau^2 * (kappa^2 * C + G)
  Q <- forceSymmetric(Q)
  diag_vals <- diag(Q)
  diag(Q)[diag_vals < jitter] <- jitter
  return(Q)
}

simulate_latent_field <- function(Q) {
  as.vector(inla.qsample(n = 1, Q = Q))
}

check_matrix_sanity <- function(Q) {
  dvals <- diag(Q)
  cat("Min diag(Q):", min(dvals), "\n")
  cat("Any NA?:", any(is.na(Q)), "\n")
  cat("Any Inf?:", any(is.infinite(Q)), "\n")
  if (any(dvals < 0)) stop("Negative diagonal entry in Q — aborting.")
}

build_projector_matrix <- function(mesh_points, simplices, eval_points) {
  d <- ncol(mesh_points)
  A <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(nrow(eval_points), nrow(mesh_points)))
  ts <- tsearchn(mesh_points, simplices, eval_points, bary = TRUE)
  valid <- which(!is.na(ts$idx))
  for (i in valid) {
    simplex_id <- ts$idx[i]
    bary <- ts$p[i, ]; if (any(is.na(bary))) next
    verts <- simplices[simplex_id, ]
    A[i, verts] <- bary
  }
  A
}

construct_likelihood_data <- function(mesh, observed_points, covariate_fn, bounds) {
  mesh_points <- mesh$points; n_mesh <- nrow(mesh_points); n_obs <- nrow(observed_points)
  total_volume <- prod(sapply(bounds, function(b) diff(b)))
  alpha_weights <- rep(total_volume / n_mesh, n_mesh)
  locations <- rbind(mesh_points, observed_points)
  A <- build_projector_matrix(mesh_points, mesh$simplices, locations)
  cov_values <- apply(locations, 1, covariate_fn)
  y <- c(rep(0, n_mesh), rep(1, n_obs))
  weight <- c(alpha_weights, rep(0, n_obs))
  idx <- c(1:n_mesh, rep(n_mesh + 1L, n_obs))
  list(y = y, weight = weight, covariate = cov_values, idx = idx, A = A)
}

run_inla_continuous <- function(likelihood_data, Q, estimate_beta = TRUE, beta = 0.0) {
  y <- as.numeric(likelihood_data$y)
  weight <- as.numeric(likelihood_data$weight)
  covariate <- as.numeric(likelihood_data$covariate)
  A <- likelihood_data$A
  offset <- if (!estimate_beta) beta * covariate else 0
  idx_latent <- 1:max(likelihood_data$idx)
  stk <- inla.stack(
    data = list(y = y, E = weight + 1e-10, offset = offset),
    A = list(A, 1),
    effects = list(idx = 1:ncol(A), covariate = covariate),
    tag = "spatial"
  )
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
  return(list(
    mesh_points = mesh$points,  
    estimated_field = result$summary.random$idx$mean,
    fixed_effects = result$summary.fixed,
    result = result
  ))
  
}

# Mesh creation using observed points
enforce_min_spacing <- function(points, min_dist = 0.03) {
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

run_spde_lgcp_for_traits <- function(coords_raw,
                                     covariate_vals,
                                     covariate_fn = NULL,
                                     d = ncol(coords_raw),
                                     beta = 1.0,
                                     m_coarse = 5,
                                     m_fine = 5,
                                     min_spacing = 0.15,
                                     scale_intensity = 50,
                                     estimate_beta = TRUE) {
  
  # Bounds in unit cube
  bounds <- replicate(d, c(0, 1), simplify = FALSE)
  
  # Scale coordinates to [0,1]^d
  coords_scaled <- scale(coords_raw, center = apply(coords_raw, 2, min),
                         scale = apply(coords_raw, 2, max) - apply(coords_raw, 2, min))
  
  # Default covariate function via nearest neighbor
  if (is.null(covariate_fn)) {
    covariate_fn <- function(x_row) {
      idx <- get.knnx(coords_scaled, matrix(x_row, nrow = 1), k = 1)$nn.index
      covariate_vals[idx]
    }
  }
  
  # Coarse mesh
  grid_axes <- replicate(d, seq(0, 1, length.out = m_coarse), simplify = FALSE)
  initial_mesh_grid <- do.call(expand.grid, grid_axes)
  initial_mesh_points <- jitter_points(as.matrix(initial_mesh_grid))
  initial_simplices <- delaunayn(initial_mesh_points, options = "QJ")
  initial_mesh <- list(points = initial_mesh_points, simplices = initial_simplices)
  
  fem_initial <- compute_fem_matrices(initial_mesh$points, initial_mesh$simplices)
  Q_initial <- assemble_precision_matrix(fem_initial$C, fem_initial$G)
  check_matrix_sanity(Q_initial)
  Y_initial <- simulate_latent_field(Q_initial)
  Y_initial <- Y_initial - mean(Y_initial)
  
  # Simulate LGCP points
  lgcp_points_initial <- simulate_lgcp_points_continuous(
    Y_initial, initial_mesh$points, covariate_fn, beta,
    bounds = bounds, scale_intensity = scale_intensity
  )
  
  # Final mesh
  fine_grid_axes <- replicate(d, seq(0, 1, length.out = m_fine), simplify = FALSE)
  regular_grid <- do.call(expand.grid, fine_grid_axes)
  combined <- rbind(lgcp_points_initial, coords_scaled, as.matrix(regular_grid))
  filtered_points <- enforce_min_spacing(combined, min_dist = min_spacing)
  mesh_points <- jitter_points(filtered_points)
  simplices <- delaunayn(mesh_points, options = "QJ")
  mesh <- list(points = mesh_points, simplices = simplices)
  
  # Final latent field
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  Q <- assemble_precision_matrix(fem$C, fem$G)
  check_matrix_sanity(Q)
  Y <- simulate_latent_field(Q)
  Y <- Y - mean(Y)
  
  # Likelihood data and INLA model
  likelihood_data <- construct_likelihood_data(mesh, coords_scaled, covariate_fn, bounds)
  result <- run_inla_continuous(likelihood_data, Q, estimate_beta = estimate_beta, beta = beta)
  
  return(list(
    dimension = d,
    mesh = mesh,
    latent_field = Y,
    estimated_field = result$estimated_field,
    fixed_effects = result$fixed_effects,
    result = result
  ))
}


coords_raw <- as.matrix(avonet[, c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width")])
trophic_covariate <- as.numeric(factor(avonet$Trophic.Level))

result3d_refined <- run_spde_lgcp_for_traits(
  coords_raw = coords_raw,
  covariate_vals = trophic_covariate,
  d = 3,
  beta = 1.0,
  estimate_beta = TRUE,
  m_coarse = 7,
  m_fine = 8,
  min_spacing = 0.10,
  scale_intensity = 200
)


result4d_fast <- run_spde_lgcp_for_traits(
  coords_raw = coords_beak_raw4d,
  covariate_vals = trophic_covariate,
  d = 4,
  beta = 0.0,
  estimate_beta = FALSE,
  m_coarse = 4,
  m_fine = 5,
  min_spacing = 0.15,
  scale_intensity = 30
)



check_quadratic_field_fit <- function(coords, field_vals) {
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  
  d <- ncol(coords)
  if (d < 1) stop("Invalid coordinate dimension.")
  
  df <- as.data.frame(coords)
  colnames(df) <- paste0("X", 1:d)
  df$Y <- field_vals
  
  # Build quadratic formula dynamically
  terms <- paste0("X", 1:d)
  squares <- paste0("I(", terms, "^2)")
  interactions <- combn(terms, 2, FUN = function(x) paste(x, collapse = ":"))
  rhs <- paste(c(terms, squares, interactions), collapse = " + ")
  formula <- as.formula(paste("Y ~", rhs))
  
  # Fit model
  model <- lm(formula, data = df)
  summ <- summary(model)
  
  return(list(
    model = model,
    summary = summ,
    r_squared = summ$r.squared,
    adj_r_squared = summ$adj.r.squared
  ))
}

quad3d <- check_quadratic_field_fit(result3d_refined$mesh$points, result3d_refined$estimated_field)
quad4d <- check_quadratic_field_fit(result4d_fast$mesh$points, result4d_fast$estimated_field)

print(quad3d)
print(quad4d)

print(result3d_refined$fixed_effects)
print(result4d_fast$fixed_effects)
