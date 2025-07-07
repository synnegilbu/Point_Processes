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

# Build mesh from AVONET beak space
filtered_points <- enforce_min_spacing(coords_beak, min_dist = 0.03)
mesh_points <- jitter_points(filtered_points)
simplices <- delaunayn(mesh_points, options = "QJ")
mesh <- list(points = mesh_points, simplices = simplices)

# Construct FEM and precision matrix
fem <- compute_fem_matrices(mesh$points, mesh$simplices)
Q <- assemble_precision_matrix(fem$C, fem$G)
check_matrix_sanity(Q)

# Simulate latent field
Y <- simulate_latent_field(Q)
Y <- Y - mean(Y)

# Likelihood data
bounds <- replicate(3, c(0, 1), simplify = FALSE)
likelihood_data <- construct_likelihood_data(mesh, coords_beak, covariate_fn, bounds)

# Run INLA model
result <- run_inla_continuous(likelihood_data, Q, estimate_beta = TRUE)

# Output
cat("Estimated beta:\n")
print(result$summary.fixed)

cat("\nLatent field RMSE:\n")
rmse <- sqrt(mean((Y - result$summary.random$idx$mean)^2))
print(rmse)

cat("\nCorrelation between true and estimated field:\n")
print(cor(Y, result$summary.random$idx$mean))


plot_field_slices <- function(mesh_points, field_values, fixed_dim = 3, n_slices = 6, output_dir = "slices") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  dims <- 1:3
  var_dims <- dims[-fixed_dim]
  
  # Safety: remove NA rows
  valid_rows <- complete.cases(mesh_points, field_values)
  mesh_points <- mesh_points[valid_rows, ]
  field_values <- field_values[valid_rows]
  
  fixed_vals <- mesh_points[, fixed_dim]
  if (all(is.na(fixed_vals))) stop("All fixed dimension values are NA")
  
  slice_vals <- seq(min(fixed_vals), max(fixed_vals), length.out = n_slices)
  
  for (i in seq_along(slice_vals)) {
    val <- slice_vals[i]
    tol <- (max(fixed_vals) - min(fixed_vals)) / (n_slices * 2)
    sel <- abs(fixed_vals - val) < tol
    
    if (sum(sel) < 10) {
      message(sprintf("Skipping slice %d: too few points", i))
      next
    }
    
    sub_points <- mesh_points[sel, ]
    sub_field <- field_values[sel]
    
    interp_result <- akima::interp(
      x = sub_points[, var_dims[1]],
      y = sub_points[, var_dims[2]],
      z = sub_field,
      linear = TRUE,
      extrap = FALSE
    )
    
    interp_df <- expand.grid(x = interp_result$x, y = interp_result$y)
    interp_df$z <- as.vector(interp_result$z)
    
    p <- ggplot(interp_df, aes(x = x, y = y, fill = z)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_viridis_c(option = "magma", na.value = "gray20") +
      coord_fixed() +
      theme_minimal() +
      labs(
        title = paste("Slice", i, "at dim", fixed_dim, "≈", round(val, 2)),
        fill = "Field"
      )
    
    ggsave(filename = sprintf("%s/slice_%02d.png", output_dir, i), plot = p, width = 6, height = 5)
  }
}

plot_field_slices(
  mesh_points = result$mesh_points,
  field_values = result$estimated_field,
  fixed_dim = 3,  
  n_slices = 10,
  output_dir = "slices"
)

