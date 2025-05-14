library(geometry)
library(FNN)
library(Matrix)
library(INLA)
library(scales)  #

data <- read.csv("/Users/Synne/Documents/Skole/H2024/master/docs/datasets/doi_10_5061_dryad_zgmsbcckk__v20250122/DB_Avonet_BirdLife.csv")
colnames(data) <- trimws(colnames(data))
data <- data[!is.na(data$beak_length_culmen) & 
               !is.na(data$wing_length) & 
               !is.na(data$body_mass), 
             c("species", "beak_length_culmen", "wing_length", "body_mass")]

# Scale all three traits
coords_scaled <- apply(data[, c("beak_length_culmen", "wing_length", "body_mass")], 2, function(x) rescale(x, to = c(0, 1)))
species_names <- data$species


# --- Mesh Construction ---
build_mesh <- function(d, m, bounds) {
  grid_axes <- lapply(bounds, function(lim) seq(lim[1], lim[2], length.out = m))
  grid_points <- expand.grid(grid_axes)
  pts <- as.matrix(grid_points)
  simplices <- delaunayn(pts, options = "QJ")
  list(points = pts, simplices = simplices)
}

# --- FEM Matrices ---
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
    GK <- vol * (t(grads) %*% grads)
    MK <- matrix(1, d+1, d+1)
    diag(MK) <- 2
    MK <- MK * vol / ((d + 1) * (d + 2))
    C[idx, idx] <- C[idx, idx] + MK
    G[idx, idx] <- G[idx, idx] + GK
  }
  list(C = C, G = G)
}

# --- Precision Matrix ---
assemble_precision_matrix <- function(C, G, tau = 0.1, kappa = 2, jitter = 1e-5) {
  Q <- tau^2 * (kappa^2 * C + G)
  Q + Diagonal(nrow(Q), jitter)
}

# --- INLA Inference ---
run_inla_inference <- function(counts, covariate, Q, estimate_beta = TRUE) {
  df <- data.frame(y = counts, idx = 1:length(counts), cov = covariate)
  fmla <- if (estimate_beta) {
    y ~ cov + f(idx, model = "generic0", Cmatrix = Q,
                hyper = list(prec = list(initial = 2, fixed = FALSE)))
  } else {
    y ~ f(idx, model = "generic0", Cmatrix = Q,
          hyper = list(prec = list(initial = 2, fixed = FALSE)))
  }
  inla(fmla,
       family = "nbinomial",
       data = df,
       control.predictor = list(compute = TRUE),
       control.inla = list(strategy = "laplace"))
}

# --- Real Data LGCP Pipeline ---
run_lgcp_on_real_data <- function(coords, covariate_fn = function(x) x[1],
                                  beta = 0.1, estimate_beta = TRUE,
                                  scale_intensity = 2000, m = 10) {
  d <- ncol(coords)
  stopifnot(d >= 1)
  
  bounds <- lapply(1:d, function(i) range(coords[, i]))
  
  mesh <- build_mesh(d = d, m = m, bounds = bounds)
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  Q <- assemble_precision_matrix(fem$C, fem$G)
  
  covariate_vals <- apply(mesh$points, 1, covariate_fn)
  
  nn <- get.knnx(mesh$points, coords, k = 1)
  counts <- tabulate(nn$nn.index, nbins = nrow(mesh$points))
  
  cat("Observed species:", nrow(coords), "\n")
  cat("Non-zero bins:", sum(counts > 0), "/", length(counts), "\n")
  
  result <- run_inla_inference(counts, covariate_vals, Q, estimate_beta)
  
  list(
    mesh_points = mesh$points,
    counts = counts,
    estimated_field = result$summary.random$idx$mean,
    beta_estimate = if (estimate_beta) result$summary.fixed else NULL,
    result = result
  )
}


# --- Run Inference ---
res <- run_lgcp_on_real_data(
  coords = coords_scaled,
  covariate_fn = function(x) x[3],  # use body mass as covariate
  beta = 0.3,
  estimate_beta = TRUE,
  m = 15  # smaller mesh size for higher dimensions
)

library(plotly)

plot_ly(
  x = coords_scaled[, 1],
  y = coords_scaled[, 2],
  z = coords_scaled[, 3],
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2),
  text = species_names
)


library(ggplot2)

# Create a data frame
df <- data.frame(
  x = res$mesh_points[, 1],
  y = res$mesh_points[, 2],
  z = res$mesh_points[, 3],
  latent = res$estimated_field
)

# Slice at a mid-range z level (e.g., ~0.5)
slice_df <- subset(df, abs(res$mesh_points[, 3] - 0.5) < 0.02)

ggplot(slice_df, aes(x = x, y = y, fill = latent)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Latent Field Slice (Body Mass ≈ 0.5)",
       x = "Beak Length", y = "Wing Length") +
  theme_minimal()



# --- Plotting ---
hist(res$estimated_field, breaks = 50, main = "Estimated Latent Field", xlab = "Y_hat", col = "skyblue")
if (!is.null(res$beta_estimate)) {
  print(res$beta_estimate)
}


# Add required library
library(ggplot2)

# Data for heatmap
plot_df <- data.frame(
  x = res$mesh_points[, 1],
  y = res$mesh_points[, 2],
  latent = res$estimated_field
)

# Data for points
points_df <- data.frame(
  x = coords_scaled[, 1],
  y = coords_scaled[, 2],
  species = species_names
)

# Plot
ggplot() +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = latent)) +
  scale_fill_viridis_c(option = "plasma") +
  #geom_point(data = points_df, aes(x = x, y = y), color = "black", alpha = 0.5, size = 1) +
  labs(
    title = "Estimated Latent Field with Species Points",
    x = "Scaled Beak Length (culmen)",
    y = "Scaled Wing Length",
    fill = "Latent\nField"
  ) +
  theme_minimal()

