# --- Load Libraries ---
library(readxl)
library(INLA)
library(FNN)
library(Matrix)
library(geometry)
library(scales)
library(ggplot2)
library(dplyr)

# --- Load Data ---
data <- read_excel("/Users/Synne/Documents/Skole/H2024/master/code/Point_Processes/datasets/16586228/AVONET Supplementary dataset 1.xlsx", sheet = "AVONET1_BirdLife")

# --- Choose Traits and Covariate ---
trait_vars <- c("Beak.Length_Culmen", "Beak.Depth")  # 2D version
covariate_var <- "Trophic.Level"

# --- Clean and Scale Data ---
model_data <- na.omit(data[, c("Species1", covariate_var, trait_vars)])
coords_scaled <- apply(model_data[, trait_vars], 2, function(x) rescale(as.numeric(x), to = c(0, 1)))
covariate <- as.factor(model_data[[covariate_var]])
species_names <- model_data$Species1

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
run_inla_inference <- function(counts, covariate, Q) {
  df <- data.frame(y = counts, idx = 1:length(counts), covariate = covariate)
  fmla <- y ~ covariate + f(idx, model = "generic0", Cmatrix = Q,
                            hyper = list(prec = list(initial = 2, fixed = FALSE)))
  inla(fmla, family = "nbinomial", data = df,
       control.predictor = list(compute = TRUE),
       control.inla = list(strategy = "laplace"))
}

# --- LGCP Wrapper ---
run_lgcp <- function(coords, covariate, m = 20) {
  d <- ncol(coords)
  bounds <- lapply(1:d, function(i) range(coords[, i]))
  mesh <- build_mesh(d, m, bounds)
  fem <- compute_fem_matrices(mesh$points, mesh$simplices)
  Q <- assemble_precision_matrix(fem$C, fem$G)
  
  nn <- get.knnx(mesh$points, coords, k = 1)
  counts <- tabulate(nn$nn.index, nbins = nrow(mesh$points))
  
  cov_on_mesh <- sapply(1:nrow(mesh$points), function(i) {
    idx <- which.min(colSums((t(coords) - mesh$points[i, ])^2))
    as.character(covariate[idx])
  })
  cov_on_mesh <- as.factor(cov_on_mesh)
  
  result <- run_inla_inference(counts, cov_on_mesh, Q)
  
  list(
    mesh_points = mesh$points,
    simplices = mesh$simplices,
    estimated_field = result$summary.random$idx$mean,
    fixed_effects = result$summary.fixed,
    result = result
  )
}

# --- Run Model ---
result <- run_lgcp(coords_scaled, covariate, m = 20)

# ----------------------------------------------------------------
# --- PLOT 1: Latent Field Visualization (heatmap)
# ----------------------------------------------------------------
field_df <- data.frame(
  x = result$mesh_points[, 1],
  y = result$mesh_points[, 2],
  latent = result$estimated_field
)

p_field <- ggplot(field_df, aes(x = x, y = y, fill = latent)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Estimated Field") +
  labs(title = "Estimated Latent Field", x = "Beak Length", y = "Beak Depth") +
  coord_fixed() +
  theme_minimal()

print(p_field)

# ----------------------------------------------------------------
# --- PLOT 2: Mesh Visualization
# ----------------------------------------------------------------
simplices <- result$simplices
mesh_points <- result$mesh_points

# Build edge list
edges_list <- lapply(1:nrow(simplices), function(i) {
  simplex <- simplices[i, ]
  rbind(
    c(simplex[1], simplex[2]),
    c(simplex[2], simplex[3]),
    c(simplex[3], simplex[1])
  )
})
edges_mat <- do.call(rbind, edges_list)
edges_df <- data.frame(
  x1 = mesh_points[edges_mat[, 1], 1],
  y1 = mesh_points[edges_mat[, 1], 2],
  x2 = mesh_points[edges_mat[, 2], 1],
  y2 = mesh_points[edges_mat[, 2], 2]
)

p_mesh <- ggplot() +
  geom_segment(data = edges_df, aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "grey40", size = 0.3) +
  labs(title = "FEM Mesh", x = "Beak Length", y = "Beak Depth") +
  coord_fixed() +
  theme_minimal()

print(p_mesh)


# Load 3D library
library(rgl)

# Extract mesh points and field
mesh_points <- result$mesh_points
estimated_field <- result$estimated_field

# 3D scatter plot of latent field
open3d()
plot3d(mesh_points[,1], mesh_points[,2], estimated_field,
       col = "lightblue", size = 5, type = "p",
       xlab = "Beak Length (scaled)", 
       ylab = "Beak Depth (scaled)", 
       zlab = "Estimated Latent Field")

# Proper FEM surface: use triangles from mesh
open3d()
for (i in 1:nrow(simplices)) {
  simplex <- simplices[i, ]
  x_coords <- mesh_points[simplex, 1]
  y_coords <- mesh_points[simplex, 2]
  z_coords <- estimated_field[simplex]
  
  triangles3d(x_coords, y_coords, z_coords, color = "skyblue", alpha = 0.9)
}

library(rgl)

mesh_points <- result$mesh_points
simplices <- result$simplices
estimated_field <- result$estimated_field

# Clip for visualization
clip_limit <- 100
clipped_field <- pmax(pmin(estimated_field, clip_limit), -clip_limit)

# BIGGER scale for visualization
scale_factor <- 0.05
z_field <- scale_factor * clipped_field

# 3D plot
open3d()
for (i in 1:nrow(simplices)) {
  simplex <- simplices[i, ]
  x_coords <- mesh_points[simplex, 1]
  y_coords <- mesh_points[simplex, 2]
  z_coords <- z_field[simplex]
  
  triangles3d(x_coords, y_coords, z_coords, color = "lightblue", alpha = 0.9)
}

# Plot mesh grid underneath
z_offset <- min(z_field) - 0.5
for (i in 1:nrow(simplices)) {
  simplex <- simplices[i, ]
  x_coords <- mesh_points[simplex, 1]
  y_coords <- mesh_points[simplex, 2]
  z_coords <- rep(z_offset, 3)
  
  triangles3d(x_coords, y_coords, z_coords, color = "grey30", alpha = 1, front = "lines")
}



library(rgl)

mesh_points <- result$mesh_points
simplices <- result$simplices
estimated_field <- result$estimated_field

# Clip field values for better visualization
clip_limit <- 100
clipped_field <- pmax(pmin(estimated_field, clip_limit), -clip_limit)

# Scale z-values for reasonable height
scale_factor <- 0.05
z_field <- scale_factor * clipped_field

open3d()

# First: plot surface (solid)
for (i in 1:nrow(simplices)) {
  simplex <- simplices[i, ]
  x_coords <- mesh_points[simplex, 1]
  y_coords <- mesh_points[simplex, 2]
  z_coords <- z_field[simplex]
  
  triangles3d(x_coords, y_coords, z_coords, color = "lightblue", alpha = 0.9)
}

# Now: plot the mesh "rising" with the surface
for (i in 1:nrow(simplices)) {
  simplex <- simplices[i, ]
  x_coords <- mesh_points[simplex, 1]
  y_coords <- mesh_points[simplex, 2]
  z_coords <- z_field[simplex]
  
  lines3d(c(x_coords, x_coords[1]),
          c(y_coords, y_coords[1]),
          c(z_coords, z_coords[1]),
          color = "black", lwd = 1)
}

aspect3d(1, 1, 1.5)

