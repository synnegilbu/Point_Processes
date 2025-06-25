# --- Load Libraries ---
library(readxl)
library(INLA)
library(FNN)
library(Matrix)
library(geometry)
library(scales)
library(ggplot2)

# --- Load Data from Excel ---
data <- read_excel("/Users/Synne/Documents/Skole/H2024/master/code/Point_Processes/datasets/16586228/AVONET Supplementary dataset 1.xlsx", sheet = "AVONET1_BirdLife")

# --- Choose Traits and Covariate ---
trait_vars_4d <- c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width", "Beak.Length_Nares")
covariate_var_4d <- "Trophic.Level"


# --- Clean and Scale Data ---
model_data_4d <- na.omit(data[, c("Species1", covariate_var_4d, trait_vars_4d)])
coords_scaled_4d <- apply(model_data_4d[, trait_vars_4d], 2, function(x) rescale(as.numeric(x), to = c(0, 1)))
covariate_4d <- as.factor(model_data_4d[[covariate_var_4d]])
species_names_4d <- model_data_4d$Species1


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

# --- INLA Inference Function ---
run_inla_inference <- function(counts, covariate, Q) {
  df <- data.frame(y = counts, idx = 1:length(counts), covariate = covariate)
  fmla <- y ~ covariate + f(idx, model = "generic0", Cmatrix = Q,
                            hyper = list(prec = list(param = c(10, 0.1))))
  inla(fmla, family = "nbinomial", data = df,
       control.predictor = list(compute = TRUE),
       control.inla = list(
         strategy = "laplace",
         diagonal = 1e-3,
         tolerance = 1e-5
       ))
}



# --- LGCP Wrapper Function ---
run_lgcp <- function(coords, covariate, m = 5) {
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
    estimated_field = result$summary.random$idx$mean,
    fixed_effects = result$summary.fixed,
    result = result
  )
}

# --- Run Model ---
result_4d <- run_lgcp(coords_scaled_4d, covariate_4d, m = 5)


result$fixed_effects

result$estimated_field

summary(result$result)


# --- Create Folder ---
dir.create("fig_4d_slices", showWarnings = FALSE)

# --- Setup ---
mesh_points <- result_4d$mesh_points
field <- result_4d$estimated_field

# --- Choose 10 evenly spaced values for dimension 3 (Beak.Width) and fix dim 4 ---
dim3_vals <- seq(0.05, 0.95, length.out = 10)
dim4_fixed <- 0.5  # Hold Beak.Length_Nares constant

# --- Generate and Save Slices ---
for (i in seq_along(dim3_vals)) {
  v3 <- dim3_vals[i]
  
  # Select points close to v3 in dimension 3 and fixed in dim 4
  indices <- which(abs(mesh_points[,3] - v3) < 0.1 & abs(mesh_points[,4] - dim4_fixed) < 0.1)
  if (length(indices) < 10) next
  
  slice_df <- data.frame(
    x = mesh_points[indices, 1],  # Trait 1
    y = mesh_points[indices, 2],  # Trait 2
    value = field[indices]
  )
  
  p <- ggplot(slice_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = paste0("Slice ", i, ": dim3 = ", round(v3, 2)),
         x = "Trait 1", y = "Trait 2") +
    theme_minimal()
  
  ggsave(
    filename = sprintf("fig_4d_slices/slice_%02d.png", i),
    plot = p, width = 5, height = 4, dpi = 300
  )
  cat(sprintf("Slice %02d: found %d points\n", i, length(indices)))
  
}





# 3D code

trait_vars_3D <- c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width")
model_data_3D <- na.omit(data[, c("Species1", covariate_var, trait_vars_3D)])
coords_scaled_3D <- apply(model_data_3D[, trait_vars_3D], 2, function(x) rescale(as.numeric(x), to = c(0, 1)))
covariate_3D <- as.factor(model_data_3D[[covariate_var]])
species_names_3D <- model_data_3D$Species1


# Run the LGCP model in 3D space
result_3D <- run_lgcp(coords_scaled_3D, covariate_3D, m = 10)  # You can experiment with m (e.g. 8, 12)

# Fixed effects
fixed_3d <- result_3D$fixed_effects
print(fixed_3d)

# Hyperparameters
hyper_3d <- result_3D$result$summary.hyperpar
print(hyper_3d)

library(ggplot2)
library(dplyr)
library(viridis)

# Extract mesh and field
mesh_points <- result_3D$mesh_points
field_values <- result_3D$estimated_field

# Ensure the folder exists
output_dir <- "fig_3d_slices"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message("Created directory: ", output_dir)
}

# Determine slice levels
z_vals <- seq(min(mesh_points[, 3]), max(mesh_points[, 3]), length.out = 10)
tol <- (max(mesh_points[, 3]) - min(mesh_points[, 3])) / 20  # Tolerance for slicing

# Loop through slices
for (i in seq_along(z_vals)) {
  z_fixed <- z_vals[i]
  slice_indices <- which(abs(mesh_points[, 3] - z_fixed) < tol)
  
  # Skip empty slices
  if (length(slice_indices) < 5) {
    message(sprintf("Skipping slice %d: too few points (%d)", i, length(slice_indices)))
    next
  }
  
  slice_data <- data.frame(
    x = mesh_points[slice_indices, 1],
    y = mesh_points[slice_indices, 2],
    z = mesh_points[slice_indices, 3],
    field = field_values[slice_indices]
  )
  
  p <- ggplot(slice_data, aes(x = x, y = y, fill = field)) +
    geom_tile() +
    scale_fill_viridis(option = "magma") +
    coord_fixed() +
    labs(
      title = paste("Slice", i, "at z ≈", round(z_fixed, 2)),
      x = "Trait 1 (scaled)", y = "Trait 2 (scaled)", fill = "Field"
    ) +
    theme_minimal()
  
  file_path <- file.path(output_dir, sprintf("slice_%02d.png", i))
  ggsave(filename = file_path, plot = p, width = 6, height = 5, dpi = 300)
  
  message("Saved: ", file_path)
}
