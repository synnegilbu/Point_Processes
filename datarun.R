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
trait_vars <- c("Beak.Length_Culmen", "Beak.Depth", "Beak.Width")
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

# --- INLA Inference Function ---
run_inla_inference <- function(counts, covariate, Q) {
  df <- data.frame(y = counts, idx = 1:length(counts), covariate = covariate)
  fmla <- y ~ covariate + f(idx, model = "generic0", Cmatrix = Q,
                            hyper = list(prec = list(initial = 2, fixed = FALSE)))
  inla(fmla, family = "nbinomial", data = df,
       control.predictor = list(compute = TRUE),
       control.inla = list(strategy = "laplace"))
}

# --- LGCP Wrapper Function ---
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
    estimated_field = result$summary.random$idx$mean,
    fixed_effects = result$summary.fixed,
    result = result
  )
}

# --- Run Model ---
result <- run_lgcp(coords_scaled, covariate, m = 20)

# --- Visualize Density Surface ---
intensity_df <- data.frame(
  x = result$mesh_points[, 1],
  y = result$mesh_points[, 2],
  z = result$mesh_points[, 3],
  log_intensity = result$estimated_field,
  intensity = exp(result$estimated_field)
)



install.packages("gganimate")
install.packages("gifski")
library(gganimate)
library(gifski)
library(dplyr)
library(ggplot2)


z_slices <- seq(0, 1, length.out = 100)
slice_thickness <- 0.001

sliced_df <- do.call(rbind, lapply(z_slices, function(zval) {
  slice <- intensity_df %>%
    filter(abs(z - zval) <= slice_thickness) %>%
    mutate(z_slice = round(zval, 2))
  return(slice)
}))


p <- ggplot(sliced_df, aes(x = x, y = y, fill = intensity)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log", name = "Density") +
  labs(
    title = 'Species Density Slice | Beak Width (z) ≈ {closest_state}',
    x = "Beak Length", y = "Beak Depth"
  ) +
  theme_minimal() +
  transition_states(z_slice, transition_length = 2, state_length = 1) +
  ease_aes('linear')

animate(p, renderer = gifski_renderer("trait_space_slices.gif"), width = 800, height = 600, fps = 5)




