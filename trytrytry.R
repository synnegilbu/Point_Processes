library(MASS)
library(Matrix)
library(INLA)
library(geometry)

# --- LGCP Simulator (your function) ---
simulate_lgcp <- function(d, bounds, m, length_scale, covariate_field, covariate_coeff, seed = 123) {
  set.seed(seed)
  xlims <- bounds
  grid_sides <- lapply(xlims, function(lim) seq(lim[1], lim[2], length.out = m))
  grid_coords <- expand.grid(grid_sides)
  X <- as.matrix(grid_coords)
  grid_step <- sapply(xlims, function(lim) diff(lim) / m)
  
  # RBF kernel
  rbf <- function(X1, X2, length_scale) {
    dists <- as.matrix(dist(rbind(X1, X2)))
    exp(-dists[1:nrow(X1), (nrow(X1)+1):(nrow(X1)+nrow(X2))]^2 / (2 * length_scale^2))
  }
  
  K <- rbf(X, X, length_scale)
  # Comment out this line if you want to skip the true latent field:
  Y <- mvrnorm(mu = rep(0, nrow(X)), Sigma = K)
  
  covariate_values <- apply(X, 1, covariate_field)
  covariate_effect <- exp(covariate_coeff * covariate_values)
  
  Z <- exp(Y) * covariate_effect
  volume <- prod(sapply(xlims, diff))
  total_rate <- sum(Z) * (volume / (m^d))
  N <- rpois(1, total_rate)
  
  probabilities <- Z / sum(Z)
  sampled_indices <- sample(1:length(Z), size = N, replace = TRUE, prob = probabilities)
  
  sampled_points <- matrix(0, nrow = N, ncol = d)
  for (i in seq_len(d)) {
    cell_centers <- grid_coords[sampled_indices, i]
    offsets <- runif(N, min = -grid_step[i] / 2, max = grid_step[i] / 2)
    sampled_points[, i] <- cell_centers + offsets
  }
  
  list(
    sampled_points = sampled_points,
    rates = Z,
    grid_coords = X,
    bounds = bounds,
    m = m,
    Y = Y  
  )
}

# --- Grid binning ---
get_grid_counts <- function(sampled_points, bounds, m) {
  d <- length(bounds)
  grid_step <- sapply(bounds, function(lim) diff(lim) / m)
  if (nrow(sampled_points) == 0) return(rep(0, m^d))
  
  point_indices <- apply(sampled_points, 1, function(pt) {
    floor((pt - sapply(bounds, `[[`, 1)) / grid_step) + 1
  })
  
  point_indices <- t(point_indices)
  point_indices <- pmin(pmax(point_indices, 1), m)
  
  linear_indices <- as.vector(
    1 + point_indices %*% cumprod(c(1, rep(m, d - 1))) - 1
  )
  
  tabulate(linear_indices, nbins = m^d)
}

# --- Precision matrix for grid GMRF ---
build_Q_grid <- function(m, d) {
  total <- m^d
  Q <- Matrix(0, total, total, sparse = TRUE)
  
  get_neighbors <- function(idx) {
    coords <- arrayInd(idx, .dim = rep(m, d))
    neighbors <- list()
    for (i in 1:d) {
      for (delta in c(-1, 1)) {
        neighbor <- coords
        neighbor[i] <- neighbor[i] + delta
        if (all(neighbor >= 1 & neighbor <= m)) {
          neighbors <- c(neighbors, list(neighbor))
        }
      }
    }
    sapply(neighbors, function(coord) {
      sum((coord - 1) * cumprod(c(1, rep(m, d - 1)[1:(d - 1)]))) + 1
    })
  }
  
  for (i in 1:total) {
    neighbors <- get_neighbors(i)
    Q[i, neighbors] <- -1
    Q[i, i] <- length(neighbors)
  }
  
  Q + Diagonal(n = total, x = 1e-5)  # Add jitter for stability
}

# --- Run the full simulation + INLA inference ---
run_lgcp_inla <- function(d = 3, m = 10, covariate_coeff = 0.05) {
  bounds <- replicate(d, c(0, 10), simplify = FALSE)
  covariate_fn <- function(x) sum(x)  # example: linear spatial trend
  
  result <- simulate_lgcp(
    d = d, bounds = bounds, m = m,
    length_scale = 2,
    covariate_field = covariate_fn,
    covariate_coeff = covariate_coeff,
    seed = 123
  )
  
  grid_counts <- get_grid_counts(result$sampled_points, bounds, m)
  Q <- build_Q_grid(m, d)
  
  data_inla <- data.frame(
    y = grid_counts,
    idx = 1:length(grid_counts)
  )
  
  result_inla <- inla(
    y ~ f(idx, model = "generic0", Cmatrix = Q),
    family = "poisson",
    data = data_inla,
    control.predictor = list(compute = TRUE),
    control.inla = list(strategy = "gaussian")
  )
  
  list(
    simulation = result,
    inla_result = result_inla
  )
}

# --- Run it ---
out <- run_lgcp_inla(d = 3, m = 10)

# --- Plot INLA estimate vs true field ---
plot(out$simulation$Y, out$inla_result$summary.random$idx$mean,
     xlab = "True latent field (Y)",
     ylab = "INLA estimated field",
     pch = 16, col = rgb(0, 0, 1, 0.3),
     main = "INLA vs True Field")
abline(0, 1, col = "red")
