# Install required package if not already installed
if (!requireNamespace("spatstat.core", quietly = TRUE)) {
  install.packages("spatstat.core")
}
library(spatstat)
library(spatstat.geom)

# Define intensity (lambda) and window
lambda <- 100        # expected number of points per unit area
window <- owin(c(0, 1), c(0, 1))  # unit square window

# Simulate homogeneous Poisson point process
ppp_sim <- rpoispp(lambda, win = window)

# Plot the result
plot(ppp_sim, main = paste(""))
