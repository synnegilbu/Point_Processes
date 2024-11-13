library(spatstat)
library(spatstat.geom)
library(ggplot2)
library(plotly)

# Function to simulate a Log-Gaussian Cox Process in 2D or 3D
generalized_lgcp <- function(d, n, sigma2, range) {
  if(d == 2) {
    # 2D Case: Simulate a 2D Gaussian Random Field
    sim2 <- grf(n, grid = "reg", cov.pars = c(sigma2, range))
    
    # Extract and process GRF values
    grf_values <- sim2$data
    grid_size <- sqrt(length(grf_values))  # Calculate grid size
    
    # Reshape GRF values into a 2D matrix
    grf_matrix <- matrix(grf_values, nrow = grid_size, ncol = grid_size)
    
    # Exponentiate to obtain inhomogeneous intensity
    intensity <- exp(grf_matrix)
    lambda <- intensity / sum(intensity) * 500  # Scale to control total points
    
    # Define grid points in 2D space
    x_coords <- seq(0, 1, length.out = grid_size)
    y_coords <- seq(0, 1, length.out = grid_size)
    
    # Generate points based on Poisson intensities in each cell
    points_x <- c()
    points_y <- c()
    
    for (i in 1:grid_size) {
      for (j in 1:grid_size) {
        # Number of points to place in this cell
        num_points <- rpois(1, lambda[i, j])
        
        if (num_points > 0) {
          # Randomly distribute points within the cell
          points_x <- c(points_x, runif(num_points, x_coords[i] - 1/(2*grid_size), x_coords[i] + 1/(2*grid_size)))
          points_y <- c(points_y, runif(num_points, y_coords[j] - 1/(2*grid_size), y_coords[j] + 1/(2*grid_size)))
        }
      }
    }
    
    # Return the points as a data frame for 2D
    pp_2d <- data.frame(x = points_x, y = points_y)
    return(pp_2d)
    
  } else if(d == 3) {
    # 3D Case: Setup grid dimensions
    dim_size <- round(n^(1/3))  # Cube dimension size
    sim3 <- grf(dim_size^3, grid = "reg", cov.pars = c(sigma2, range))
    
    # Extract and reshape GRF values into 3D array
    grf_values_3d <- sim3$data
    grf_array <- array(grf_values_3d, dim = c(dim_size, dim_size, dim_size))
    
    # Exponentiate to obtain inhomogeneous intensity
    intensity_3d <- exp(grf_array)
    lambda_3d <- intensity_3d / sum(intensity_3d) * 500  # Scale to control total points
    
    # Define grid points in 3D space
    x_coords <- seq(0, 1, length.out = dim_size)
    y_coords <- seq(0, 1, length.out = dim_size)
    z_coords <- seq(0, 1, length.out = dim_size)
    
    # Generate points based on Poisson intensities in each cell
    points_x <- c()
    points_y <- c()
    points_z <- c()
    
    for (i in 1:dim_size) {
      for (j in 1:dim_size) {
        for (k in 1:dim_size) {
          # Number of points to place in this cell
          num_points <- rpois(1, lambda_3d[i, j, k])
          
          if (num_points > 0) {
            # Randomly distribute points within the voxel
            points_x <- c(points_x, runif(num_points, x_coords[i] - 1/(2*dim_size), x_coords[i] + 1/(2*dim_size)))
            points_y <- c(points_y, runif(num_points, y_coords[j] - 1/(2*dim_size), y_coords[j] + 1/(2*dim_size)))
            points_z <- c(points_z, runif(num_points, z_coords[k] - 1/(2*dim_size), z_coords[k] + 1/(2*dim_size)))
          }
        }
      }
    }
    
    # Return the points as a data frame for 3D
    pp_3d <- data.frame(x = points_x, y = points_y, z = points_z)
    return(pp_3d)
  } else {
    stop("Dimension must be either 2 or 3.")
  }
}

# Example usage for 2D
result_2d <- generalized_lgcp(d = 2, n = 1000, sigma2 = 1, range = 0.1)
ggplot(result_2d, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, color = "blue") +
  labs(title = "2D Log-Gaussian Cox Process", x = "X", y = "Y") +
  theme_minimal()

# Example usage for 3D
result_3d <- generalized_lgcp(d = 3, n = 1000, sigma2 = 1, range = 0.1)
plot_ly(data = result_3d, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 2, color = 'blue', opacity = 0.5)) %>%
  layout(
    title = "3D Log-Gaussian Cox Process",
    scene = list(
      xaxis = list(title = 'X'),
      yaxis = list(title = 'Y'),
      zaxis = list(title = 'Z')
    )
  )
