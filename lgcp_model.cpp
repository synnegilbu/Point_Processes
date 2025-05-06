#include <TMB.hpp>
using namespace density;  

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_MATRIX(grid_coords);
  DATA_VECTOR(counts);
  DATA_VECTOR(covariate_values);
  DATA_SCALAR(grid_volume);

  // Parameters
  PARAMETER(log_length_scale);
  PARAMETER(log_sigma);
  PARAMETER(beta);
  PARAMETER_VECTOR(field);

  int n = grid_coords.rows();
  Type length_scale = exp(log_length_scale);
  Type sigma = exp(log_sigma);

  Type nll = 0.0;

  // Build RBF covariance matrix
  matrix<Type> cov(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      Type dist2 = ((grid_coords.row(i) - grid_coords.row(j)).array().square()).sum();
      cov(i, j) = sigma * sigma * exp(-dist2 / (2 * length_scale * length_scale));
      if (i != j) cov(j, i) = cov(i, j);  // symmetry
    }
  }

  // MVN prior
  MVNORM_t<Type> mvn(cov);
  nll += mvn(field);

  // Poisson likelihood
  for (int i = 0; i < n; i++) {
    Type lambda = exp(field(i)) * exp(beta * covariate_values(i));
    nll -= dpois(counts(i), lambda * grid_volume, true);
  }

  return nll;
}
