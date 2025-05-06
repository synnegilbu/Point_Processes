data {
  int<lower=1> N_nodes;          // Number of mesh nodes
  int<lower=1> N_obs;            // Number of observed points
  int<lower=1, upper=N_nodes> obs_idx[N_obs];  // Indices of mesh points for observed locations
  matrix[N_nodes, N_nodes] C;    // Mass matrix (symmetric)
  matrix[N_nodes, N_nodes] Q;    // Precision matrix for Z
}

parameters {
  vector[N_nodes] Z;             // Latent log-intensity field
}

model {
  // Prior: SPDE-based GP
  Z ~ multi_normal_prec(rep_vector(0, N_nodes), Q);

  // Likelihood: Approximate Poisson likelihood using mass matrix
  target += -dot_product(exp(Z), C * rep_vector(1, N_nodes));  // Integral term
  target += sum(Z[obs_idx]);  // Point-wise log intensity at observed locations
}
