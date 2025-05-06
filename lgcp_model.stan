data {
    int<lower=1> N;              // Number of mesh nodes
    int<lower=1> M;              // Number of observed points
    matrix[N, N] Q;              // SPDE precision matrix (precision of latent field)
    int<lower=0> Y[M];           // Observed count data (Poisson distributed)
    int<lower=1, upper=N> mesh_index[M];  // Index mapping observations to mesh nodes
    real log_tau_prior_mean;     // Prior mean for log(tau)
    real log_kappa_prior_mean;   // Prior mean for log(kappa)
}

parameters {
    vector[N] Z;                // Latent Gaussian field
    real<lower=0> sigma;        // Variance parameter for GRF
    real log_tau;               // SPDE parameter (log-transformed)
    real log_kappa;             // SPDE parameter (log-transformed)
}

transformed parameters {
    real tau = exp(log_tau);
    real kappa = exp(log_kappa);
}

model {
    // Priors for SPDE parameters based on domain knowledge
    log_tau ~ normal(log_tau_prior_mean, 0.5);
    log_kappa ~ normal(log_kappa_prior_mean, 0.5);

    // SPDE prior: Gaussian field with precision matrix Q
    Z ~ multi_normal_prec(rep_vector(0, N), Q);

    // Log-Gaussian Cox Process likelihood
    for (m in 1:M) {
        Y[m] ~ poisson_log(Z[mesh_index[m]]);
    }

    // Prior for sigma
    sigma ~ normal(0, 1);
}

generated quantities {
    vector[N] lambda_pred;
    for (n in 1:N) {
        lambda_pred[n] = exp(Z[n]);  // Convert log-Gaussian field to intensity scale
    }
}
