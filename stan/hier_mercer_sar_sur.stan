functions {
  matrix kron_prod(matrix A, matrix B) {
    int p = rows(A), q = rows(B);
    matrix[p*q, p*q] C;
    for (i in 1:p) for (j in 1:p)
      C[(i-1)*q+1:i*q, (j-1)*q+1:j*q] = A[i,j] * B;
    return C;
  }
}
data {
  int<lower=1> C; // counties
  int<lower=1> M; // equations (2 here)
  int<lower=1> Kbeta; // coeffs per county block
  int<lower=1> max_n; // max neighborhoods across counties
  int<lower=1> max_N; // = max_n * M
  int<lower=1> n[C];  // neighborhoods in county c
  int<lower=1> N[C];  // = n[c] * M
  matrix[max_n, max_n] D2[C]; // squared centroid distances (km^2)
  matrix[max_N, Kbeta] X[C];  // block-diagonal designs
  vector[max_N] Y[C];  // responses
  real<lower=0> eps; // kernel jitter
  vector[C] log_r0; // log median NN distance per region (km)
  real<lower=0> eta_lkj; // LKJ shape 
}
parameters {
  // Hierarchical betas 
  vector[Kbeta] mu_beta;   // population mean for betas
  vector[Kbeta] log_tau_beta;  // log SDs for population betas
  cholesky_factor_corr[Kbeta] Lcorr_beta;  // population correlation
  matrix[Kbeta, C] z_beta;  // county-specific standard normals

  // length-scale hierarchy with per-region centering
  real mu_ell;  // population shift from log r0
  real<lower=0> tau_ell;  // population SD on log-scale
  vector[C] z_ell; // county-specific standard normals

  // SUR covariance (shared across counties)
  vector[M]  log_tau;  // log SDs per equation
  cholesky_factor_corr[M] L_corr;  // across-equation correlation
}
transformed parameters {
  matrix[Kbeta, C] beta_c;
  vector[Kbeta] tau_beta = exp(log_tau_beta);
  vector[C] log_ell;
  vector<lower=0>[M] tau = exp(log_tau);

  {
    matrix[Kbeta, Kbeta] Lbeta = diag_pre_multiply(tau_beta, Lcorr_beta);
    for (c in 1:C)
      beta_c[, c] = mu_beta + Lbeta * z_beta[, c];
  }
  log_ell = log_r0 + mu_ell + tau_ell * z_ell;   // per-region centering
}
model { // priors matches to sinlge fit
  // Priors: coefficients 
  mu_beta ~ normal(0, 1);
  log_tau_beta ~ normal(0, 0.5);  
  Lcorr_beta ~ lkj_corr_cholesky(eta_lkj);
  to_vector(z_beta) ~ normal(0, 1);

  // Priors: log length-scale
  mu_ell ~ normal(0, 0.50);  // population shift from log r0
  tau_ell ~ normal(0, 0.30);  // half-Normal via <lower=0>
  z_ell ~ normal(0, 1);

  // Priors: SUR 
  log_tau ~ normal(log(0.5), 0.35);
  L_corr ~ lkj_corr_cholesky(eta_lkj);

  // Likelihood
  {
    matrix[M, M] Ls = diag_pre_multiply(tau, L_corr);
    matrix[M, M] Sig = multiply_lower_tri_self_transpose(Ls);
    matrix[M, M] Sig_inv = inverse_spd(Sig);

    for (c in 1:C) {
      int nc = n[c];
      int Nc = N[c];

      // Mercer kernel per county 
      matrix[nc, nc] Kc;
      {
        real ell = exp(log_ell[c]);
        Kc = exp( - D2[c][1:nc, 1:nc] / (ell * ell) );
        for (i in 1:nc) Kc[i, i] = 1.0;
        Kc = Kc + diag_matrix(rep_vector(eps, nc));
      }

      // Kronecker-structured precision
      {
        matrix[Nc, Nc] Prec = kron_prod(Sig_inv, Kc);
        vector[Nc] mu = X[c][1:Nc, 1:Kbeta] * beta_c[, c];
        target += multi_normal_prec_lpdf( Y[c][1:Nc] | mu, Prec );
      }
    }
  }
}
generated quantities {
  vector[C] ell_km;
  vector[C] r10_km;
  vector[M] sigma = exp(log_tau);
  corr_matrix[M] R = multiply_lower_tri_self_transpose(L_corr);
  real rho12 = (M >= 2) ? R[1, 2] : 0;

  for (c in 1:C) {
    ell_km[c] = exp(log_ell[c]);
    r10_km[c] = 1.517 * ell_km[c];    
  }
}