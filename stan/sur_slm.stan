functions {
  real log_det_I_minus_rhoW(matrix W, real rho) {
    int n = rows(W);
    matrix[n,n] A = diag_matrix(rep_vector(1.0, n)) - rho * W;
    return log_determinant(A);
  }
}
data {
  int<lower=1> n;
  int<lower=1> K1;
  int<lower=1> K2;
  matrix[n,K1] X1;
  matrix[n,K2] X2;
  vector[n] y1;
  vector[n] y2;
  matrix[n,n] W;     // row-standardized
  real rho_lower;
  real rho_upper;
  matrix[2,2] Omega;
  real<lower=2> nu;
}
parameters {
  vector[K1] beta1;
  vector[K2] beta2;
  cov_matrix[2] Sinv;  // precision of SUR errors
  real<lower=0,upper=1> rho_raw1;
  real<lower=0,upper=1> rho_raw2;
}
transformed parameters {
  real rho1 = rho_lower + (rho_upper - rho_lower) * rho_raw1;
  real rho2 = rho_lower + (rho_upper - rho_lower) * rho_raw2;
}
model {
  matrix[n,n] A1 = diag_matrix(rep_vector(1.0, n)) - rho1 * W;
  matrix[n,n] A2 = diag_matrix(rep_vector(1.0, n)) - rho2 * W;

  vector[n] r1 = A1 * y1 - X1 * beta1;
  vector[n] r2 = A2 * y2 - X2 * beta2;

  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  Sinv ~ wishart(nu, Omega);

  target += log_det_I_minus_rhoW(W, rho1) + log_det_I_minus_rhoW(W, rho2);

  {
    matrix[2,2] S = inverse_spd(Sinv);
    for (i in 1:n) {
      vector[2] e_i;
      e_i[1] = r1[i];
      e_i[2] = r2[i];
      target += multi_normal_lpdf(e_i | rep_vector(0.0,2), S);
    }
  }
}
generated quantities {
  real log_det1 = log_det_I_minus_rhoW(W, rho1);
  real log_det2 = log_det_I_minus_rhoW(W, rho2);
}