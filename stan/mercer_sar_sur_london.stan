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
  int<lower=1> n; 
  int<lower=1> M;  
  int<lower=1> N;  
  int<lower=1> Kbeta;
  matrix[n,n] D2; 
  matrix[N,Kbeta] X;  
  vector[N] Y;
  real<lower=0> eps; // tiny value
  
  // priors
  real mu_ell;        
  real<lower=0> sd_ell;
  vector[Kbeta] u_beta;
  matrix[Kbeta,Kbeta] V_beta; // SPD
  vector[M] m_tau;    
  vector<lower=0>[M] s_tau;
  real<lower=0> eta_lkj; // LKJ shape
}
parameters {
  vector[Kbeta] beta;
  real log_ell; // log length-scale (km)
  vector[M] log_tau;     
  cholesky_factor_corr[M] L_corr; // Cholesky of correlation
}
transformed parameters {
  real<lower=0> ell = exp(log_ell);
  matrix[n,n] K = exp( - D2 / (ell*ell) );
  for (i in 1:n) K[i,i] = 1.0;
  K = K + diag_matrix(rep_vector(eps, n));
  
  vector<lower=0>[M] tau = exp(log_tau);  
  matrix[M,M] L_S = diag_pre_multiply(tau, L_corr);
  matrix[M,M] Sigma = multiply_lower_tri_self_transpose(L_S);
  matrix[M,M] Sigma_inv = inverse_spd(Sigma);

  matrix[N,N] Prec = kron_prod(K, Sigma_inv);
}
model {
  // priors
  beta ~ multi_normal(u_beta, V_beta);
  log_ell ~ normal(mu_ell, sd_ell);
  log_tau ~ normal(m_tau, s_tau);
  L_corr ~ lkj_corr_cholesky(eta_lkj);
  
  // likelihood
  Y ~ multi_normal_prec(X*beta, Prec);
}
generated quantities {
  real ell_km = exp(log_ell);
  corr_matrix[M] R = multiply_lower_tri_self_transpose(L_corr);
  real rho12 = (M>=2) ? R[1,2] : 0;
  vector[M] sigma = exp(log_tau);
  real r10 = 1.517 * ell_km; // sqrt(log 10)
}
