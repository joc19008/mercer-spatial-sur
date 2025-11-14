# Mercer SAR–SUR simulation
# Following Lopez's sim. structure - two eqns, intercept+one regressor.

set.seed(2025)
library(Matrix)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)


# reparameterization
n  <- 124
M <- 2
ell_true_km  <- 10              
Sigma_true <- diag(M)         
beta1_true <- c(1,  1)
beta2_true <- c(1, -1)
eps_nugget <- 1e-6

# Coordinates in km: 0..25 km 
coords_km <- cbind(runif(n, 0, 25), runif(n, 0, 25))   # n x 2, km
D_km <- as.matrix(dist(coords_km))
D2_km2 <- D_km^2

# Characteristic length r0 = median nearest-neighbor distance (km) 
nn_d <- apply(replace(D_km, row(D_km)==col(D_km), Inf), 1, min)
r0 <- median(nn_d)

# Mercer kernel on precision side with chosen ell_true_km 
K <- exp( - D2_km2 / (ell_true_km^2) )
diag(K) <- diag(K) + eps_nugget

# Factors: A A' = K^{-1}, B'B = Sigma 
R_K <- chol(K) 
A <- backsolve(R_K, diag(n))
R_S <- chol(Sigma_true)
B <- R_S

#  regressors and means 
x1 <- runif(n); x2 <- runif(n)
X1 <- cbind(1, x1); X2 <- cbind(1, x2)
mu1 <- as.vector(X1 %*% beta1_true)
mu2 <- as.vector(X2 %*% beta2_true)

# Errors with Cov(vec(E)) 
G <- matrix(rnorm(n*2L), n, 2L)
E <- A %*% G %*% B

# Responses
y1 <- mu1 + E[,1]
y2 <- mu2 + E[,2]

# Stan data
X_blk <- as.matrix(bdiag(X1, X2))         # (2n) × (k1+k2)
perm <- as.vector(rbind(1:n, (n+1):(2*n)))  
X_loc <- X_blk[perm, , drop = FALSE]
Y_loc <- c(rbind(y1, y2))

N <- n*M
Kbeta <- ncol(X_loc)


# Mercer SAR–SUR (precision form) in Stan
stan_code_mercer <- '
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
  // real<lower=0> r0; // prior center for ell (km)
  real<lower=0> eps; // tiny jitter
  
  // priors
  real mu_ell;        real<lower=0> sd_ell;
  vector[Kbeta] u_beta;
  matrix[Kbeta,Kbeta] V_beta; // SPD
  vector[M] m_tau;    vector<lower=0>[M] s_tau;
  real<lower=0> eta_lkj; // LKJ shape, e.g., 2
}
parameters {
  vector[Kbeta] beta;
  real log_ell; // log length-scale (km)
  // vector<lower=0>[M] tau; // per-equation SDs
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
  vector[M] tau = exp(log_tau);
  real r10 = 1.517 * ell_km; // sqrt(log 10)
}
'

mod_mercer <- stan_model(model_code = stan_code_mercer)


stan_data_mercer <- list(
  n = n, M = M, N = N, Kbeta = Kbeta,
  D2 = D2_km2, X = X_loc, Y = Y_loc,
  eps = 1e-6,
  mu_ell = log(r0), sd_ell = 0.5,
  u_beta = rep(0, Kbeta),
  V_beta = diag(1, Kbeta),
  m_tau = rep(log(0.5), M),
  s_tau = rep(0.35, M),
  eta_lkj = 2.0
)


fit.sim <- sampling(
  mod_mercer, data = stan_data_mercer,
  chains = 2, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

print(fit.sim, 
      pars = c("beta", "Sigma", "ell_km", "rho12", "sigma", "r10"),
      probs = c(0.025,0.5,0.975))

# saveRDS(fit.sim, "Output/251024_sim_mercer_fit.Rds")

fit.sim <- readRDS("Output/251024_sim_mercer_fit.Rds")









