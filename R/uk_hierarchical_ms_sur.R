# Hierarchical Mercer SAR–SUR across 4 UK regions
# (London, Bristol, Edinburgh, Greater Manchester)

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

library(data.table)
library(dplyr)
library(FNN)
library(Matrix)
library(rstan)
library(loo)
library(mvtnorm)
library(ggplot2)
library(bayesplot)
library(sf)
library(viridis)
library(purrr)
library(tidyr)

sf::sf_use_s2(FALSE)


CRS_METRIC <- 27700  # OSGB 1936 / British National Grid (meters)

region_ids <- c("london", "bristol", "edinburgh", "manchester")
listings_path_fun <- function(id) sprintf("Data/listings_%s.csv", id)
neigh_path_fun <- function(id) sprintf("Data/neighbourhoods_%s.geojson", id)

form_X1 <- ~ 1 + med_host_cnt + prop_super + med_services + med_accomm +
  med_beds_home + med_baths_home + med_min_n + med_rev_m
form_X2 <- ~ 1 + med_host_cnt + prop_super + med_services + med_accomm +
  med_baths_room + med_min_n + med_rev_m 

# Replace NA medians in covariates by region-wise medians 
fill_na_with_med <- function(df, cols) {
  for (cl in cols) {
    if (cl %in% names(df)) {
      med <- suppressWarnings(median(df[[cl]], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      df[[cl]][is.na(df[[cl]])] <- med
    }
  }
  df
}

# Read and aggregate one region
process_region <- function(reg_id) {
  listings_path <- listings_path_fun(reg_id)
  neigh_path <- neigh_path_fun(reg_id)
  
  listings <- fread(listings_path, na.strings = c("", "NA")) %>%
    mutate(
      price = as.numeric(gsub("[$,]", "", price)),
      log_price = ifelse(!is.na(price), log(price), NA_real_),
      is_home = room_type == "Entire home/apt",
      is_room = room_type == "Private room",
      superhost = host_is_superhost == "t",
      n_amenities = lengths(strsplit(gsub("^\\[|\\]$", "", amenities), ",")),
      reviews_m = as.numeric(reviews_per_month),
      bathrooms = as.numeric(bathrooms),
      bedrooms = as.integer(bedrooms),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      minimum_nights= as.integer(minimum_nights),
      accommodates = as.integer(accommodates)
    ) %>%
    filter(room_type %in% c("Entire home/apt","Private room")) %>%
    filter(!is.na(price), !is.na(latitude), !is.na(longitude))
  
  neigh_raw <- st_read(neigh_path, quiet = TRUE) %>% st_transform(4326)
  nm <- tolower(names(neigh_raw))
  name_col <- if ("neighbourhood" %in% nm) names(neigh_raw)[which(nm=="neighbourhood")[1]] else
    if ("name" %in% nm) names(neigh_raw)[which(nm=="name")[1]] else names(neigh_raw)[1]
  neigh <- neigh_raw %>% rename(neighbourhood = !!name_col)
  
  pts <- st_as_sf(listings, coords = c("longitude","latitude"), crs = 4326, remove = FALSE) %>%
    st_join(neigh["neighbourhood"], left = FALSE) %>%
    dplyr::mutate(neighbourhood = neighbourhood.y) %>%
    dplyr::select(-neighbourhood.x, -neighbourhood.y)
  
  agg <- pts %>%
    st_drop_geometry() %>%
    group_by(neighbourhood) %>%
    summarise(
      n_listings = n(),
      y_home = median(log_price[is_home], na.rm=TRUE),
      y_room = median(log_price[is_room], na.rm=TRUE),
      med_accomm = median(accommodates, na.rm=TRUE),
      med_beds_home = median(bedrooms[is_home], na.rm=TRUE),
      med_baths_home= median(bathrooms[is_home], na.rm=TRUE),
      med_beds_room = median(bedrooms[is_room], na.rm=TRUE),
      med_baths_room = median(bathrooms[is_room], na.rm=TRUE),
      med_min_n = median(minimum_nights, na.rm=TRUE),
      med_host_cnt = median(host_total_listings_count, na.rm=TRUE),
      prop_super = mean(superhost, na.rm=TRUE),
      med_rev_m = median(reviews_m, na.rm=TRUE),
      med_services = median(n_amenities, na.rm=TRUE)
    ) %>%
    filter(n_listings>0, !is.na(y_home), !is.na(y_room))
  
  # Keep geometry for aggregated areas
  neigh_sf <- neigh %>%
    right_join(agg %>% select(neighbourhood), by = "neighbourhood") %>%
    st_as_sf() %>%
    st_transform(CRS_METRIC) %>%
    mutate(
      area_km2 = as.numeric(st_area(geometry))/1e6,
      centroid = st_centroid(geometry)
    )
  
  # attach density
  agg <- agg %>%
    left_join(neigh_sf %>% st_drop_geometry() %>% select(neighbourhood, area_km2), by = "neighbourhood") 
  
  # fill NAs in covariates if any
  # covar_cols <- c("med_host_cnt","prop_super","med_services","med_accomm",
  #                 "med_beds_home","med_baths_home","med_beds_room","med_baths_room",
  #                 "med_min_n","med_rev_m","dens_airbnb")
  # agg <- fill_na_with_med(agg, covar_cols)
  
  # design matrices
  X1 <- model.matrix(form_X1, data = agg)
  X2 <- model.matrix(form_X2, data = agg)

  y1 <- agg$y_home
  y2 <- agg$y_room
  
  valid <- !is.na(y1) & !is.na(y2)
  X1 <- X1[valid, , drop=FALSE]; X2 <- X2[valid, , drop=FALSE]
  y1 <- y1[valid]; y2 <- y2[valid]
  neigh_sf <- neigh_sf[valid,]


  X_blk <- as.matrix(bdiag(X1, X2))   # 2n × (k1+k2),
  Y <- c(y1, y2)

  # n <- nrow(X1)
  # perm <- as.vector(rbind(1:n, (n+1):(2*n))) # 1, n+1, 2, n+2, ...
  # X_loc <- X_blk[perm, , drop = FALSE]
  # Y_loc <- c(rbind(y1, y2))  # y_home1, y_room1, y_home2, y_room2, ...
  # X_blk <- X_loc
  # Y <- Y_loc
  
  # centroids -> distances (km)
  coords_m  <- st_coordinates(neigh_sf$centroid)
  coords_km <- coords_m/1000
  D_km <- as.matrix(dist(coords_km))
  D2 <- D_km^2
  n <- nrow(D2)
  M <- 2
  N <- n*M
  
  # median NN distance r0 (km)
  D_off <- D_km + diag(Inf, n)
  r0 <- median(apply(D_off, 1, min), na.rm = TRUE)
  
  city_area_km2 <- sum(neigh_sf$area_km2)            
  city_listings <- sum(agg$n_listings)
  city_density  <- city_listings / city_area_km2
  
  list(
    id = reg_id, n = n, N = N, M = M,
    Kbeta = ncol(X_blk), Xblk = X_blk, Y = Y,
    D2 = D2, r0 = r0,
    city_stats = list(area_km2 = city_area_km2, listings = city_listings, density = city_density)
  )
}

regions <- lapply(region_ids, process_region)


Kbeta <- unique(sapply(regions, `[[`, "Kbeta"))
M <- 2
C <- length(regions)

n_vec <- sapply(regions, `[[`, "n")
N_vec <- sapply(regions, `[[`, "N")
max_n <- max(n_vec)
max_N <- max(N_vec)

# arrays for Stan 
D2_arr <- array(0.0, dim = c(C, max_n, max_n))
X_arr <- array(0.0, dim = c(C, max_N, Kbeta))
Y_arr <- array(0.0, dim = c(C, max_N))
for (c in seq_len(C)) {
  nc <- n_vec[c]; Nc <- N_vec[c]
  D2_arr[c, 1:nc, 1:nc] <- regions[[c]]$D2
  X_arr[c, 1:Nc, 1:Kbeta] <- regions[[c]]$Xblk
  Y_arr[c, 1:Nc] <- regions[[c]]$Y
}
r0_vec <- sapply(regions, `[[`, "r0")
# mu_log_r0 <- mean(log(r0_vec)) # is this okay? centered for log l_c


stan_code_hier <- '
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
'

mod_hier <- stan_model(model_code = stan_code_hier)

stan_data_hier <- list(
  C = C, M = M, Kbeta = Kbeta,
  max_n = max_n, max_N = max_N,
  n = as.integer(n_vec), N = as.integer(N_vec),
  D2 = D2_arr, X = X_arr, Y = Y_arr,
  eps = 1e-3,                 
  # mu_log_r0 = mu_log_r0,
  log_r0 = log(r0_vec),
  eta_lkj = 2.0
)


fit.hier <- sampling(
  mod_hier, data = stan_data_hier,
  chains = 2, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.995, max_treedepth = 18)
)



# saveRDS(fit.hier, "Output/251110_uk_hier_fit.Rds")


fit.hier <- readRDS("Output/251031_uk_hier_fit.Rds")


print(fit.hier, pars = c("ell_km","r10_km","gamma_km2","rho12",
                         "beta_c","L_corr","sigma","z_ell","tau"),
      probs = c(0.025,0.5,0.975))





# compute LONO conditional mean under Mercer kernel ----
lono_region <- function(fit.hier, reg, region_ids, jitter = 1e-6) {
  ex <- rstan::extract(fit.hier, pars = c("beta_c","ell_km","R","sigma"))
  S <- dim(ex$ell_km)[1]     
  
  reg_idx <- match(reg$id, region_ids)
  n <- reg$n
  y_home <- reg$Y[1:n]
  y_room <- reg$Y[(n+1):(2*n)]
  
  mu_acc <- matrix(0.0, n, 2)
  loglik <- matrix(NA_real_, S, n)
  ypred <- array(NA_real_, dim = c(S, n, 2))
  
  for (s in seq_len(S)) {
    bvec <- ex$beta_c[s, , reg_idx]           
    mu_vec <- as.vector(reg$Xblk %*% bvec)         
    Mu <- cbind(mu_vec[1:n], mu_vec[(n+1):(2*n)])  # n x 2
    Rres <- cbind(y_home, y_room) - Mu
    
    Sigma_s <- diag(ex$sigma[s, ], 2) %*% ex$R[s, , ] %*% diag(ex$sigma[s, ], 2)
    ell <- ex$ell_km[s, reg_idx]
    K <- exp(- reg$D2 / (ell^2))
    diag(K) <- 1.0
    K <- K + diag(jitter, n)
    Kii <- diag(K)
    
    KR <- K %*% Rres
    svec <- KR - Kii * Rres
    mu_c <- Mu - sweep(svec, 1, Kii, "/")
    mu_acc <- mu_acc + mu_c
    
    for (i in 1:n) {
      loglik[s, i] <- mvtnorm::dmvnorm(c(y_home[i], y_room[i]), mean = mu_c[i, ], sigma = Sigma_s / Kii[i], log = TRUE)
    }
    Z <- mvtnorm::rmvnorm(n, sigma = Sigma_s) # marginal predictive
    Z <- sweep(Z, 1, sqrt(Kii), "/")
    ypred[s, , ] <- mu_c + Z
  }
  
  mu_hat <- mu_acc / S
  res <- cbind(y_home - mu_hat[,1], y_room - mu_hat[,2])
  rmse <- sqrt(colMeans(res^2))
  
  lpd_i <- apply(loglik, 2, function(v){ m <- max(v); m + log(mean(exp(v - m))) })
  elpd <- sum(lpd_i)
  
  qlo_home <- apply(ypred[ , , 1], 2, quantile, 0.025)
  qhi_home <- apply(ypred[ , , 1], 2, quantile, 0.975)
  qlo_room <- apply(ypred[ , , 2], 2, quantile, 0.025)
  qhi_room <- apply(ypred[ , , 2], 2, quantile, 0.975)
  cov_home <- mean(y_home >= qlo_home & y_home <= qhi_home)
  cov_room <- mean(y_room >= qlo_room & y_room <= qhi_room)
  
  list(yhat = mu_hat, res = res, rmse = rmse, elpd = elpd, cov = c(home = cov_home, room = cov_room))
}

results <- lapply(regions, function(reg) lono_region(fit.hier, reg, region_ids))


#  Moran's I 
make_listw_kernel <- function(reg, fit.hier, region_ids, tau = 0.1) {
  # compute kernel weights from posterior mean ell
  ex <- rstan::extract(fit.hier, pars = "ell_km")
  ell_hat <- mean(ex$ell_km[, match(reg$id, region_ids)])
  K <- exp(-reg$D2 / (ell_hat^2))
  diag(K) <- 0
  W <- (K > tau) * K
  row_empty <- which(rowSums(W) == 0)
  for (i in row_empty) {
    j <- which.max(K[i, ])
    W[i, j] <- K[i, j]
  }
  W <- sweep(W, 1, rowSums(W), "/")
  spdep::mat2listw(W, style = "W")
}

make_listw_knn <- function(reg, k = 6) {
  coords <- cmdscale(as.dist(reg$D2), k = 2)
  nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  spdep::nb2listw(nb, style = "W")
}


moran_tbl <- purrr::map2_dfr(regions, results, function(reg, resobj){
  lwK <- make_listw_kernel(reg, fit.hier, region_ids)
  lwNN <- make_listw_knn(reg, k = 6)
  mhK <- spdep::moran.mc(resobj$res[,1], lwK,  nsim = 999, alternative = "two.sided")
  mrK <- spdep::moran.mc(resobj$res[,2], lwK,  nsim = 999, alternative = "two.sided")
  mhN <- spdep::moran.mc(resobj$res[,1], lwNN, nsim = 999, alternative = "two.sided")
  mrN <- spdep::moran.mc(resobj$res[,2], lwNN, nsim = 999, alternative = "two.sided")
  tibble(
    region = reg$id,
    eq = rep(c("home","room"), each=2),
    W = rep(c("kernel","kNN"), times=2),
    I = c(mhK$statistic, mhN$statistic, mrK$statistic, mrN$statistic),
    p = c(mhK$p.value, mhN$p.value, mrK$p.value, mrN$p.value)
  )
})
moran_tbl











# City-level density
city_density <- sapply(regions, function(r) r$city_stats$density)
names(city_density) <- vapply(regions, `[[`, "", "id")

exH <- rstan::extract(fit.hier, pars = "ell_km")
ell_mat <- exH$ell_km 

center <- function(x) x - mean(x)
x  <- center(log(as.numeric(city_density)))
vx <- sum(x^2)

slope_draws <- apply(ell_mat, 1, function(ell) {
  y <- center(log(pmax(ell, 1e-6)))
  sum(x * y) / vx # OLS slope each draw
})
print(c(
  `2.5%` = quantile(slope_draws, 0.025),
  `50%` = quantile(slope_draws, 0.50),
  `97.5%` = quantile(slope_draws, 0.975),
  `Pr(>0)`= mean(slope_draws > 0)
))

# posterior medians and 95% CrIs of ell_c vs density 
ell_summ <- as.data.frame(ell_mat) %>%
  setNames(region_ids) %>%
  pivot_longer(everything(), names_to = "city", values_to = "ell") %>%
  group_by(city) %>%
  summarize(
    ell_med = median(ell),
    ell_lo = quantile(ell, 0.025),
    ell_hi = quantile(ell, 0.975)
  ) %>%
  mutate(density = as.numeric(city_density[city]),
         city_lab = tools::toTitleCase(city))

# posterior-median regression line on log–log scale
slope_med <- as.numeric(quantile(slope_draws, 0.5))
b0_med <- mean(log(ell_summ$ell_med)) - slope_med * mean(log(ell_summ$density))

p_scatter <- ggplot(ell_summ,
                    aes(x = log(density), y = log(ell_med), color = city)) +
  geom_pointrange(aes(ymin = log(pmax(ell_lo, 1e-6)), ymax = log(ell_hi)), show.legend = TRUE) +
  scale_color_discrete(name = "City", labels = function(x) tools::toTitleCase(x)) +
  labs(
    x = "log listing density",
    y = expression(log~"\u2113"[c])
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_line(color = scales::alpha("grey70", 0.3), linewidth = 0.2),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 13),
    axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 5)),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  )

p_scatter

ggsave("Figures/251108_hier_density.png", p_scatter, width = 10, height = 5, dpi = 300, device = "png")




