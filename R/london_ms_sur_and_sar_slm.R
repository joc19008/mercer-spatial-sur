# Mercer SAR–SUR & SUR–SLM on London boroughs

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


sf::sf_use_s2(FALSE)

listings_path <- "Data/listings_london.csv"
neigh_path <- "Data/neighbourhoods_london.geojson"

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

# Read neighborhoods
neigh_raw <- st_read(neigh_path, quiet = TRUE) %>% st_transform(4326)
nm <- tolower(names(neigh_raw))
name_col <- if ("neighbourhood" %in% nm) names(neigh_raw)[which(nm=="neighbourhood")[1]] else
  if ("name" %in% nm) names(neigh_raw)[which(nm=="name")[1]] else names(neigh_raw)[1]
neigh <- neigh_raw %>% rename(neighbourhood = !!name_col)

pts <- st_as_sf(listings, coords = c("longitude","latitude"), crs = 4326, remove = FALSE) %>%
  st_join(neigh["neighbourhood"], left = FALSE) %>%
  dplyr::mutate(neighbourhood = neighbourhood.y) %>%
  dplyr::select(-neighbourhood.x, -neighbourhood.y)


# Aggregate to borough
agg <- pts %>%
  st_drop_geometry() %>%
  group_by(neighbourhood) %>%
  summarise(
    n_listings = n(),
    # dependent vars: median log-price
    y_home = median(log_price[is_home], na.rm=TRUE),
    y_room = median(log_price[is_room], na.rm=TRUE),
    # covariates (medians / proportions)
    med_accomm = median(accommodates, na.rm=TRUE),
    med_beds_home = median(bedrooms[is_home], na.rm=TRUE),
    med_baths_home = median(bathrooms[is_home], na.rm=TRUE),
    med_beds_room = median(bedrooms[is_room], na.rm=TRUE),
    med_baths_room = median(bathrooms[is_room], na.rm=TRUE),
    med_min_n  = median(minimum_nights, na.rm=TRUE),
    med_host_cnt = median(host_total_listings_count, na.rm=TRUE),
    prop_super = mean(superhost, na.rm=TRUE),
    med_rev_m  = median(reviews_m, na.rm=TRUE),
    med_services = median(n_amenities, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  filter(n_listings>0, !is.na(y_home), !is.na(y_room))

# Keep geometry only for those boroughs present after aggregation
neigh_sf <- neigh %>%
  right_join(agg %>% select(neighbourhood), by = "neighbourhood") %>%
  st_as_sf()

# Project to metric CRS and compute areas, centroids
CRS_METRIC <- 27700  # OSGB 1936 / British National Grid (meters)
neigh_sf <- st_transform(neigh_sf, CRS_METRIC) %>%
  mutate(
    area_km2 = as.numeric(st_area(geometry))/1e6,
    centroid = st_centroid(geometry)
  )

coords_m <- st_coordinates(neigh_sf$centroid) # meters
coords_km <- coords_m / 1000  # km
n_boros <- nrow(neigh_sf)

# Distances to POIs (km): Trafalgar Sq, Wembley, Heathrow T5
# trafalgar <- st_sfc(st_point(c(-0.1278, 51.5074)), crs = 4326) %>% st_transform(CRS_METRIC)
# wembley   <- st_sfc(st_point(c(-0.2796, 51.5560)), crs = 4326) %>% st_transform(CRS_METRIC)
# heathrow  <- st_sfc(st_point(c(-0.4890, 51.4700)), crs = 4326) %>% st_transform(CRS_METRIC)
# 
# dist_km <- function(P) as.numeric(st_distance(neigh_sf$centroid, P))/1000
# 
# agg2 <- agg %>%
#   mutate(
#     dist_center  = dist_km(trafalgar),
#     dist_stadium = dist_km(wembley),
#     dist_airport = dist_km(heathrow)
#   ) %>%
#   mutate(dens_airbnb = n_listings / neigh_sf$area_km2)

agg <- agg %>%
  left_join(neigh_sf, by = "neighbourhood") %>%
  mutate(dens_airbnb = n_listings / area_km2)

# Design matrices & responses
X1 <- model.matrix(~ 1 + med_host_cnt + prop_super + med_services + med_accomm +
                     med_beds_home + med_baths_home + med_min_n + med_rev_m, # + dens_airbnb, 
                   # + dist_airport + dist_stadium + dist_center + dens_airbnb,
                   data = agg)

X2 <- model.matrix(~ 1 + med_host_cnt + prop_super + med_services + med_accomm +
                     med_baths_room + med_min_n + med_rev_m, # + dens_airbnb,
                     # + dist_airport + dist_stadium + dist_center + dens_airbnb,
                   data = agg)

y1 <- agg$y_home
y2 <- agg$y_room

valid <- !is.na(y1) & !is.na(y2)
X1 <- X1[valid,]; X2 <- X2[valid,]
y1 <- y1[valid];  y2 <- y2[valid]
neigh_sf <- neigh_sf[valid,]

# Distances (km) and D^2 (km^2) on centroids
coords_km_valid <- coords_km[valid, , drop=FALSE]
D_km <- as.matrix(dist(coords_km_valid))
D2 <- D_km^2
n <- nrow(D2)
M <- 2
N <- n*M

# Median nearest-neighbour distance r0 (km) for log-ell prior center
D_off <- D_km + diag(Inf, n)
r0 <- median(apply(D_off, 1, min), na.rm = TRUE)

X_blk <- as.matrix(bdiag(X1, X2)) # 2n x (k1+k2)
perm <- as.vector(rbind(1:n, (n+1):(2*n)))    
X_loc <- X_blk[perm,]
Y_loc <- c(rbind(y1, y2))

Kbeta <- ncol(X_blk)


# Fit Mercer SAR SUR
mod_mercer <- rstan::stan_model(file = "stan/mercer_sar_sur_london.stan")


stan_data_mercer <- list(
  n = n, M = M, N = N, Kbeta = Kbeta,
  D2 = D2, X = X_loc, Y = Y_loc,
  eps = 1e-6,
  mu_ell = log(r0), sd_ell = 0.5,            
  u_beta = rep(0, Kbeta), V_beta = diag(1, Kbeta),
  m_tau = rep(log(0.5), M), s_tau = rep(0.35, M),
  eta_lkj = 2.0
)


fit.mercer <- sampling(
  mod_mercer, data = stan_data_mercer,
  chains = 2, iter = 2000, warmup = 1000, 
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)


# saveRDS(fit.mercer, "Output/251024_london_mercer_fit.Rds")

fit.mercer <- readRDS("Output/251024_london_mercer_fit.Rds")

print(fit.mercer, pars = c("beta","Sigma","ell_km","r10","rho12","sigma","tau"),
      probs = c(0.025,0.5,0.975))






# SUR–SLM comparator in Stan

## kNN W (row-standardized) for SUR–SLM 
make_W <- function(coords, k) {
  n <- nrow(coords)
  idx <- FNN::get.knn(coords, k = k + 1)$nn.index
  idx <- t(vapply(seq_len(n), function(i) idx[i, idx[i,] != i][1:k], integer(k)))
  ii <- rep(seq_len(n), each = k)
  jj <- as.vector(t(idx))
  W <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))
  diag(W) <- 0
  rs <- rowSums(W)
  stopifnot(all(rs > 0))
  W <- Diagonal(x = 1/rs) %*% W
  as.matrix(W)
}
W <- make_W(coords_km_valid, 6)
# W <- make_W(coords_km_valid, 4)

# Bounds for rho (directed W allowed)
# ev <- eigen(W6, only.values = TRUE)$values
# lam_min <- min(Re(ev))
# rho_lower <- 1/lam_min  # negative
# rho_upper <- 1.0

ev <- eigen(W, only.values = TRUE)$values
rho_upper <- 0.999 / max(Mod(ev))
rho_lower <- -rho_upper




# Fit SUR–SLM
stan_sur <- rstan::stan_model(file = "stan/sur_slm.stan")


stan_data_sar <- list(
  n  = nrow(X1),
  K1 = ncol(X1), K2 = ncol(X2),
  X1 = X1, X2 = X2, y1 = y1, y2 = y2,
  W  = W,
  rho_lower = rho_lower, rho_upper = rho_upper,
  Omega = diag(2), nu = 4
)

fit.sar <- sampling(stan_sur, data = stan_data_sar,
                       iter=2000, warmup=1000, chains=2,
                       control = list(adapt_delta = 0.9, max_treedepth = 13))

# saveRDS(fit.sar, "Output/251009_london_sar_W6_fit.Rds")

fit.sar <- readRDS("Output/251009_london_sar_W6_fit.Rds")

print(fit.sar)


# saveRDS(fit.sar, "Output/250902_london_sar_W6_fit.Rds")






# Pointwise log-likelihoods for ELPD LOO
## Mercer pointwise conditional log-lik
exM <- rstan::extract(fit.mercer, pars = c("beta","ell_km","R","sigma"))  # 'sigma' here equals tau
S <- length(exM$ell_km)
K1c <- ncol(X1); K2c <- ncol(X2)
split_beta <- function(b) list(b1=b[1:K1c], b2=b[(K1c+1):(K1c+K2c)])

loglik_mercer <- matrix(NA, S, n)
for (s in seq_len(S)) {
  b <- split_beta(exM$beta[s,])
  Mu <- cbind(as.vector(X1 %*% b$b1), as.vector(X2 %*% b$b2))
  Rres <- cbind(y1,y2) - Mu
  
  Sigma_s <- (diag(exM$sigma[s,], 2) %*% exM$R[s,,] %*% diag(exM$sigma[s,], 2))
  
  ell <- exM$ell_km[s]
  Kmat <- exp( - D2 / (ell*ell) ) + diag(1e-6, n)
  Kii <- diag(Kmat)
  
  KR <- Kmat %*% Rres
  svec <- KR - Kii*Rres        
  mu_c <- Mu - svec/Kii     
  
  for (i in 1:n) {
    loglik_mercer[s,i] <- mvtnorm::dmvnorm(
      x = c(y1[i], y2[i]), mean = mu_c[i,], sigma = Sigma_s/Kii[i], log = TRUE
    )
  }
}
loo_mercer <- loo::loo(loglik_mercer)
print(loo_mercer)



## SUR–SLM pointwise conditional log-lik
exS <- rstan::extract(fit.sar, pars = c("beta1","beta2","rho1","rho2","Sinv"))
S_s <- length(exS$rho1)
I_n <- Matrix::Diagonal(n)

loglik_sar <- matrix(NA_real_, S_s, n)
for (s in seq_len(S_s)) {
  A1 <- as.matrix(I_n - exS$rho1[s] * W)
  A2 <- as.matrix(I_n - exS$rho2[s] * W)
  mu1 <- as.vector(solve(A1, X1 %*% exS$beta1[s,]))
  mu2 <- as.vector(solve(A2, X2 %*% exS$beta2[s,]))
  ytilde1 <- y1 - mu1; ytilde2 <- y2 - mu2
  
  Sinv_s <- matrix(exS$Sinv[s,,], 2, 2); S_s2 <- solve(Sinv_s)
  s11 <- Sinv_s[1,1]; s12 <- Sinv_s[1,2]; s22 <- Sinv_s[2,2]
  
  B11 <- crossprod(A1)            
  B12 <- t(A1) %*% A2          
  B22 <- crossprod(A2)     
  
  d11 <- s11*diag(B11); d12 <- s12*diag(B12); d22 <- s22*diag(B22)
  b1_full <- s11*(B11 %*% ytilde1) + s12*(B12 %*% ytilde2)
  b2_full <- s12*(t(B12) %*% ytilde1) + s22*(B22 %*% ytilde2)
  
  for (i in 1:n) {
    b1_exc <- b1_full[i] - s11*B11[i,i]*ytilde1[i] - s12*B12[i,i]*ytilde2[i]
    b2_exc <- b2_full[i] - s12*B12[i,i]*ytilde1[i] - s22*B22[i,i]*ytilde2[i]
    Qii <- matrix(c(d11[i], d12[i], d12[i], d22[i]), 2, 2)
    mi  <- c(mu1[i], mu2[i]) - solve(Qii, c(b1_exc, b2_exc))
    loglik_sar[s,i] <- mvtnorm::dmvnorm(c(y1[i], y2[i]), mean = mi, sigma = solve(Qii), log = TRUE)
  }
}
loo_sar <- loo::loo(loglik_sar)
print(loo_sar)


loo::loo_compare(loo_mercer, loo_sar)



# distance–decay 
ex  <- rstan::extract(fit.mercer, pars = "ell_km")
ell_draws <- as.numeric(ex$ell_km)

r_grid <- seq(0, max(D_km), length.out = 300)
S_plot <- min(1000, length(ell_draws))
ids <- sample(seq_along(ell_draws), S_plot)

Kgrid <- sapply(ids, function(s) exp( - (r_grid^2) / (ell_draws[s]^2) ))
df_k  <- data.frame(
  r = r_grid,
  mean = rowMeans(Kgrid),
  lo = apply(Kgrid, 1, quantile, 0.025),
  hi = apply(Kgrid, 1, quantile, 0.975)
)

decay_plot <- ggplot(df_k, aes(r, mean)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, fill = "#6BAED6") +
  geom_line(color = "#08519C", linewidth = 0.9) +
  xlim(0, 30) + 
  labs(x = "Distance r (km)", y = "Kernel weight K(r)") +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
decay_plot

# ggsave("Figures/251030_london_decay.png", decay_plot, width = 8, height = 6, dpi = 300, device = "png")


# Effective neighbors (> tau)
k_thr <- 0.1 
S_eff <- min(500, length(ell_draws))
ids_e <- sample(seq_along(ell_draws), S_eff)

Eff <- sapply(ids_e, function(s) {
  K <- exp( - D2 / (ell_draws[s]^2) )
  diag(K) <- 0
  rowSums(K > k_thr)
})
neigh_sf$eff_neigh_mean <- as.numeric(rowMeans(Eff))

cor_plot <- ggplot(neigh_sf) +
  geom_sf(aes(fill = eff_neigh_mean), color = NA) +
  viridis::scale_fill_viridis(
    name = "Effective neighbors",
    breaks = c(0.0, 1.0, 2.0)
  ) +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
  )
cor_plot

# ggsave("Figures/251030_london_map.png", cor_plot, width = 8, height = 6, dpi = 300, device = "png")



bothplot <- gridExtra::grid.arrange(decay_plot, cor_plot, ncol=2)


bothplot <- plot_grid(
  decay_plot, cor_plot,
  ncol = 2# ,
  # rel_widths = c(0.9, 1)#,
  #labels = c("", "(C)")
)


ggsave("Figures/251030_london.png", bothplot, width = 9, height = 4, dpi = 300, device = "png")







# EDA
library(cowplot)
## Log-price densities
lab_rt <- c("Entire home/apt"="Homes","Private room"="Rooms")
cols_rt <- c("Entire home/apt"="darkorange","Private room"="steelblue")

p_dens <- listings %>%
  ggplot(aes(x = log_price, color = room_type, fill = room_type)) +
  geom_density(alpha = 0.15, linewidth = 1) +
  scale_color_manual(values = cols_rt, labels = lab_rt, name = "Type") +
  scale_fill_manual(values = cols_rt, labels = lab_rt, name = "Type") +
  labs(x = "log(price)", y = "Density") +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ) + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0, "pt"),
    legend.box.margin = margin(t = 2, b = -5),
    plot.margin = margin(2, 2, 2, 2)
  )

p_dens



## NN distance histogram with r0
D_off <- D_km + diag(Inf, nrow(D_km))
nn1 <- apply(D_off, 1, min)
r0 <- median(nn1, na.rm = TRUE)

p_nnhist <- data.frame(nn_km = nn1) %>%
  ggplot(aes(nn_km)) +
  geom_histogram(bins = 12, fill = "steelblue1", color = "white") +
  geom_vline(xintercept = r0, linetype = 2, color = "black") +
  # annotate("text", x = r0, y = Inf, vjust = 1.6, label = "median NN = r0", size = 3) +
  labs(x = "Nearest-neighbour distance (km)", y = "Count") +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )
p_nnhist



## borough maps (homes & rooms medians)
lims <- range(c(agg$y_home, agg$y_room), na.rm = TRUE)

map_home <- neigh_sf %>%
  left_join(agg %>% select(neighbourhood, y_home), by = "neighbourhood") %>%
  ggplot() +
  geom_sf(aes(fill = y_home), color = NA) +
  scale_fill_viridis(option = "C", limits = lims, name = "Median log-price") +
  # labs(title = "Homes: borough median log-price") +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ) + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.box.spacing = unit(0, "pt"),
    legend.box.margin = margin(t = -5, b = -5),
    plot.margin = margin(2, 2, 2, 2)
  )

map_room <- neigh_sf %>%
  left_join(agg %>% select(neighbourhood, y_room), by = "neighbourhood") %>%
  ggplot() +
  geom_sf(aes(fill = y_room), color = NA) +
  scale_fill_viridis(option = "C", limits = lims, name = "Median log-price") +
  # labs(title = "Rooms: borough median log-price") +
  theme_bw(base_size = 18) +
  theme(
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  ) + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.4, "cm"),
    legend.box.spacing = unit(0, "pt"),
    legend.box.margin = margin(t = -5, b = -5),
    plot.margin = margin(2, 2, 2, 2)
  )


map_home
map_room


legend <- get_legend(map_home + theme(legend.box.margin = margin(0, 0, 0, 0)))

map_home_noleg <- map_home + theme(legend.position = "none")
map_room_noleg <- map_room + theme(legend.position = "none")

maps_combined <- plot_grid(map_home_noleg, map_room_noleg, ncol = 2, align = "hv")
maps_combined <- plot_grid(maps_combined, legend, ncol = 1, rel_heights = c(1, 0.20))

# top_row <- plot_grid(p_dens, p_nnhist, ncol = 2, labels = c("(A)", "(B)"))
top_row <- plot_grid(p_dens, p_nnhist, ncol = 2)


london_eda <- plot_grid(
  top_row, maps_combined,
  ncol = 1,
  rel_heights = c(0.9, 1.2)#,
  #labels = c("", "(C)")
)
london_eda


ggsave("Figures/251030_london_eda.png", london_eda, width = 15, height = 10, dpi = 300)



# Residual map

# LOO posterior-predictive mean (same conditional math we used for elpd-LOO)
loo_mu_mercer_london <- function(fit_mercer, X1, X2, y1, y2, D2, jitter = 1e-6) {
  ex <- rstan::extract(fit_mercer, pars = c("beta","ell_km"))
  S_all <- length(ex$ell_km)  # number of posterior draws
  K1c <- ncol(X1); K2c <- ncol(X2); n <- nrow(X1)
  mu_acc <- matrix(0.0, n, 2)
  
  split_beta <- function(b) list(b1 = b[1:K1c], b2 = b[(K1c+1):(K1c+K2c)])
  
  for (s in seq_len(S_all)) {
    b <- split_beta(ex$beta[s, ])
    Mu <- cbind(as.vector(X1 %*% b$b1), as.vector(X2 %*% b$b2))  # n x 2
    Rres <- cbind(y1, y2) - Mu
    
    ell <- ex$ell_km[s]
    Kmat <- exp(- D2 / (ell^2))
    diag(Kmat) <- 1.0
    Kmat <- Kmat + diag(jitter, n) # numerical stability
    Kii <- diag(Kmat)
    
    # conditional mean y_i | y_-i (for both equations jointly)
    KR <- Kmat %*% Rres
    svec <- KR - Kii * Rres      
    mu_cv <- Mu - sweep(svec, 1, Kii, "/")
    mu_acc <- mu_acc + mu_cv
  }
  mu_acc / S_all
}

mu_cv_lon <- loo_mu_mercer_london(fit.mercer, X1, X2, y1, y2, D2)
res_home <- y1 - mu_cv_lon[,1]
res_room <- y2 - mu_cv_lon[,2]

lon_res_sf <- neigh_sf %>%
  st_as_sf() %>%
  mutate(res_home = res_home, res_room = res_room)


Lmax <- max(abs(c(res_home, res_room)), na.rm = TRUE)

residual_map <- function(df, column, title) {
  ggplot(df) +
    geom_sf(aes(fill = .data[[column]]), color = NA) +
    scale_fill_gradient2(
      name = "Residual",
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = c(-Lmax, Lmax)
    ) +
    theme_bw(base_size = 18) +
    theme(
      axis.title = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    ) + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.key.width = unit(1.0, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(t = -5, b = -5),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p_home <- residual_map(lon_res_sf, "res_home", "Homes")
p_room <- residual_map(lon_res_sf, "res_room", "Rooms")


leg <- cowplot::get_legend(p_home +
                             theme(legend.position = "bottom",
                                 legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0)))

p_home_nl <- p_home + theme(legend.position = "none")
p_room_nl <- p_room + theme(legend.position = "none")

panel <- cowplot::plot_grid(p_home_nl, p_room_nl, ncol = 2, align = "hv")
p_resid <- cowplot::plot_grid(panel, leg, ncol = 1, rel_heights = c(1, 0.20))
p_resid

ggsave("Figures/251110_london_resid.png", p_resid, width = 10, height = 5, dpi = 300)
