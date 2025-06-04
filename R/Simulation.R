## --- 0.  Packages ------------------------------------------------------------
library(fda) # basis objects, if you need them later
#library(MHPCA) # your hybrid‑PCA estimator

## --- 1.  Basic settings ------------------------------------------------------
set.seed(66)
N <- 600
k <- 3 # hybrid rank
p_fun <- 150 # functional variables (3 * d)

## --- 2.  Shared scores  Q  ---------------------------------------------------
generate_coefficient_matrix <- function(N,
                                        sd_vec = c(0, 0.9, 0.5), # <<— non‑zero SD1
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N1 <- as.integer(N * 0.40)
  N2 <- as.integer(N * 0.35)
  N3 <- N - N1 - N2 # guarantee exact N

  a <- rnorm(N1, 0, sd_vec[1])
  b <- rnorm(N1, 0, sd_vec[2])
  g <- rnorm(N1, 0, sd_vec[3])
  cc <- rnorm(N2, 0, sd_vec[1])
  d <- rnorm(N2, 0, sd_vec[2])
  hh <- rnorm(N2, 0, sd_vec[3])
  e <- rnorm(N3, 0, sd_vec[1])
  f <- rnorm(N3, 0, sd_vec[2])
  l <- rnorm(N3, 0, sd_vec[3])

  Q <- rbind(
    cbind(b, a, g),
    cbind(hh, d, cc),
    cbind(e, l, f)
  )
  qr.Q(qr(Q)) # N × 3
}
Q <- generate_coefficient_matrix(N)

## --- 3.  Functional loadings  V_fun  -----------------------------------------
d <- 50
t <- seq(0, 1, length.out = d)

f11 <- c(sin(pi * t), rep(0, d), sin(2 * pi * t))
f12 <- c(sin(2 * pi * t), rep(0, d), sin(pi * t))
f13 <- c(rep(0, d), sin(3 * pi * t), rep(0, d))

V_fun <- cbind(f11, f12, f13) # 150 × 3
V_fun <- qr.Q(qr(V_fun))

## --- 3.  Vector loadings  V_vec  (100 × 3, 50% sparsity) --------------------


p_vec <- 10 # new vector‑block size
set.seed(123)
v_tilde_1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
v_tilde_2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
v_tilde_3 <-  c(1, 1, -1, -1, 1, -1, 0, 0, 0, 0)

# v_tilde <- matrix(rnorm(300,5,sd = 2),ncol=3)
# n_rows <- nrow(v_tilde)
# zero_fraction <- 0.4
#
# for (j in 1:ncol(v_tilde)) {
#   zero_indices <- sample(1:n_rows, size = floor(n_rows * zero_fraction), replace = FALSE)
#   v_tilde[zero_indices, j] <- 0
# }
v_tilde <- cbind(v_tilde_1,v_tilde_2,v_tilde_3)
V_vec <- qr.Q(qr(v_tilde))

## --- 5.  Data generation -----------------------------------------------------
## Noise SD can be tuned; keep 1 for both blocks
eps_fun <- matrix(rnorm(N * p_fun, sd = .1), N, p_fun)
eps_vec <- matrix(rnorm(N * p_vec, sd = .1), N,p_vec)## p_vec)

eig_val <- c(200,100,50)
# variations
100*(eig_val/sum(eig_val))
D <- diag(sqrt(eig_val))

Y <- Q %*% D %*% t(V_fun) + eps_fun # N × 150
X <- Q %*% D %*% t(V_vec) + eps_vec # N × 10

## --- 6.  Quick sanity checks -------------------------------------------------
# (i) Per‑variable signal variance should be similar
#sig_var_fun <- var(as.vector((Q %*% D %*% t(V_fun))))
#sig_var_vec <- var(as.vector((Q %*% D %*% t(V_vec))))
#round(c(fun = sig_var_fun, vec = sig_var_vec), 3)
# (ii) Optional: standardise X and Z columns before fitting HPCA



## --- 7.  Object definition ---------------------------------------------------

# library(ReMPCA)
# Yfd <- fdClass(Y, Smoothing_parameter = 0, Sparsity_parameter = 0)
# Xrd <- rdClass(X, Sparsity_parameter = 0)
#
# hdobj <- hdClass(list(Xrd,Yfd), Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
#
# SimulTest <- ReMPCA(hd = hdobj,
#                      centerhds = TRUE,
#                      num_pcs = 3,
#                      nfolds_u = 5,
#                      nfolds_v = NULL,
#                      thresh = 1e-10,
#                      maxit = 100,
#                      tuning_iter = 1,
#                      parallel = FALSE,
#                      weights = 0,
#                      smoothness_type = "Second_order",
#                      sparse_tuning_type = "soft",
#                      tuning_order = "Sparsity",
#                      cv.pick = "1se",
#                      sparse_tuning_u = NULL,
#                      sparse_tuning_v = NULL,
#                      smooth_tuning_u = NULL,
#                      smooth_tuning_v = NULL)
#
#
# matplot(SimulTest$PCFunctions[[1]][[1]], type = 'l')
# matplot(SimulTest$PCFunctions[[1]][[2]], type = 'l')
# matplot(SimulTest$PCFunctions[[2]][[1]], type = 'l')
# matplot(SimulTest$PCFunctions[[2]][[2]], type = 'l')
# matplot(SimulTest$PCFunctions[[3]][[1]], type = 'l')
# matplot(SimulTest$PCFunctions[[3]][[2]], type = 'l')
#
# hd = hdobj
# centerhds = TRUE
# num_pcs = 3
# nfolds_u = 5
# nfolds_v = NULL
# thresh = 1e-10
# maxit = 100
# tuning_iter = 1
# parallel = FALSE
# weights = 0
# smoothness_type = "Second_order"
# sparse_tuning_type = "soft"
# tuning_order = "Sparsity"
# cv.pick = "1se"
# sparse_tuning_u = NULL
# sparse_tuning_v = NULL
# smooth_tuning_u = NULL
# smooth_tuning_v = NULL
