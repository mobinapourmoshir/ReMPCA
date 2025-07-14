# ## --- 0.  Packages ------------------------------------------------------------
# library(fda) # basis objects, if you need them later
# #library(MHPCA) # your hybrid‑PCA estimator
#
# ## --- 1.  Basic settings ------------------------------------------------------
# set.seed(66)
# N <- 600
# k <- 3 # hybrid rank
# p_fun <- 150 # functional variables (3 * d)
#
# ## --- 2.  Shared scores  Q  ---------------------------------------------------
# generate_coefficient_matrix <- function(N,
#                                         sd_vec = c(0, 0.9, 0.5), # <<— non‑zero SD1
#                                         seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#   N1 <- as.integer(N * 0.40)
#   N2 <- as.integer(N * 0.35)
#   N3 <- N - N1 - N2 # guarantee exact N
#
#   a <- rnorm(N1, 0, sd_vec[1])
#   b <- rnorm(N1, 0, sd_vec[2])
#   g <- rnorm(N1, 0, sd_vec[3])
#   cc <- rnorm(N2, 0, sd_vec[1])
#   d <- rnorm(N2, 0, sd_vec[2])
#   hh <- rnorm(N2, 0, sd_vec[3])
#   e <- rnorm(N3, 0, sd_vec[1])
#   f <- rnorm(N3, 0, sd_vec[2])
#   l <- rnorm(N3, 0, sd_vec[3])
#
#   Q <- rbind(
#     cbind(b, a, g),
#     cbind(hh, d, cc),
#     cbind(e, l, f)
#   )
#   qr.Q(qr(Q)) # N × 3
# }
# Q <- generate_coefficient_matrix(N)
#
# ## --- 3.  Functional loadings  V_fun  -----------------------------------------
# d <- 50
# t <- seq(0, 1, length.out = d)
#
# f11 <- c(sin(pi * t), rep(0, d), sin(2 * pi * t))
# f12 <- c(sin(2 * pi * t), rep(0, d), sin(pi * t))
# f13 <- c(rep(0, d), sin(3 * pi * t), rep(0, d))
#
# V_fun <- cbind(f11, f12, f13) # 150 × 3
# V_fun <- qr.Q(qr(V_fun))
#
# ## --- 3.  Vector loadings  V_vec  (100 × 3, 50% sparsity) --------------------
#
#
# p_vec <- 10 # new vector‑block size
# set.seed(123)
# v_tilde_1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
# v_tilde_2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
# v_tilde_3 <-  c(1, 1, -1, -1, 1, -1, 0, 0, 0, 0)
#
# # v_tilde <- matrix(rnorm(300,5,sd = 2),ncol=3)
# # n_rows <- nrow(v_tilde)
# # zero_fraction <- 0.4
# #
# # for (j in 1:ncol(v_tilde)) {
# #   zero_indices <- sample(1:n_rows, size = floor(n_rows * zero_fraction), replace = FALSE)
# #   v_tilde[zero_indices, j] <- 0
# # }
# v_tilde <- cbind(v_tilde_1,v_tilde_2,v_tilde_3)
# V_vec <- qr.Q(qr(v_tilde))
#
# ## --- 5.  Data generation -----------------------------------------------------
# ## Noise SD can be tuned; keep 1 for both blocks
# eps_fun <- matrix(rnorm(N * p_fun, sd = .1), N, p_fun)
# eps_vec <- matrix(rnorm(N * p_vec, sd = .1), N,p_vec)## p_vec)
#
# eig_val <- c(200,100,50)
# # variations
# 100*(eig_val/sum(eig_val))
# D <- diag(sqrt(eig_val))
#
# Y <- Q %*% D %*% t(V_fun) + eps_fun # N × 150
# X <- Q %*% D %*% t(V_vec) + eps_vec # N × 10
#
# ## --- 6.  Quick sanity checks -------------------------------------------------
# #(i) Per‑variable signal variance should be similar
# #sig_var_fun <- var(as.vector((Q %*% D %*% t(V_fun))))
# #sig_var_vec <- var(as.vector((Q %*% D %*% t(V_vec))))
# #round(c(fun = sig_var_fun, vec = sig_var_vec), 3)
# #(ii) Optional: standardise X and Z columns before fitting HPCA
#
#
#
# ## --- 7.  Object definition ---------------------------------------------------
#
# library(ReMPCA)
# Yfd <- fdClass(Y, Smoothing_parameter = NULL,
#                Sparsity_parameter = c(0, 1, 2, 3, 4, 8, 16, 32, 64, 75, 85, 95, 100, 110, 128, 149))
# Xrd <- rdClass(X, Sparsity_parameter = NULL)
#
# hdobj <- hdClass(list(Xrd,Yfd), Smoothing_parameter = 0,
#                  Sparsity_parameter = NULL)
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
# plot_pc_functions(SimulTest)
# plot_pc_scores(SimulTest)
#
# SimulTest$VarianceExplained
#
# par(mfrow = c(3,2))
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
#
#
#
#
# X_temp =  X_temp
# n_var = n_var
# ncol = ncol
# n = n
# GridPoints_u = GridPoints_u
# GridPoints_v = GridPoints_v
# smooth_tuning_v = smooth_tuning_col
# smooth_tuning_u = smooth_tuning_row
# sparse_tuning_u = sparsity_row_list
# sparse_tuning_v = sparsity_col_list
# sparse_tuning_type = sparse_tuning_type
# nfolds_u = nfolds_u
# nfolds_v = nfolds_v
# parallel = parallel
# S_alpha_list_v = S_alpha_list_v
# S_alpha_list_u = S_alpha_list_u
# Omegas_u = Omegas_u
# Omegas_v = Omegas_v
# tuning_iter = tuning_iter
# tuning_order = tuning_order
# thresh = thresh
# maxit = maxit
# cv.pick = cv.pick
# smoothness_type = smoothness_type
#
#
#
#
#
#
# X = X_temp
# n_var = n_var
# ncol = ncol
# thresh = thresh
# maxit = maxit
# conditional = TRUE
# S_alphas_v = S_alpha_list_v
# S_alphas_u = opt_u$opt_s.alpha_u
# alphas_v = smooth_tuning_v
# alpha_u = alpha_u
# Omega_u = opt_u$Omega_u
# Omegas_v = Omegas_v
# sparse_tuning_result_u = gamma_u
# sparse_tuning_result_v = gamma_v
# sparse_tuning_type = sparse_tuning_type
#
#
# ############## Inner Product
# v1 <- c(1.2,3.4,-2.1)
# v2 <- c(0.5,-1,4)
#
# t(v1) %*% v2
#
# sum(v1 *v2)
# # dx = 1 is used because if no spacing is specified, and it’s the default assumption for equal, unit spacing.
#
# # Having grid points:
# grid <- seq(0, 1, length.out = 3)
# dx <- diff(grid)[1]  # = 0.5
#
# sum(v1 * v2) * dx
#
#
# # Functional and approximation
# f <- function(x) x
# g <- function(x) sin(x)
# product_fg <- function(x) f(x) * g(x)
# integrate(product_fg, lower = 0, upper = pi)$value
#
#
# # Approximate
# n <- 1000                         # number of intervals
# a <- 0                           # start of interval
# b <- pi                          # end of interval
# x_vals <- seq(a, b, length.out = n)
# h <- (b - a) / (n - 1)           # step size
#
# # Riemann sum approximation
# approx_inner_product <- sum(f(x_vals) * g(x_vals)) * h
# cat("Approximate inner product: ", approx_inner_product, "\n")
#
#
#
#
#
# # Inner product test
# par(mfrow = c(3,2))
# x <- seq(0,1,len = 100)*2*pi
# v <- sin(x)
# norm_vec(v)
# dx <- diff(x)[1]
# sqrt(t(v)%*%v * dx)
# W1 = diag(100)*dx
# sqrt(t(v)%*% W1%*%v)
#
# plot(v, main = "v1")
# # v <- v/norm_vec(v)
# # plot(v)
# # norm_vec(v)
#
# x2 <- seq(0,1,len = 200)*2*pi
# v2 <- sin(x2)
# dx2 <- diff(x2)[1]
# sqrt(t(v2)%*%v2 * dx2)
# W2 = diag(200)*dx2
# sqrt(t(v2)%*% W2%*%v2)
#
# norm_vec(v2)
# plot(v2, main = "v2")
# #v2 <- v2/norm_vec(v2)
# #plot(v2)
# #norm_vec(v2)
#
# set.seed(11)
# u <- rnorm(50)
# X <- u%*%t(v)
# X2 <- u%*%t(v2)
#
#
# matplot(svd(X)$v[,1], type = 'l', main ="SVD")
# matplot(svd(X2)$v[,1], type = 'l', main ="SVD")
#
# xfd <- fdClass(X, Smoothing_parameter = 0)
# xobj <-hdClass(list(xfd))
# xpca <- ReMPCA(xobj)
#
# matplot(xpca$PCFunctions[[1]][[1]], type = 'l', main = "MPCA v1")
#
#
# xfd2 <- fdClass(X2, Smoothing_parameter = 0)
# xobj2 <-hdClass(list(xfd2))
# xpca2 <- ReMPCA(xobj2)
#
# matplot(xpca2$PCFunctions[[1]][[1]], type = 'l', main = "MPCA v2")
#
#
# #multivariate
# hdxx2obj <-hdClass(list(xfd, xfd2), Smoothing_parameter = 0)
# xx2pca <- ReMPCA(hdxx2obj)
#
#
