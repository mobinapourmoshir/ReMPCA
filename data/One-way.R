# library(ReMPCA)
# ############### Simulation for two functional variables - one-way ###############
# OneWaySimulation <- function(N, sigma,random_seed){
#   norm_vec <- function(x) sqrt(sum(x^2))
#   rescale <- function(x, to = c(0, 1)) {
#     rng <- range(x)
#     (x - rng[1]) / diff(rng) * (to[2] - to[1]) + to[1]
#   }
#   mse_L2 <- function(true, est) {
#     # true <- true / sqrt(sum(true^2))
#     # est <- est / sqrt(sum(est^2))
#     # err1 <- sum((true - est)^2)
#     # err2 <- sum((true + est)^2)
#     # min(err1, err2) / length(true)  # approximate L2 by Riemann sum
#     mean((true - est)^2)
#   }
#   set.seed(random_seed)
#
#   sigma1 <- 20
#   sigma2 <- 10
#   sigma <- sigma
#   n <- m <- N
#   t <- seq(-1,1,length.out = m)
#
#   # generate scores
#   u11 <- rnorm(n, 0, sigma1)
#   u12 <- rnorm(n, 0, sigma2)
#
#   # generate FPCs for first variable
#   v11 <- t + sin(pi*t)
#   v11 <- v11 / norm_vec(v11)
#
#   v12 <- cos(3*pi*t)
#   v12 <- v12 / norm_vec(v12)
#
#   # generate FPCs for second variable
#   block1 <- t <= -1/3
#   block2 <- (t > -1/3) & (t < 1/3)
#   block3 <- t >= 1/3
#
#   # Define three sparse functions
#   f1 <- numeric(m)
#   f1[block1] <- sin(pi * rescale(t[block1], to = c(0, 1)))
#   f1[block3] <- sin(2 * pi * rescale(t[block3], to = c(0, 1)))
#
#   f2 <- numeric(m)
#   f2[block1] <- sin(2 * pi * rescale(t[block1], to = c(0, 1)))
#   f2[block3] <- sin(pi * rescale(t[block3], to = c(0, 1)))
#
#   f3 <- numeric(m)
#   f3[block2] <- sin(3 * pi * rescale(t[block2], to = c(0, 1)))
#
#   # Normalize
#   v21 <- f3 / sqrt(sum(f3^2))
#   v22 <- f2 / sqrt(sum(f2^2))
#
#   # generate noise
#   eps1 <- matrix(rnorm(n*m, 0, sigma), n, m)
#   eps2 <- matrix(rnorm(n*m, 0, sigma), n, m)
#
#   # data matrices
#   X1 <- u11 %*% t(v11) + u12 %*% t(v12) + eps1
#   X2 <- u11 %*% t(v21) + u12 %*% t(v22) + eps2
#
#   # combine side by side
#   X_all <- cbind(X1, X2)
#
#   ###### SVD ######
#
#   ###### Smooth Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0)
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0)
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#   print("Multivariate SVD")
#   SimulTest_svd_multi <- ReMPCA(hd = Xobj,
#                                    centerhds = TRUE,
#                                    num_pcs = 2,
#                                    nfolds_u = 5,
#                                    nfolds_v = NULL,
#                                    thresh = 1e-10,
#                                    maxit = 100,
#                                    tuning_iter = 1,
#                                    parallel = FALSE,
#                                    weights = 0,
#                                    smoothness_type = "Second_order",
#                                    sparse_tuning_type = "SCAD",
#                                    tuning_order = "Sparsity",
#                                    cv.pick = "min",
#                                    sparse_tuning_u = NULL,
#                                    sparse_tuning_v = NULL,
#                                    smooth_tuning_u = NULL,
#                                    smooth_tuning_v = NULL)
#
#
#   # Results
#   # 1st PC, 1st variable
#   v11_est_svd_multi <- SimulTest_svd_multi$PCFunctions[[1]][[1]]
#   v11_est_svd_multi <- v11_est_svd_multi / norm_vec(v11_est_svd_multi)
#
#   # 2nd PC, 1st variable
#   v12_est_svd_multi <- SimulTest_svd_multi$PCFunctions[[2]][[1]]
#   v12_est_svd_multi <- v12_est_svd_multi / norm_vec(v12_est_svd_multi)
#
#   # 1st PC, 2nd variable
#   v21_est_svd_multi <- SimulTest_svd_multi$PCFunctions[[1]][[2]]
#   v21_est_svd_multi <- v21_est_svd_multi / norm_vec(v21_est_svd_multi)
#
#   # 2nd PC, 2nd variable
#   v22_est_svd_multi <- SimulTest_svd_multi$PCFunctions[[2]][[2]]
#   v22_est_svd_multi <- v22_est_svd_multi / norm_vec(v22_est_svd_multi)
#
#
#   # Align signs
#   if (mean((v11 - v11_est_svd_multi)^2) > mean((v11 + v11_est_svd_multi)^2)) v11_est_svd_multi <- -v11_est_svd_multi
#   if (mean((v12 - v12_est_svd_multi)^2) > mean((v12 + v12_est_svd_multi)^2)) v12_est_svd_multi <- -v12_est_svd_multi
#   if (mean((v21 - v21_est_svd_multi)^2) > mean((v21 + v21_est_svd_multi)^2)) v21_est_svd_multi <- -v21_est_svd_multi
#   if (mean((v22 - v22_est_svd_multi)^2) > mean((v22 + v22_est_svd_multi)^2)) v22_est_svd_multi <- -v22_est_svd_multi
#
#   # X1
#   mse11_svd_multi <- mse_L2(v11, v11_est_svd_multi)
#   mse12_svd_multi <- mse_L2(v12, v12_est_svd_multi)
#   mse_multivariate_svdX1 <- (mse11_svd_multi + mse12_svd_multi) / 2
#
#   # X2
#   mse21_svd_multi <- mse_L2(v21, v21_est_svd_multi)
#   mse22_svd_multi <- mse_L2(v22, v22_est_svd_multi)
#   mse_multivariate_svdX2 <- (mse21_svd_multi + mse22_svd_multi) / 2
#
#   # Final average
#   mse_multivariate_svd <- (mse_multivariate_svdX1 + mse_multivariate_svdX2) / 2
#
#
#
#   ###### Smooth Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#   print("Multivariate Smooth")
#   SimulTest_smooth_multi <- ReMPCA(hd = Xobj,
#                                    centerhds = TRUE,
#                                    num_pcs = 2,
#                                    nfolds_u = 5,
#                                    nfolds_v = NULL,
#                                    thresh = 1e-10,
#                                    maxit = 100,
#                                    tuning_iter = 1,
#                                    parallel = FALSE,
#                                    weights = 0,
#                                    smoothness_type = "Second_order",
#                                    sparse_tuning_type = "SCAD",
#                                    tuning_order = "Sparsity",
#                                    cv.pick = "min",
#                                    sparse_tuning_u = NULL,
#                                    sparse_tuning_v = NULL,
#                                    smooth_tuning_u = NULL,
#                                    smooth_tuning_v = NULL)
#
#
#   # Results
#   v11_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[1]][[1]]
#   v11_est_sm_multi <- v11_est_sm_multi/norm_vec(v11_est_sm_multi)
#   v12_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[2]][[1]]
#   v12_est_sm_multi <- v12_est_sm_multi/norm_vec(v12_est_sm_multi)
#
#   v21_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[1]][[2]]
#   v21_est_sm_multi <- v21_est_sm_multi/norm_vec(v21_est_sm_multi)
#   v22_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[2]][[2]]
#   v22_est_sm_multi <- v22_est_sm_multi/norm_vec(v22_est_sm_multi)
#
#   # Align signs
#   if (mean((v11 - v11_est_sm_multi)^2) > mean((v11 + v11_est_sm_multi)^2)) v11_est_sm_multi <- -v11_est_sm_multi
#   if (mean((v12 - v12_est_sm_multi)^2) > mean((v12 + v12_est_sm_multi)^2)) v12_est_sm_multi <- -v12_est_sm_multi
#   if (mean((v21 - v21_est_sm_multi)^2) > mean((v21 + v21_est_sm_multi)^2)) v21_est_sm_multi <- -v21_est_sm_multi
#   if (mean((v22 - v22_est_sm_multi)^2) > mean((v22 + v22_est_sm_multi)^2)) v22_est_sm_multi <- -v22_est_sm_multi
#
#   # X1
#   mse11_sm_multi <- mse_L2(v11, v11_est_sm_multi)
#   mse12_sm_multi <- mse_L2(v12, v12_est_sm_multi)
#   mse_multivariate_smX1 <- (mse11_sm_multi + mse12_sm_multi) / 2
#
#   # X2
#   mse21_sm_multi <- mse_L2(v21, v21_est_sm_multi)
#   mse22_sm_multi <- mse_L2(v22, v22_est_sm_multi)
#   mse_multivariate_smX2 <- (mse21_sm_multi + mse22_sm_multi) / 2
#
#   # Final average
#   mse_multivariate_sm <- (mse_multivariate_smX1 + mse_multivariate_smX2) / 2
#
#
#   ###### Sparse Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = c(0,15,30,35,50,60,69))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = c(0,15,30,35,50,60,69))
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#
#   print("Multivariate Sparse")
#   SimulTest_spaese_multi <- ReMPCA(hd = Xobj,
#                                    centerhds = TRUE,
#                                    num_pcs = 2,
#                                    nfolds_u = 5,
#                                    nfolds_v = NULL,
#                                    thresh = 1e-10,
#                                    maxit = 100,
#                                    tuning_iter = 1,
#                                    parallel = FALSE,
#                                    weights = 0,
#                                    smoothness_type = "Second_order",
#                                    sparse_tuning_type = "SCAD",
#                                    tuning_order = "Sparsity",
#                                    cv.pick = "min",
#                                    sparse_tuning_u = NULL,
#                                    sparse_tuning_v = NULL,
#                                    smooth_tuning_u = NULL,
#                                    smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[1]][[1]]
#   v11_est_sp_multi <- v11_est_sp_multi/norm_vec(v11_est_sp_multi)
#   v12_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[2]][[1]]
#   v12_est_sp_multi <- v12_est_sp_multi/norm_vec(v12_est_sp_multi)
#
#   v21_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[1]][[2]]
#   v21_est_sp_multi <- v21_est_sp_multi/norm_vec(v21_est_sp_multi)
#   v22_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[2]][[2]]
#   v22_est_sp_multi <- v22_est_sp_multi/norm_vec(v22_est_sp_multi)
#
#   # Align signs
#   if (mean((v11 - v11_est_sp_multi)^2) > mean((v11 + v11_est_sp_multi)^2)) v11_est_sp_multi <- -v11_est_sp_multi
#   if (mean((v12 - v12_est_sp_multi)^2) > mean((v12 + v12_est_sp_multi)^2)) v12_est_sp_multi <- -v12_est_sp_multi
#   if (mean((v21 - v21_est_sp_multi)^2) > mean((v21 + v21_est_sp_multi)^2)) v21_est_sp_multi <- -v21_est_sp_multi
#   if (mean((v22 - v22_est_sp_multi)^2) > mean((v22 + v22_est_sp_multi)^2)) v22_est_sp_multi <- -v22_est_sp_multi
#
#   # X1
#   mse11_sp_multi <- mse_L2(v11, v11_est_sp_multi)
#   mse12_sp_multi <- mse_L2(v12, v12_est_sp_multi)
#   mse_multivariate_spX1 <- (mse11_sp_multi + mse12_sp_multi) / 2
#
#   # X2
#   mse21_sp_multi <- mse_L2(v21, v21_est_sp_multi)
#   mse22_sp_multi <- mse_L2(v22, v22_est_sp_multi)
#   mse_multivariate_spX2 <- (mse21_sp_multi + mse22_sp_multi) / 2
#
#   # Final average
#   mse_multivariate_sp <- (mse_multivariate_spX1 + mse_multivariate_spX2) / 2
#
#
#
#   ###### Smooth + Sparse Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = c(0,15,30,35,50,60,69)) #round(seq(0,N-1, length.out = round(N/5, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = c(0,15,30,35,50,60,69)) #round(seq(0,N-1, length.out = round(N/5, digits = 0)), digits = 0))
#
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#
#   print("Multivariate Smooth and Sparse")
#   SimulTest_smsp_multi <- ReMPCA(hd = Xobj,
#                                  centerhds = TRUE,
#                                  num_pcs = 2,
#                                  nfolds_u = 5,
#                                  nfolds_v = NULL,
#                                  thresh = 1e-10,
#                                  maxit = 100,
#                                  tuning_iter = 1,
#                                  parallel = FALSE,
#                                  weights = 0,
#                                  smoothness_type = "Second_order",
#                                  sparse_tuning_type = "SCAD",
#                                  tuning_order = "Sparsity",
#                                  cv.pick = "min",
#                                  sparse_tuning_u = NULL,
#                                  sparse_tuning_v = NULL,
#                                  smooth_tuning_u = NULL,
#                                  smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[1]][[1]]
#   v11_est_ss_multi <- v11_est_ss_multi/norm_vec(v11_est_ss_multi)
#   v12_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[2]][[1]]
#   v12_est_ss_multi <- v12_est_ss_multi/norm_vec(v12_est_ss_multi)
#
#   v21_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[1]][[2]]
#   v21_est_ss_multi <- v21_est_ss_multi/norm_vec(v21_est_ss_multi)
#   v22_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[2]][[2]]
#   v22_est_ss_multi <- v22_est_ss_multi/norm_vec(v22_est_ss_multi)
#
#   # Align signs
#   if (mean((v11 - v11_est_ss_multi)^2) > mean((v11 + v11_est_ss_multi)^2)) v11_est_ss_multi <- -v11_est_ss_multi
#   if (mean((v12 - v12_est_ss_multi)^2) > mean((v12 + v12_est_ss_multi)^2)) v12_est_ss_multi <- -v12_est_ss_multi
#   if (mean((v21 - v21_est_ss_multi)^2) > mean((v21 + v21_est_ss_multi)^2)) v21_est_ss_multi <- -v21_est_ss_multi
#   if (mean((v22 - v22_est_ss_multi)^2) > mean((v22 + v22_est_ss_multi)^2)) v22_est_ss_multi <- -v22_est_ss_multi
#
#   # X1
#   mse11_ss_multi <- mse_L2(v11, v11_est_ss_multi)
#   mse12_ss_multi <- mse_L2(v12, v12_est_ss_multi)
#   mse_multivariate_ssX1 <- (mse11_ss_multi + mse12_ss_multi) / 2
#
#   # X2
#   mse21_ss_multi <- mse_L2(v21, v21_est_ss_multi)
#   mse22_ss_multi <- mse_L2(v22, v22_est_ss_multi)
#   mse_multivariate_ssX2 <- (mse21_ss_multi + mse22_ss_multi) / 2
#
#   # Final average
#   mse_multivariate_ss <- (mse_multivariate_ssX1 + mse_multivariate_ssX2) / 2
#
#
#   # assemble per‑PC MSE’s
#   mse_svd_pc    <- c(
#     PC1 = (mse11_svd_multi + mse21_svd_multi)/2,
#     PC2 = (mse12_svd_multi + mse22_svd_multi)/2
#   )
#   mse_smooth_pc <- c(
#     PC1 = (mse11_sm_multi  + mse21_sm_multi)/2,
#     PC2 = (mse12_sm_multi  + mse22_sm_multi)/2
#   )
#   mse_sparse_pc <- c(
#     PC1 = (mse11_sp_multi  + mse21_sp_multi)/2,
#     PC2 = (mse12_sp_multi  + mse22_sp_multi)/2
#   )
#   mse_ss_pc     <- c(
#     PC1 = (mse11_ss_multi  + mse21_ss_multi)/2,
#     PC2 = (mse12_ss_multi  + mse22_ss_multi)/2
#   )
#
#
#   return(list(N = N, sigma = sigma, v11 = v11, v12 = v12,
#               v21 = v21, v22 = v22, u11 = u11, u12 = u12,
#
#               mse_svd           = mse_svd_pc,
#               mse_smooth        = mse_smooth_pc,
#               mse_sparse        = mse_sparse_pc,
#               mse_smooth_sparse = mse_ss_pc,
#
#
#               v11_est_SVD = v11_est_svd_multi,
#               v12_est_SVD = v12_est_svd_multi,
#               v21_est_SVD = v21_est_svd_multi,
#               v22_est_SVD = v22_est_svd_multi,
#               mse_multivariate_svd = mse_multivariate_svd,
#
#               v11_est_smooth_multi = v11_est_sm_multi,
#               v12_est_smooth_multi = v12_est_sm_multi,
#               v21_est_smooth_multi = v21_est_sm_multi,
#               v22_est_smooth_multi = v22_est_sm_multi,
#               mse_multivariate_smooth = mse_multivariate_sm,
#
#               v11_est_sparse_multi = v11_est_sp_multi,
#               v12_est_sparse_multi = v12_est_sp_multi,
#               v21_est_sparse_multi = v21_est_sp_multi,
#               v22_est_sparse_multi = v22_est_sp_multi,
#               mse_multivariate_sparse = mse_multivariate_sp,
#
#               v11_est_smooth_sparse_multi = v11_est_ss_multi,
#               v12_est_smooth_sparse_multi = v12_est_ss_multi,
#               v21_est_smooth_sparse_multi = v21_est_ss_multi,
#               v22_est_smooth_sparse_multi = v22_est_ss_multi,
#               mse_multivariate_smooth_sparse = mse_multivariate_ss))
# }
#
# result1 <- OneWaySimulation(N = 101, sigma = 4, random_seed = 20)
# result2 <- OneWaySimulation(N = 101, sigma = 4, random_seed = 51)
#
# # Plots
# # Plot 1: PCs True vs estimated
# par(mfrow = c(2,2))
# matplot(result1$v11, type = 'l', main = "TRUE")
# matplot(result1$v12, type = 'l')
# matplot(result1$v21, type = 'l')
# matplot(result1$v22, type = 'l')
#
#
# matplot(result1$v11_est_SVD, type = 'l', main = "SVD")
# points(result1$v11, type = 'l', col = 'red')
# matplot(result1$v12_est_SVD, type = 'l')
# points(result1$v12, type = 'l', col = 'red')
# matplot(result1$v21_est_SVD, type = 'l')
# points(result1$v21, type = 'l', col = 'red')
# matplot(result1$v22_est_SVD, type = 'l')
# points(result1$v22, type = 'l', col = 'red')
#
# matplot(result1$v11_est_smooth_multi, type = 'l', main = "Smooth - Multivariate")
# points(result1$v11, type = 'l', col = 'red')
# matplot(result1$v12_est_smooth_multi, type = 'l')
# points(result1$v12, type = 'l', col = 'red')
# matplot(result1$v21_est_smooth_multi, type = 'l')
# points(result1$v21, type = 'l', col = 'red')
# matplot(result1$v22_est_smooth_multi, type = 'l')
# points(result1$v22, type = 'l', col = 'red')
#
# matplot(result1$v11_est_sparse_multi, type = 'l', main = "Sparse - Multivariate")
# points(result1$v11, type = 'l', col = 'red')
# matplot(result1$v12_est_sparse_multi, type = 'l')
# points(result1$v12, type = 'l', col = 'red')
# matplot(result1$v21_est_sparse_multi, type = 'l')
# points(result1$v21, type = 'l', col = 'red')
# matplot(result1$v22_est_sparse_multi, type = 'l')
# points(result1$v22, type = 'l', col = 'red')
#
# matplot(result1$v11_est_smooth_sparse_multi, type = 'l', main = "Sparse + Smooth - Multivariate")
# points(result1$v11, type = 'l', col = 'red')
# matplot(result1$v12_est_smooth_sparse_multi, type = 'l')
# points(result1$v12, type = 'l', col = 'red')
# matplot(result1$v21_est_smooth_sparse_multi, type = 'l')
# points(result1$v21, type = 'l', col = 'red')
# matplot(result1$v22_est_smooth_sparse_multi, type = 'l')
# points(result1$v22, type = 'l', col = 'red')
#
#
# sprintf("%.7f",result1$mse_multivariate_svd)
# sprintf("%.7f",result1$mse_multivariate_smooth)
# sprintf("%.7f",result1$mse_multivariate_sparse)
# sprintf("%.7f",result1$mse_multivariate_smooth_sparse)
#
# # Plot 2: PC scores True vs estimated
#
#
#
# # Plot 3: Box plot of MSE
# nrep <- 50
# seeds <- sample(200:500, size = 50, replace = FALSE)
# results_list1 <- lapply(seeds, function(s) {
#   set.seed(s)  # initialize RNG
#   OneWaySimulation(N = 101, sigma = 4, random_seed = s)
# })
#
# nrep <- 50
# seeds <- sample(1:199, size = 50, replace = FALSE)
# results_list <- lapply(seeds, function(s) {
#   set.seed(s)  # initialize RNG
#   OneWaySimulation(N = 101, sigma = 4, random_seed = s)
# })
#
# combined_list <- append(results_list, results_list1)
#
# save(combined_list, file = "combined_list.RData")
#
# mse_svd1    <- sapply(results_list, `[[`, "mse_multivariate_svd")
# mse_smooth1 <- sapply(results_list, `[[`, "mse_multivariate_smooth")
# mse_sparse1 <- sapply(results_list, `[[`, "mse_multivariate_sparse")
# mse_ss1     <- sapply(results_list, `[[`, "mse_multivariate_smooth_sparse")
#
# mse_svd0    <- sapply(results_list1, `[[`, "mse_multivariate_svd")
# mse_smooth0 <- sapply(results_list1, `[[`, "mse_multivariate_smooth")
# mse_sparse0 <- sapply(results_list1, `[[`, "mse_multivariate_sparse")
# mse_ss0     <- sapply(results_list1, `[[`, "mse_multivariate_smooth_sparse")
#
# mse_svd <- c(mse_svd1,mse_svd0)
# mse_smooth <- c(mse_smooth0, mse_smooth1)
# mse_sparse <- c(mse_sparse0, mse_sparse1)
# mse_smooth_sparse <- c(mse_ss0,mse_ss1)
#
# boxplot(
#   mse_svd0, mse_smooth0, mse_sparse0, mse_ss0,
#   names = c("SVD","Smooth","Sparse","Smooth+Sparse"),
#   ylab  = "Multivariate MSE"
# )
#
#
# # Plot 4: Summary of MSE
# library(dplyr)
# library(tibble)
# library(knitr)
#
# # 1) put them in a named list
# mse_list <- list(
#   SVD           = mse_svd1,
#   Smooth        = mse_smooth1,
#   Sparse        = mse_sparse1,
#   `Smooth+Sparse` = mse_smooth_sparse1
# )
#
# # 2) compute the four summaries for each method
# summary_tbl <- map_df(mse_list,
#                       ~ tibble(
#                         Q1     = quantile(.x, .25),
#                         Median = median(.x),
#                         Mean   = mean(.x),
#                         Q3     = quantile(.x, .75)
#                       ),
#                       .id = "Method"
# ) %>%
#   select(Method, Q1, Median, Mean, Q3)
#
# # 3) print as a nice table
# kable(summary_tbl,
#       digits = 4,
#       caption = "Simulation summary of multivariate MSE by method")
#
#
