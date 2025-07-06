# library(ReMPCA)
# ############### Simulation for two functional variables - one-way ###############
# OneWaySimulation <- function(N, sigma){
#   norm_vec <- function(x) sqrt(sum(x^2))
#   rescale <- function(x, to = c(0, 1)) {
#     rng <- range(x)
#     (x - rng[1]) / diff(rng) * (to[2] - to[1]) + to[1]
#   }
#   mse_L2 <- function(true, est) {
#     true <- true / sqrt(sum(true^2))
#     est <- est / sqrt(sum(est^2))
#     err1 <- sum((true - est)^2)
#     err2 <- sum((true + est)^2)
#     min(err1, err2) / length(true)  # approximate L2 by Riemann sum
#   }
#
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
#   ######################
#   # UNIVARIATE ESTIMATION
#   ######################
#
#   ###### SVD ######
#
#   # X1
#   svd1 <- svd(X1)
#   v11_est_SVD <- svd1$v[,1]; v11_est_SVD <- v11_est_SVD / norm_vec(v11_est_SVD)
#   v12_est_SVD <- svd1$v[,2]; v12_est_SVD <- v12_est_SVD / norm_vec(v12_est_SVD)
#   mse11_SVD <- mse_L2(v11, v11_est_SVD)
#   mse21_SVD <- mse_L2(v12, v12_est_SVD)
#   mse_univariate_SVDX1 <- (mse11_SVD + mse21_SVD) / 2
#
#   # X2
#   svd2 <- svd(X2)
#   v21_est_SVD <- svd2$v[,1]; v21_est_SVD <- v21_est_SVD / norm_vec(v21_est_SVD)
#   v22_est_SVD <- svd2$v[,2]; v22_est_SVD <- v22_est_SVD / norm_vec(v22_est_SVD)
#   mse21_SVD <- mse_L2(v21, v21_est_SVD)
#   mse22_SVD <- mse_L2(v22, v22_est_SVD)
#   mse_univariate_SVDX2 <- (mse21_SVD + mse22_SVD) / 2
#
#   # Final average
#   mse_univariate_SVD <- (mse_univariate_SVDX1 + mse_univariate_SVDX2) / 2
#
#   # MRAE
#   # Reconstruct X1 and X2 using top 2 PCs
#   X1reconstructed_SVD <- svd1$d[1] * svd1$u[,1] %*% t(v11_est_SVD) +
#     svd1$d[2] * svd1$u[,2] %*% t(v12_est_SVD)
#
#   X2reconstructed_SVD <- svd2$d[1] * svd2$u[,1] %*% t(v21_est_SVD) +
#     svd2$d[2] * svd2$u[,2] %*% t(v22_est_SVD)
#
#   # Compute MRAE for X1
#   numerator_X1 <- sqrt(rowSums((X1 - X1reconstructed_SVD)^2))
#   denominator_X1 <- sqrt(rowSums(X1^2))
#   mrae_X1 <- mean(numerator_X1 / denominator_X1)
#
#   # Compute MRAE for X2
#   numerator_X2 <- sqrt(rowSums((X2 - X2reconstructed_SVD)^2))
#   denominator_X2 <- sqrt(rowSums(X2^2))
#   mrae_X2 <- mean(numerator_X2 / denominator_X2)
#
#   # Final average MRAE over both datasets
#   mrae_total <- (mrae_X1 + mrae_X2) / 2
#
#
#   ###### Smooth Univariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   hdobj1 <- hdClass(list(X1obj), Smoothing_parameter = 0,
#                    Sparsity_parameter = 0)
#   hdobj2 <- hdClass(list(X2obj), Smoothing_parameter = 0,
#                     Sparsity_parameter = 0)
#   print("Univariate Smooth")
#   SimulTest_smoothuni_obj1 <- ReMPCA(hd = hdobj1,
#                                      centerhds = TRUE,
#                                      num_pcs = 2,
#                                      nfolds_u = 5,
#                                      nfolds_v = NULL,
#                                      thresh = 1e-10,
#                                      maxit = 100,
#                                      tuning_iter = 1,
#                                      parallel = FALSE,
#                                      weights = 0,
#                                      smoothness_type = "Second_order",
#                                      sparse_tuning_type = "soft",
#                                      tuning_order = "Sparsity",
#                                      cv.pick = "min",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#
#   SimulTest_smoothuni_obj2 <- ReMPCA(hd = hdobj2,
#                                      centerhds = TRUE,
#                                      num_pcs = 2,
#                                      nfolds_u = 5,
#                                      nfolds_v = NULL,
#                                      thresh = 1e-10,
#                                      maxit = 100,
#                                      tuning_iter = 1,
#                                      parallel = FALSE,
#                                      weights = 0,
#                                      smoothness_type = "Second_order",
#                                      sparse_tuning_type = "soft",
#                                      tuning_order = "Sparsity",
#                                      cv.pick = "min",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_sm_uni <- SimulTest_smoothuni_obj1$PCFunctions[[1]][[1]]
#   v11_est_sm_uni <-  v11_est_sm_uni/ norm_vec(v11_est_sm_uni)
#   v12_est_sm_uni <- SimulTest_smoothuni_obj1$PCFunctions[[2]][[1]]
#   v12_est_sm_uni <-  v12_est_sm_uni/ norm_vec(v12_est_sm_uni)
#
#   v11_est_sm_uniX2 <- SimulTest_smoothuni_obj2$PCFunctions[[1]][[1]]
#   v11_est_sm_uniX2 <-  v11_est_sm_uniX2/ norm_vec(v11_est_sm_uniX2)
#   v12_est_sm_uniX2 <- SimulTest_smoothuni_obj2$PCFunctions[[2]][[1]]
#   v12_est_sm_uniX2 <- v12_est_sm_uniX2/ norm_vec(v12_est_sm_uniX2)
#
#   # Align signs
#   if (mean((v11 - v11_est_sm_uni)^2) > mean((v11 + v11_est_sm_uni)^2)) v11_est_sm_uni <- -v11_est_sm_uni
#   if (mean((v12 - v12_est_sm_uni)^2) > mean((v12 + v12_est_sm_uni)^2)) v12_est_sm_uni <- -v12_est_sm_uni
#   if (mean((v21 - v11_est_sm_uniX2)^2) > mean((v21 + v11_est_sm_uniX2)^2)) v11_est_sm_uniX2 <- -v11_est_sm_uniX2
#   if (mean((v22 - v12_est_sm_uniX2)^2) > mean((v22 + v12_est_sm_uniX2)^2)) v12_est_sm_uniX2 <- -v12_est_sm_uniX2
#
#   # X1
#   mse11_sm_uni <- mse_L2(v11, v11_est_sm_uni)
#   mse21_sm_uni <- mse_L2(v12, v12_est_sm_uni)
#   mse_univariate_smX1 <- (mse11_sm_uni + mse21_sm_uni) / 2
#
#   # X2
#   mse11_sm_uni <- mse_L2(v21, v11_est_sm_uniX2)
#   mse21_sm_uni <- mse_L2(v22, v12_est_sm_uniX2)
#   mse_univariate_smX2 <- (mse11_sm_uni + mse21_sm_uni) / 2
#
#   # Final average
#   mse_univariate_sm <- (mse_univariate_smX1 + mse_univariate_smX2) / 2
#
#   # # MRAE
#   # # Reconstruct X1 and X2 using top 2 PCs
#   # X1reconstructed_sm_uni <- svd1$d[1] * svd1$u[,1] %*% t(v11_est_sm_uni) +
#   #   svd1$d[2] * svd1$u[,2] %*% t(v12_est_sm_uni)
#   #
#   # X2reconstructed_sm_uni <- svd2$d[1] * svd2$u[,1] %*% t(v11_est_sm_uniX2) +
#   #   svd2$d[2] * svd2$u[,2] %*% t(v12_est_sm_uniX2)
#   #
#   # # Compute MRAE for X1
#   # numerator_X1 <- sqrt(rowSums((X1 - X1reconstructed_sm_uni)^2))
#   # denominator_X1 <- sqrt(rowSums(X1^2))
#   # mrae_X1 <- mean(numerator_X1 / denominator_X1)
#   #
#   # # Compute MRAE for X2
#   # numerator_X2 <- sqrt(rowSums((X2 - X2reconstructed_sm_uni)^2))
#   # denominator_X2 <- sqrt(rowSums(X2^2))
#   # mrae_X2 <- mean(numerator_X2 / denominator_X2)
#   #
#   # # Final average MRAE over both datasets
#   # mrae_total <- (mrae_X1 + mrae_X2) / 2
#
#   ###### Sparse Univariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   hdobj1 <- hdClass(list(X1obj), Smoothing_parameter = 0,
#                     Sparsity_parameter = 0)
#   hdobj2 <- hdClass(list(X2obj), Smoothing_parameter = 0,
#                     Sparsity_parameter = 0)
#   print("Univariate Sparse")
#   SimulTest_sparse_uni_obj1 <- ReMPCA(hd = hdobj1,
#                                      centerhds = TRUE,
#                                      num_pcs = 2,
#                                      nfolds_u = 5,
#                                      nfolds_v = NULL,
#                                      thresh = 1e-10,
#                                      maxit = 100,
#                                      tuning_iter = 1,
#                                      parallel = FALSE,
#                                      weights = 0,
#                                      smoothness_type = "Second_order",
#                                      sparse_tuning_type = "soft",
#                                      tuning_order = "Sparsity",
#                                      cv.pick = "min",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#
#   SimulTest_sparse_uni_obj2 <- ReMPCA(hd = hdobj2,
#                                      centerhds = TRUE,
#                                      num_pcs = 2,
#                                      nfolds_u = 5,
#                                      nfolds_v = NULL,
#                                      thresh = 1e-10,
#                                      maxit = 100,
#                                      tuning_iter = 1,
#                                      parallel = FALSE,
#                                      weights = 0,
#                                      smoothness_type = "Second_order",
#                                      sparse_tuning_type = "soft",
#                                      tuning_order = "Sparsity",
#                                      cv.pick = "min",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_sp_uni <- SimulTest_sparse_uni_obj1$PCFunctions[[1]][[1]]
#   v11_est_sp_uni <- v11_est_sp_uni/ norm_vec(v11_est_sp_uni)
#   v12_est_sp_uni <- SimulTest_sparse_uni_obj1$PCFunctions[[2]][[1]]
#   v12_est_sp_uni <- v12_est_sp_uni/ norm_vec(v12_est_sp_uni)
#
#   v11_est_sp_uniX2 <- SimulTest_sparse_uni_obj2$PCFunctions[[1]][[1]]
#   v11_est_sp_uniX2 <- v11_est_sp_uniX2/ norm_vec(v11_est_sp_uniX2)
#   v12_est_sp_uniX2 <- SimulTest_sparse_uni_obj2$PCFunctions[[2]][[1]]
#   v12_est_sp_uniX2 <- v12_est_sp_uniX2/ norm_vec(v12_est_sp_uniX2)
#
#   # Align signs
#   if (mean((v11 - v11_est_sp_uni)^2) > mean((v11 + v11_est_sp_uni)^2)) v11_est_sp_uni <- -v11_est_sp_uni
#   if (mean((v12 - v12_est_sp_uni)^2) > mean((v12 + v12_est_sp_uni)^2)) v12_est_sp_uni <- -v12_est_sp_uni
#   if (mean((v21 - v11_est_sp_uniX2)^2) > mean((v21 + v11_est_sp_uniX2)^2)) v11_est_sp_uniX2 <- -v11_est_sp_uniX2
#   if (mean((v22 - v12_est_sp_uniX2)^2) > mean((v22 + v12_est_sp_uniX2)^2)) v12_est_sp_uniX2 <- -v12_est_sp_uniX2
#
#   # X1
#   mse11_sm_uni <- mse_L2(v11, v11_est_sp_uni)
#   mse21_sm_uni <- mse_L2(v12, v12_est_sp_uni)
#   mse_univariate_spX1 <- (mse11_sm_uni + mse21_sm_uni) / 2
#
#   # X2
#   mse11_sm_uni <- mse_L2(v21, v11_est_sp_uniX2)
#   mse21_sm_uni <- mse_L2(v22, v12_est_sp_uniX2)
#   mse_univariate_spX2 <- (mse11_sm_uni + mse21_sm_uni) / 2
#
#   # Final average
#   mse_univariate_sp <- (mse_univariate_spX1 + mse_univariate_spX2) / 2
#
#
#   ##########################
#   # Multivariate ESTIMATION
#   ##########################
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
#                                    sparse_tuning_type = "soft",
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
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
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
#                                    sparse_tuning_type = "soft",
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
#   mse21_sp_multi <- mse_L2(v12, v12_est_sp_multi)
#   mse_multivariate_spX1 <- (mse11_sp_multi + mse21_sp_multi) / 2
#
#   # X2
#   mse11_sp_multi <- mse_L2(v21, v21_est_sp_multi)
#   mse21_sp_multi <- mse_L2(v22, v22_est_sp_multi)
#   mse_multivariate_spX2 <- (mse11_sp_multi + mse21_sp_multi) / 2
#
#   # Final average
#   mse_multivariate_sp <- (mse_multivariate_spX1 + mse_multivariate_spX2) / 2
#
#
#
#   ###### Smooth + Sparse Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
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
#                                  sparse_tuning_type = "soft",
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
#   mse21_ss_multi <- mse_L2(v12, v12_est_ss_multi)
#   mse_multivariate_ssX1 <- (mse11_ss_multi + mse21_ss_multi) / 2
#
#   # X2
#   mse11_ss_multi <- mse_L2(v21, v21_est_ss_multi)
#   mse21_ss_multi <- mse_L2(v22, v22_est_ss_multi)
#   mse_multivariate_ssX2 <- (mse11_ss_multi + mse21_ss_multi) / 2
#
#   # Final average
#   mse_multivariate_ss <- (mse_multivariate_ssX1 + mse_multivariate_ssX2) / 2
#
#
#
#   return(list(N = N, sigma = sigma, v11 = v11, v12 = v12,
#               v21 = v21, v22 = v22, u11 = u11, u12 = u12,
#
#               v11_est_SVD = v11_est_SVD,
#               v12_est_SVD = v12_est_SVD,
#               v21_est_SVD = v21_est_SVD,
#               v22_est_SVD = v22_est_SVD,
#               mse_univariate_SVD = mse_univariate_SVD,
#
#               v11_est_uni_smooth = v11_est_sm_uni,
#               v12_est_uni_smooth = v12_est_sm_uni,
#               v21_est_uni_smooth = v11_est_sm_uniX2,
#               v22_est_uni_smooth = v12_est_sm_uniX2,
#               mse_univariate_smooth = mse_univariate_sm,
#
#               v11_est_uni_sparse = v11_est_sp_uni,
#               v12_est_uni_sparse = v12_est_sp_uni,
#               v21_est_uni_sparse = v11_est_sp_uniX2,
#               v22_est_uni_sparse = v12_est_sp_uniX2,
#               mse_univariate_sparse = mse_univariate_sp,
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
# result1 <- OneWaySimulation(N = 200, sigma = 3)
#
# # par(mfrow = c(2,2))
# # matplot(v11, type = 'l', main = "TRUE")
# # matplot(v12, type = 'l')
# # matplot(v21, type = 'l')
# # matplot(v22, type = 'l')
# #
# #
# # matplot(v11_est_SVD, type = 'l', main = "SVD")
# # points(v11, type = 'l', col = 'red')
# # matplot(v12_est_SVD, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(v21_est_SVD, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(v22_est_SVD, type = 'l')
# # points(v22, type = 'l', col = 'red')
# #
# # matplot(result1$v11_est_uni_smooth, type = 'l', main = "Smooth - univariate")
# # points(v11, type = 'l', col = 'red')
# # matplot(result1$v12_est_uni_smooth, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(result1$v21_est_uni_smooth, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(result1$v22_est_uni_smooth, type = 'l')
# # points(v22, type = 'l', col = 'red')
# #
# # matplot(result1$v11_est_uni_sparse, type = 'l', main = "Sparse - Univariate")
# # points(v11, type = 'l', col = 'red')
# # matplot(result1$v12_est_uni_sparse, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(result1$v21_est_uni_sparse, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(result1$v22_est_uni_sparse, type = 'l')
# # points(v22, type = 'l', col = 'red')
# #
# # matplot(result1$v11_est_smooth_multi, type = 'l', main = "Smooth - Multivariate")
# # points(v11, type = 'l', col = 'red')
# # matplot(result1$v12_est_smooth_multi, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(result1$v21_est_smooth_multi, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(result1$v22_est_smooth_multi, type = 'l')
# # points(v22, type = 'l', col = 'red')
# #
# # matplot(result1$v11_est_sparse_multi, type = 'l', main = "Sparse - Multivariate")
# # points(v11, type = 'l', col = 'red')
# # matplot(result1$v12_est_sparse_multi, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(result1$v21_est_sparse_multi, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(result1$v22_est_sparse_multi, type = 'l')
# # points(v22, type = 'l', col = 'red')
# #
# # matplot(result1$v11_est_smooth_sparse_multi, type = 'l', main = "Sparse + Smooth - Multivariate")
# # points(v11, type = 'l', col = 'red')
# # matplot(result1$v12_est_smooth_sparse_multi, type = 'l')
# # points(v12, type = 'l', col = 'red')
# # matplot(result1$v21_est_smooth_sparse_multi, type = 'l')
# # points(v21, type = 'l', col = 'red')
# # matplot(result1$v22_est_smooth_sparse_multi, type = 'l')
# # points(v22, type = 'l', col = 'red')
#
#
# sprintf("%.7f",result1$mse_univariate_SVD)
# sprintf("%.7f",result1$mse_univariate_smooth)
# sprintf("%.7f",result1$mse_univariate_sparse)
# sprintf("%.7f",result1$mse_multivariate_smooth)
# sprintf("%.7f",result1$mse_multivariate_sparse)
# sprintf("%.7f",result1$mse_multivariate_smooth_sparse)
#
