# library(ReMPCA)
# ############### Simulation for two functional variables - one-way ###############
# OneWaySimulation <- function(N, sigma){
#   norm_vec <- function(x) sqrt(sum(x^2))
#   rescale <- function(x, to = c(0, 1)) {
#     rng <- range(x)
#     (x - rng[1]) / diff(rng) * (to[2] - to[1]) + to[1]
#   }
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
#
#
#   ######################
#   # UNIVARIATE ESTIMATION
#   ######################
#
#   ###### SVD ######
#   svd1 <- svd(X1)
#   v11_est_SVD <- svd1$v[,1]
#   v11_est_SVD <- v11_est_SVD / norm_vec(v11_est_SVD)
#   svd2 <- svd(X1)
#   v12_est_SVD <- svd2$v[,2]
#   v12_est_SVD <- v12_est_SVD / norm_vec(v12_est_SVD)
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_SVD)^2 )
#   mse11_neg <- mean( (v11 + v11_est_SVD)^2 )
#   mse11_SVD <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_SVD)^2 )
#   mse21_neg <- mean( (v12 + v12_est_SVD)^2 )
#   mse21_SVD <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_univariate_SVDX1 <- (mse11_SVD + mse21_SVD) / 2
#
#   svd1 <- svd(X2)
#   v21_est_SVD <- svd1$v[,1]
#   v21_est_SVD <- v21_est_SVD / norm_vec(v21_est_SVD)
#   svd2 <- svd(X2)
#   v22_est_SVD <- svd2$v[,2]
#   v22_est_SVD <- v22_est_SVD / norm_vec(v22_est_SVD)
#   # align sign
#   mse21_pos <- mean( (v21 - v21_est_SVD)^2 )
#   mse21_neg <- mean( (v21 + v21_est_SVD)^2 )
#   mse21_SVD <- min(mse21_pos, mse21_neg)
#   mse22_pos <- mean( (v22 - v22_est_SVD)^2 )
#   mse22_neg <- mean( (v22 + v22_est_SVD)^2 )
#   mse22_SVD <- min(mse22_pos, mse22_neg)
#   # average univariate MSE
#   mse_univariate_SVDX2 <- (mse21_SVD + mse22_SVD) / 2
#
#   mse_univariate_SVD <- (mse_univariate_SVDX1 + mse_univariate_SVDX2)/2
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
#
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
#                                      cv.pick = "1se",
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
#                                      cv.pick = "1se",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_sm_uni <- SimulTest_smoothuni_obj1$PCFunctions[[1]][[1]]
#   v12_est_sm_uni <- SimulTest_smoothuni_obj1$PCFunctions[[2]][[1]]
#
#   v11_est_sm_uniX2 <- SimulTest_smoothuni_obj2$PCFunctions[[1]][[1]]
#   v12_est_sm_uniX2 <- SimulTest_smoothuni_obj2$PCFunctions[[2]][[1]]
#
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sm_uni)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sm_uni)^2 )
#   mse11_sm_uni <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sm_uni)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sm_uni)^2 )
#   mse21_sm_uni <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_univariate_smX1 <- (mse11_sm_uni + mse21_sm_uni) / 2
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sm_uniX2)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sm_uniX2)^2 )
#   mse11_sm_uni <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sm_uniX2)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sm_uniX2)^2 )
#   mse21_sm_uni <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_univariate_smX2 <- (mse11_sm_uni + mse21_sm_uni) / 2
#
#   mse_univariate_sm <- (mse_univariate_smX1 + mse_univariate_smX2) /2
#
#   ###### Sparse Univariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#   hdobj1 <- hdClass(list(X1obj), Smoothing_parameter = 0,
#                     Sparsity_parameter = 0)
#   hdobj2 <- hdClass(list(X2obj), Smoothing_parameter = 0,
#                     Sparsity_parameter = 0)
#
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
#                                      cv.pick = "1se",
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
#                                      cv.pick = "1se",
#                                      sparse_tuning_u = NULL,
#                                      sparse_tuning_v = NULL,
#                                      smooth_tuning_u = NULL,
#                                      smooth_tuning_v = NULL)
#
#
#   # Results
#   v11_est_sp_uni <- SimulTest_sparse_uni_obj1$PCFunctions[[1]][[1]]
#   v12_est_sp_uni <- SimulTest_sparse_uni_obj1$PCFunctions[[2]][[1]]
#
#   v11_est_sp_uniX2 <- SimulTest_sparse_uni_obj2$PCFunctions[[1]][[1]]
#   v12_est_sp_uniX2 <- SimulTest_sparse_uni_obj2$PCFunctions[[2]][[1]]
#
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sp_uni)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sp_uni)^2 )
#   mse11_sp_uni <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sp_uni)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sp_uni)^2 )
#   mse21_sp_uni <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_univariate_spX1 <- (mse11_sp_uni + mse21_sp_uni) / 2
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sp_uniX2)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sp_uniX2)^2 )
#   mse11_sp_uni <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sp_uniX2)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sp_uniX2)^2 )
#   mse21_sp_uni <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_univariate_spX2 <- (mse11_sp_uni + mse21_sp_uni) / 2
#
#   mse_univariate_sp <- (mse_univariate_spX1 + mse_univariate_spX2) /2
#
#
#   ###### Smooth Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0)
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#
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
#                                    cv.pick = "1se",
#                                    sparse_tuning_u = NULL,
#                                    sparse_tuning_v = NULL,
#                                    smooth_tuning_u = NULL,
#                                    smooth_tuning_v = NULL)
#
#
#   # Results
#   v11_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[1]][[1]]
#   v12_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[2]][[1]]
#
#   v21_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[1]][[2]]
#   v22_est_sm_multi <- SimulTest_smooth_multi$PCFunctions[[2]][[2]]
#
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sm_multi)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sm_multi)^2 )
#   mse11_sm_multi <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sm_multi)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sm_multi)^2 )
#   mse21_sm_multi <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_multivariate_smX1 <- (mse11_sm_multi + mse21_sm_multi) / 2
#   # align sign
#   mse21_pos <- mean( (v21 - v21_est_sm_multi)^2 )
#   mse21_neg <- mean( (v11 + v21_est_sm_multi)^2 )
#   mse21_sm_multi <- min(mse21_pos, mse21_neg)
#   mse22_pos <- mean( (v12 - v22_est_sm_multi)^2 )
#   mse22_neg <- mean( (v12 + v22_est_sm_multi)^2 )
#   mse22_sm_multi <- min(mse22_pos, mse22_neg)
#   # average univariate MSE
#   mse_multivariate_smX2 <- (mse21_sm_multi + mse22_sm_multi) / 2
#
#   mse_multivariate_sm <- (mse_multivariate_smX1 + mse_multivariate_smX2) /2
#
#
#   ###### Sparse Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#
#
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
#                                    cv.pick = "1se",
#                                    sparse_tuning_u = NULL,
#                                    sparse_tuning_v = NULL,
#                                    smooth_tuning_u = NULL,
#                                    smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[1]][[1]]
#   v12_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[2]][[1]]
#
#   v21_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[1]][[2]]
#   v22_est_sp_multi <- SimulTest_spaese_multi$PCFunctions[[2]][[2]]
#
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_sp_multi)^2 )
#   mse11_neg <- mean( (v11 + v11_est_sp_multi)^2 )
#   mse11_sp_multi <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_sp_multi)^2 )
#   mse21_neg <- mean( (v12 + v12_est_sp_multi)^2 )
#   mse21_sp_multi <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_multivariate_spX1 <- (mse11_sp_multi + mse21_sp_multi) / 2
#   # align sign
#   mse21_pos <- mean( (v21 - v21_est_sp_multi)^2 )
#   mse21_neg <- mean( (v11 + v21_est_sp_multi)^2 )
#   mse21_sp_multi <- min(mse21_pos, mse21_neg)
#   mse22_pos <- mean( (v12 - v22_est_sp_multi)^2 )
#   mse22_neg <- mean( (v12 + v22_est_sp_multi)^2 )
#   mse22_sp_multi <- min(mse22_pos, mse22_neg)
#   # average univariate MSE
#   mse_multivariate_spX2 <- (mse21_sp_multi + mse22_sp_multi) / 2
#
#   mse_multivariate_sp <- (mse_multivariate_spX1 + mse_multivariate_spX2) /2
#
#
#
#   ###### Smooth + Sparse Multivariate ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                            Sparsity_parameter = 0)
#
#
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
#                                  cv.pick = "1se",
#                                  sparse_tuning_u = NULL,
#                                  sparse_tuning_v = NULL,
#                                  smooth_tuning_u = NULL,
#                                  smooth_tuning_v = NULL)
#
#   # Results
#   v11_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[1]][[1]]
#   v12_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[2]][[1]]
#
#   v21_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[1]][[2]]
#   v22_est_ss_multi <- SimulTest_smsp_multi$PCFunctions[[2]][[2]]
#
#   # align sign
#   mse11_pos <- mean( (v11 - v11_est_ss_multi)^2 )
#   mse11_neg <- mean( (v11 + v11_est_ss_multi)^2 )
#   mse11_ss_multi <- min(mse11_pos, mse11_neg)
#   mse21_pos <- mean( (v12 - v12_est_ss_multi)^2 )
#   mse21_neg <- mean( (v12 + v12_est_ss_multi)^2 )
#   mse21_ss_multi <- min(mse21_pos, mse21_neg)
#   # average univariate MSE
#   mse_multivariate_ssX1 <- (mse11_ss_multi + mse21_ss_multi) / 2
#   # align sign
#   mse21_pos <- mean( (v21 - v21_est_ss_multi)^2 )
#   mse21_neg <- mean( (v11 + v21_est_ss_multi)^2 )
#   mse21_ss_multi <- min(mse21_pos, mse21_neg)
#   mse22_pos <- mean( (v12 - v22_est_ss_multi)^2 )
#   mse22_neg <- mean( (v12 + v22_est_ss_multi)^2 )
#   mse22_ss_multi <- min(mse22_pos, mse22_neg)
#   # average univariate MSE
#   mse_multivariate_ssX2 <- (mse21_ss_multi + mse22_ss_multi) / 2
#
#   mse_multivariate_ss <- (mse_multivariate_ssX1 + mse_multivariate_ssX2) /2
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
#               v21_est_uni_smooth = v21_est_sm_uniX2,
#               v22_est_uni_smooth = v22_est_sm_uniX2,
#               mse_univariate_smooth = mse_univariate_sm,
#
#               v11_est_uni_sparse = v11_est_sp_uni,
#               v12_est_uni_sparse = v12_est_sp_uni,
#               v21_est_uni_sparse = v21_est_sp_uniX2,
#               v22_est_uni_sparse = v22_est_sp_uniX2,
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
#
