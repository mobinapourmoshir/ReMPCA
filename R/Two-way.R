# library(ReMPCA)
# ############### Simulation for two functional variables - two-way ###############
# TwoWaySimulation <- function(N, sigma){
#   norm_vec <- function(x) sqrt(sum(x^2))
#
#   rescale <- function(x, to = c(0, 1)) {
#     rng <- range(x)
#     (x - rng[1]) / diff(rng) * (to[2] - to[1]) + to[1]
#   }
#
#   mse_L2 <- function(true, est) {
#     true <- true / sqrt(sum(true^2))
#     est <- est / sqrt(sum(est^2))
#     err1 <- sum((true - est)^2)
#     err2 <- sum((true + est)^2)
#     min(err1, err2) / length(true)  # approximate L2 by Riemann sum
#   }
#
#   ise <- function(true, estimated, t = NULL) {
#     # Compute squared errors
#     se <- (true - estimated)^2
#
#     # If no grid provided, assume equal spacing
#     if (is.null(t)) {
#       return(mean(se))
#     }
#
#     # Otherwise integrate using trapezoidal rule
#     dt <- diff(t)
#     midpts <- (se[-1] + se[-length(se)]) / 2
#     integral <- sum(midpts * dt)
#     return(integral)
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
#   s <- seq(-1,1,length.out = n)
#   # define two disjoint blocks
#   block1 <- s <= 0      # left half
#   block2 <- s >  0      # right half
#
#   # build u1 nonzero only on block1
#   u1 <- numeric(n)
#   u1[block1] <- sin(2*pi*s[block1])
#   u1 <- u1 / norm_vec(u1)
#   u2 <- numeric(n)
#   mask2 <- (s >= 0)
#   u2[mask2] <- sin(pi * s[mask2])         # sin π·s runs from 0→π over [0,1]
#   u2 <- u2 / norm_vec(u2)
#
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
#   Q <- cbind(u1, u2)
#   X1 <- Q %*% t(cbind(v11, v12)) + eps1
#   X2 <- Q %*% t(cbind(v21, v22)) + eps2
#
#   # combine side by side
#   X_all <- cbind(X1, X2)
#
#   ###### SVD ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                   Sparsity_parameter = 0)
#
#   print("SVD")
#   SimulTest_SVD <- ReMPCA(hd = Xobj,
#                           centerhds = TRUE,
#                           num_pcs = 2,
#                           nfolds_u = 5,
#                           nfolds_v = NULL,
#                           thresh = 1e-10,
#                           maxit = 100,
#                           tuning_iter = 1,
#                           parallel = FALSE,
#                           weights = 0,
#                           smoothness_type = "Second_order",
#                           sparse_tuning_type = "soft",
#                           tuning_order = "Sparsity",
#                           cv.pick = "min",
#                           sparse_tuning_u = NULL,
#                           sparse_tuning_v = NULL,
#                           smooth_tuning_u = NULL,
#                           smooth_tuning_v = NULL)
#
#   # Results
#   algo <- SimulTest_SVD
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_SVD <- u1_est; u2_est_SVD  <- u2_est
#   v11_est_SVD  <- v11_est; v12_est_SVD  <- v12_est
#   v21_est_SVD  <- v21_est; v22_est_SVD <- v22_est
#
#   # MSE
#   mse11_SVD  <- mse_L2(v11, v11_est)
#   mse12_SVD  <- mse_L2(v12, v12_est)
#   mse21_SVD  <- mse_L2(v21, v21_est)
#   mse22_SVD  <- mse_L2(v22, v22_est)
#   mse_V_SVD  <- (mse11_SVD  + mse11_SVD  + mse21_SVD  + mse22_SVD ) / 4
#
#   mse_u1_SVD  <- mse_L2(u1, u1_est)
#   mse_u2_SVD  <- mse_L2(u2, u2_est)
#   mse_U_SVD  <- (mse_u1_SVD  + mse_u2_SVD) / 2
#
#   # ISE
#   ISE_u1_SVD  <- ise(u1, u1_est)
#   ISE_u2_SVD  <- ise(u2, u2_est)
#   ISE_v11_SVD  <- ise(v11, v11_est)
#   ISE_v12_SVD <- ise(v12, v12_est)
#   ISE_v21_SVD <- ise(v21, v21_est)
#   ISE_v22_SVD <- ise(v22, v22_est)
#
#
#   ###### Smoothness on u only ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
#                            Sparsity_parameter = 0)
#
#   print("Smoothness on u")
#   SimulTest_sm_u <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 2,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Sparsity",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = NULL)
#
#
#   # Results
#   algo <- SimulTest_sm_u
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_sm_u <- u1_est; u2_est_sm_u <- u2_est
#   v11_est_sm_u <- v11_est; v12_est_sm_u <- v12_est
#   v21_est_sm_u <- v21_est; v22_est_sm_u <- v22_est
#
#   # MSE
#   mse11_sm_u <- mse_L2(v11, v11_est)
#   mse12_sm_u <- mse_L2(v12, v12_est)
#   mse21_sm_u <- mse_L2(v21, v21_est)
#   mse22_sm_u <- mse_L2(v22, v22_est)
#   mse_V_sm_u <- (mse11_sm_u + mse11_sm_u + mse21_sm_u + mse22_sm_u) / 4
#
#   mse_u1_sm_u <- mse_L2(u1, u1_est)
#   mse_u2_sm_u <- mse_L2(u2, u2_est)
#   mse_U_sm_u <- (mse_u1_sm_u + mse_u2_sm_u) / 2
#
#   # ISE
#   ISE_u1_sm_u <- ise(u1, u1_est)
#   ISE_u2_sm_u <- ise(u2, u2_est)
#   ISE_v11_sm_u <- ise(v11, v11_est)
#   ISE_v12_sm_u <- ise(v12, v12_est)
#   ISE_v21_sm_u <- ise(v21, v21_est)
#   ISE_v22_sm_u <- ise(v22, v22_est)
#
#
#   # Ratio Relative to baseline (SVD)
#   R_ISE_u1_sm_u <- ISE_u1_sm_u/ISE_u1_SVD
#   R_ISE_u2_sm_u <- ISE_u2_sm_u/ISE_u2_SVD
#   R_ISE_v11_sm_u <- ISE_v11_sm_u/ISE_v11_SVD
#   R_ISE_v12_sm_u <- ISE_v12_sm_u/ISE_v12_SVD
#   R_ISE_v21_sm_u <- ISE_v21_sm_u/ISE_v21_SVD
#   R_ISE_v22_sm_u <- ISE_v22_sm_u/ISE_v22_SVD
#
#
#   ###### Sparsity and Smoothness on u only ######
#   X1obj <- fdClass(X1, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = 0,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
#                   Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   print("Smoothness & Sparsity on u")
#   SimulTest_ss_u <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 2,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Sparsity",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = NULL)
#
#
#   # Results
#   algo <- SimulTest_ss_u
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_ss_u <- u1_est; u2_est_ss_u <- u2_est
#   v11_est_ss_u <- v11_est; v12_est_ss_u <- v12_est
#   v21_est_ss_u <- v21_est; v22_est_ss_u <- v22_est
#
#   # MSE
#   mse11_ss_u <- mse_L2(v11, v11_est)
#   mse12_ss_u <- mse_L2(v12, v12_est)
#   mse21_ss_u <- mse_L2(v21, v21_est)
#   mse22_ss_u <- mse_L2(v22, v22_est)
#   mse_V_ss_u <- (mse11_ss_u + mse11_ss_u + mse21_ss_u + mse22_ss_u) / 4
#
#   mse_u1_ss_u <- mse_L2(u1, u1_est)
#   mse_u2_ss_u <- mse_L2(u2, u2_est)
#   mse_U_ss_u <- (mse_u1_ss_u + mse_u2_ss_u) / 2
#
#   # ISE
#   ISE_u1_ss_u <- ise(u1, u1_est)
#   ISE_u2_ss_u <- ise(u2, u2_est)
#   ISE_v11_ss_u <- ise(v11, v11_est)
#   ISE_v12_ss_u <- ise(v12, v12_est)
#   ISE_v21_ss_u <- ise(v21, v21_est)
#   ISE_v22_ss_u <- ise(v22, v22_est)
#
#   # Ratio Relative to baseline (SVD)
#   R_ISE_u1_ss_u <- ISE_u1_ss_u/ISE_u1_SVD
#   R_ISE_u2_ss_u <- ISE_u2_ss_u/ISE_u2_SVD
#   R_ISE_v11_ss_u <- ISE_v11_ss_u/ISE_v11_SVD
#   R_ISE_v12_ss_u <- ISE_v12_ss_u/ISE_v12_SVD
#   R_ISE_v21_ss_u <- ISE_v21_ss_u/ISE_v21_SVD
#   R_ISE_v22_ss_u <- ISE_v22_ss_u/ISE_v22_SVD
#
#
#   ###### Smoothness on v only ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                   Sparsity_parameter = 0)
#
#   print("Smoothness on v")
#   SimulTest_sm_v <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 2,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Sparsity",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = NULL)
#
#
#   # Results
#   algo <- SimulTest_sm_v
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_sm_v <- u1_est; u2_est_sm_v <- u2_est
#   v11_est_sm_v <- v11_est; v12_est_sm_v <- v12_est
#   v21_est_sm_v <- v21_est; v22_est_sm_v <- v22_est
#
#   # MSE
#   mse11_sm_v <- mse_L2(v11, v11_est)
#   mse12_sm_v <- mse_L2(v12, v12_est)
#   mse21_sm_v <- mse_L2(v21, v21_est)
#   mse22_sm_v <- mse_L2(v22, v22_est)
#   mse_V_sm_v <- (mse11_sm_v + mse11_sm_v + mse21_sm_v + mse22_sm_v) / 4
#
#   mse_u1_sm_v <- mse_L2(u1, u1_est)
#   mse_u2_sm_v <- mse_L2(u2, u2_est)
#   mse_U_sm_v <- (mse_u1_sm_v + mse_u2_sm_v) / 2
#
#   # ISE
#   ISE_u1_sm_v <- ise(u1, u1_est)
#   ISE_u2_sm_v <- ise(u2, u2_est)
#   ISE_v11_sm_v <- ise(v11, v11_est)
#   ISE_v12_sm_v <- ise(v12, v12_est)
#   ISE_v21_sm_v <- ise(v21, v21_est)
#   ISE_v22_sm_v <- ise(v22, v22_est)
#
#   # Ratio Relative to baseline (SVD)
#   R_ISE_u1_sm_v <- ISE_u1_sm_v/ISE_u1_SVD
#   R_ISE_u2_sm_v <- ISE_u2_sm_v/ISE_u2_SVD
#   R_ISE_v11_sm_v <- ISE_v11_sm_v/ISE_v11_SVD
#   R_ISE_v12_sm_v <- ISE_v12_sm_v/ISE_v12_SVD
#   R_ISE_v21_sm_v <- ISE_v21_sm_v/ISE_v21_SVD
#   R_ISE_v22_sm_v <- ISE_v22_sm_v/ISE_v22_SVD
#
#
#
#   ###### Sparsity and Smoothness on v only ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                   Sparsity_parameter = 0)
#
#   print("Smoothness & Sparsity on v")
#   SimulTest_ss_v <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 2,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Sparsity",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = NULL)
#
#
#   # Results
#   algo <- SimulTest_ss_v
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_ss_v <- u1_est; u2_est_ss_v <- u2_est
#   v11_est_ss_v <- v11_est; v12_est_ss_v <- v12_est
#   v21_est_ss_v <- v21_est; v22_est_ss_v <- v22_est
#
#   # MSE
#   mse11_ss_v <- mse_L2(v11, v11_est)
#   mse12_ss_v <- mse_L2(v12, v12_est)
#   mse21_ss_v <- mse_L2(v21, v21_est)
#   mse22_ss_v <- mse_L2(v22, v22_est)
#   mse_V_ss_v <- (mse11_ss_v + mse11_ss_v + mse21_ss_v + mse22_ss_v) / 4
#
#   mse_u1_ss_v <- mse_L2(u1, u1_est)
#   mse_u2_ss_v <- mse_L2(u2, u2_est)
#   mse_U_ss_v <- (mse_u1_ss_v + mse_u2_ss_v) / 2
#
#   # ISE
#   ISE_u1_ss_v <- ise(u1, u1_est)
#   ISE_u2_ss_v <- ise(u2, u2_est)
#   ISE_v11_ss_v <- ise(v11, v11_est)
#   ISE_v12_ss_v <- ise(v12, v12_est)
#   ISE_v21_ss_v <- ise(v21, v21_est)
#   ISE_v22_ss_v <- ise(v22, v22_est)
#
#   # Ratio Relative to baseline (SVD)
#   R_ISE_u1_ss_v <- ISE_u1_ss_v/ISE_u1_SVD
#   R_ISE_u2_ss_v <- ISE_u2_ss_v/ISE_u2_SVD
#   R_ISE_v11_ss_v <- ISE_v11_ss_v/ISE_v11_SVD
#   R_ISE_v12_ss_v <- ISE_v12_ss_v/ISE_v12_SVD
#   R_ISE_v21_ss_v <- ISE_v21_ss_v/ISE_v21_SVD
#   R_ISE_v22_ss_v <- ISE_v22_ss_v/ISE_v22_SVD
#
#
#   ###### Sparsity and Smoothness on u and v ######
#   X1obj <- fdClass(X1, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#   X2obj <- fdClass(X2, Smoothing_parameter = NULL,
#                    Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
#                   Sparsity_parameter = round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
#
#   print("Smoothness & Sparsity on u and v")
#   SimulTest_ss_uv <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 2,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Sparsity",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = NULL)
#
#
#   # Results
#   algo <- SimulTest_ss_uv
#
#
#   # Don't change!
#   u1_est <- algo$PCScores[,1]
#   u2_est <- algo$PCScores[,2]
#
#   v11_est <- algo$PCFunctions[[1]][[1]]
#   v11_est <- v11_est/norm_vec(v11_est)
#   v12_est <- algo$PCFunctions[[2]][[1]]
#   v12_est <- v12_est/norm_vec(v12_est)
#
#   v21_est <- algo$PCFunctions[[1]][[2]]
#   v21_est <- v21_est/norm_vec(v21_est)
#   v22_est <- algo$PCFunctions[[2]][[2]]
#   v22_est <- v22_est/norm_vec(v22_est)
#
#   # Align signs
#   if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
#   if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
#   if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
#   if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
#   if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
#   if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est
#
#   # Change
#   u1_est_ss_uv <- u1_est; u2_est_ss_uv <- u2_est
#   v11_est_ss_uv <- v11_est; v12_est_ss_uv <- v12_est
#   v21_est_ss_uv <- v21_est; v22_est_ss_uv <- v22_est
#
#   # MSE
#   mse11_ss_uv <- mse_L2(v11, v11_est)
#   mse12_ss_uv <- mse_L2(v12, v12_est)
#   mse21_ss_uv <- mse_L2(v21, v21_est)
#   mse22_ss_uv <- mse_L2(v22, v22_est)
#   mse_V_ss_uv <- (mse11_ss_uv + mse11_ss_uv + mse21_ss_uv + mse22_ss_uv) / 4
#
#   mse_u1_ss_uv <- mse_L2(u1, u1_est)
#   mse_u2_ss_uv <- mse_L2(u2, u2_est)
#   mse_U_ss_uv <- (mse_u1_ss_uv + mse_u1_ss_uv) / 2
#
#   # ISE
#   ISE_u1_ss_uv <- ise(u1, u1_est)
#   ISE_u2_ss_uv <- ise(u2, u2_est)
#   ISE_v11_ss_uv <- ise(v11, v11_est)
#   ISE_v12_ss_uv <- ise(v12, v12_est)
#   ISE_v21_ss_uv <- ise(v21, v21_est)
#   ISE_v22_ss_uv <- ise(v22, v22_est)
#
#   # Ratio Relative to baseline (SVD)
#   R_ISE_u1_ss_uv <- ISE_u1_ss_uv/ISE_u1_SVD
#   R_ISE_u2_ss_uv <- ISE_u2_ss_uv/ISE_u2_SVD
#   R_ISE_v11_ss_uv <- ISE_v11_ss_uv/ISE_v11_SVD
#   R_ISE_v12_ss_uv <- ISE_v12_ss_uv/ISE_v12_SVD
#   R_ISE_v21_ss_uv <- ISE_v21_ss_uv/ISE_v21_SVD
#   R_ISE_v22_ss_uv <- ISE_v22_ss_uv/ISE_v22_SVD
#
#
#   results.u1.table <- data.frame(param = rep("u1", 6),
#                               method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                          "Smooth v", "Smooth & Sparse v",
#                                          "Smooth & Sparse u & v"),
#                               MSE = c(mse_u1_SVD,
#                                       mse_u1_sm_u,
#                                       mse_u1_ss_u,
#                                       mse_u1_sm_v,
#                                       mse_u1_ss_v,
#                                       mse_u1_ss_uv),
#                               ISE = c(ISE_u1_SVD,
#                                       ISE_u1_sm_u,
#                                       ISE_u1_ss_u,
#                                       ISE_u1_sm_v,
#                                       ISE_u1_ss_v,
#                                       ISE_u1_ss_uv),
#                               R_ISE = c(1,
#                                         R_ISE_u1_sm_u,
#                                         R_ISE_u1_ss_u,
#                                         R_ISE_u1_sm_v,
#                                         R_ISE_u1_ss_v,
#                                         R_ISE_u1_ss_uv))
#
#   results.u2.table <- data.frame(param = rep("u2", 6),
#                                  method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                             "Smooth v", "Smooth & Sparse v",
#                                             "Smooth & Sparse u & v"),
#                                  MSE = c(mse_u2_SVD,
#                                          mse_u2_sm_u,
#                                          mse_u2_ss_u,
#                                          mse_u2_sm_v,
#                                          mse_u2_ss_v,
#                                          mse_u2_ss_uv),
#                                  ISE = c(ISE_u2_SVD,
#                                          ISE_u2_sm_u,
#                                          ISE_u2_ss_u,
#                                          ISE_u2_sm_v,
#                                          ISE_u2_ss_v,
#                                          ISE_u2_ss_uv),
#                                  R_ISE = c(1,
#                                            R_ISE_u2_sm_u,
#                                            R_ISE_u2_ss_u,
#                                            R_ISE_u2_sm_v,
#                                            R_ISE_u2_ss_v,
#                                            R_ISE_u2_ss_uv))
#
#   results.v11.table <- data.frame(param = rep("v11", 6),
#                                  method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                             "Smooth v", "Smooth & Sparse v",
#                                             "Smooth & Sparse u & v"),
#                                  MSE = c(mse11_SVD,
#                                          mse11_sm_u,
#                                          mse11_ss_u,
#                                          mse11_sm_v,
#                                          mse11_ss_v,
#                                          mse11_ss_uv),
#                                  ISE = c(ISE_v11_SVD,
#                                          ISE_v11_sm_u,
#                                          ISE_v11_ss_u,
#                                          ISE_v11_sm_v,
#                                          ISE_v11_ss_v,
#                                          ISE_v11_ss_uv),
#                                  R_ISE = c(1,
#                                            R_ISE_v11_sm_u,
#                                            R_ISE_v11_ss_u,
#                                            R_ISE_v11_sm_v,
#                                            R_ISE_v11_ss_v,
#                                            R_ISE_v11_ss_uv))
#
#   results.v12.table <- data.frame(param = rep("v12", 6),
#                                   method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                              "Smooth v", "Smooth & Sparse v",
#                                              "Smooth & Sparse u & v"),
#                                   MSE = c(mse12_SVD,
#                                           mse12_sm_u,
#                                           mse12_ss_u,
#                                           mse12_sm_v,
#                                           mse12_ss_v,
#                                           mse12_ss_uv),
#                                   ISE = c(ISE_v12_SVD,
#                                           ISE_v12_sm_u,
#                                           ISE_v12_ss_u,
#                                           ISE_v12_sm_v,
#                                           ISE_v12_ss_v,
#                                           ISE_v12_ss_uv),
#                                   R_ISE = c(1,
#                                             R_ISE_v12_sm_u,
#                                             R_ISE_v12_ss_u,
#                                             R_ISE_v12_sm_v,
#                                             R_ISE_v12_ss_v,
#                                             R_ISE_v12_ss_uv))
#
#   results.v21.table <- data.frame(param = rep("v21", 6),
#                                   method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                              "Smooth v", "Smooth & Sparse v",
#                                              "Smooth & Sparse u & v"),
#                                   MSE = c(mse21_SVD,
#                                           mse21_sm_u,
#                                           mse21_ss_u,
#                                           mse21_sm_v,
#                                           mse21_ss_v,
#                                           mse21_ss_uv),
#                                   ISE = c(ISE_v21_SVD,
#                                           ISE_v21_sm_u,
#                                           ISE_v21_ss_u,
#                                           ISE_v21_sm_v,
#                                           ISE_v21_ss_v,
#                                           ISE_v21_ss_uv),
#                                   R_ISE = c(1,
#                                             R_ISE_v21_sm_u,
#                                             R_ISE_v21_ss_u,
#                                             R_ISE_v21_sm_v,
#                                             R_ISE_v21_ss_v,
#                                             R_ISE_v21_ss_uv))
#
#   results.v22.table <- data.frame(param = rep("v22", 6),
#                                   method = c("SVD", "Smooth u", "Smooth & Sparse u",
#                                              "Smooth v", "Smooth & Sparse v",
#                                              "Smooth & Sparse u & v"),
#                                   MSE = c(mse22_SVD,
#                                           mse22_sm_u,
#                                           mse22_ss_u,
#                                           mse22_sm_v,
#                                           mse22_ss_v,
#                                           mse22_ss_uv),
#                                   ISE = c(ISE_v22_SVD,
#                                           ISE_v22_sm_u,
#                                           ISE_v22_ss_u,
#                                           ISE_v22_sm_v,
#                                           ISE_v22_ss_v,
#                                           ISE_v22_ss_uv),
#                                   R_ISE = c(1,
#                                             R_ISE_v22_sm_u,
#                                             R_ISE_v22_ss_u,
#                                             R_ISE_v22_sm_v,
#                                             R_ISE_v22_ss_v,
#                                             R_ISE_v22_ss_uv))
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
