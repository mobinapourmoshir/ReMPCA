library(ReMPCA)
#N = 101; sigma = 0.05; random_seed =20
############### Simulation for two functional variables - two-way ###############
TwoWaySimulation <- function(N, sigma,random_seed){
  norm_vec <- function(x) sqrt(sum(x^2))

  rescale <- function(x, to = c(0, 1)) {
    rng <- range(x)
    (x - rng[1]) / diff(rng) * (to[2] - to[1]) + to[1]
  }

  mse_L2 <- function(true, est) {
    true <- true / sqrt(sum(true^2))
    est <- est / sqrt(sum(est^2))
    err1 <- sum((true - est)^2)
    err2 <- sum((true + est)^2)
    min(err1, err2) / length(true)  # approximate L2 by Riemann sum
  }

  ise <- function(true, estimated, t = NULL) {
    # Compute squared errors
    se <- (true - estimated)^2

    # If no grid provided, assume equal spacing
    if (is.null(t)) {
      return(mean(se))
    }

    # Otherwise integrate using trapezoidal rule
    dt <- diff(t)
    midpts <- (se[-1] + se[-length(se)]) / 2
    integral <- sum(midpts * dt)
    return(integral)
  }

  set.seed(random_seed)
  sigma1 <- 20
  sigma2 <- 10
  sigma <- sigma
  n <- m <- N
  t <- seq(-1,1,length.out = m)

  # generate scores
  s <- seq(-1,1,length.out = n)
  # define two disjoint blocks
  block1 <- s <= 0      # left half
  block2 <- s >  0      # right half

  # build u1 nonzero only on block1
  u1 <- numeric(n)
  u1[block1] <- sin(2*pi*s[block1])
  u1[block2] <- sin(2*pi*s[block2])
  #u1 <- u1 / norm_vec(u1)
  u2 <- numeric(n)
  mask2 <- (s >= 0)
  u2[mask2] <- sin(pi * s[mask2])
  #u2 <- u2 / norm_vec(u2)


  # generate FPCs for first variable
  v11 <- t + sin(pi*t)
  v11 <- v11 / norm_vec(v11)

  v12 <- cos(3*pi*t)
  v12 <- v12 / norm_vec(v12)

  # generate FPCs for second variable
  block1 <- t <= -1/3
  block2 <- (t > -1/3) & (t < 1/3)
  block3 <- t >= 1/3

  # Define three sparse functions
  f1 <- numeric(m)
  f1[block1] <- sin(pi * rescale(t[block1], to = c(0, 1)))
  f1[block3] <- sin(2 * pi * rescale(t[block3], to = c(0, 1)))

  f2 <- numeric(m)
  f2[block1] <- sin(2 * pi * rescale(t[block1], to = c(0, 1)))
  f2[block3] <- sin(pi * rescale(t[block3], to = c(0, 1)))

  f3 <- numeric(m)
  f3[block2] <- sin(3 * pi * rescale(t[block2], to = c(0, 1)))

  # Normalize
  v21 <- f3 / sqrt(sum(f3^2))
  v21 <- v21/norm_vec(v21)
  v22 <- f2 / sqrt(sum(f2^2))
  v22 <- v22/ norm_vec(v22)

  # generate noise
  eps1 <- matrix(rnorm(n*m, 0, sigma), n, m)
  eps2 <- matrix(rnorm(n*m, 0, sigma), n, m)

  # data matrices
  Q <- cbind(u1, u2)
  X1 <- u1 %*% t(v11) + u2 %*% t(v12) + eps1
  X2 <- u1 %*% t(v21) + u2 %*% t(v22) + eps2

  # combine side by side
  X_all <- cbind(X1, X2)

  ###### Sparsity and Smoothness on u and v ######
  X1obj <- fdClass(X1, Smoothing_parameter = NULL,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))
  X2obj <- fdClass(X2, Smoothing_parameter = NULL,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
                  Sparsity_parameter = c(0,15,30,35,40,42,45,50,55))

  print("Smoothness & Sparsity on u and v")
  SimulTest_ss_uv <- ReMPCA(hd = Xobj,
                            centerhds = FALSE,
                            num_pcs = 2,
                            nfolds_u = 5,
                            nfolds_v = NULL,
                            thresh = 1e-10,
                            maxit = 100,
                            tuning_iter = 1,
                            parallel = FALSE,
                            weights = 0,
                            smoothness_type = "Second_order",
                            sparse_tuning_type = "SCAD",
                            tuning_order = "Sparsity",
                            cv.pick = "min",
                            sparse_tuning_u = NULL,
                            sparse_tuning_v = NULL,
                            smooth_tuning_u = NULL,
                            smooth_tuning_v = NULL)


  # Results
  algo <- SimulTest_ss_uv


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_ss_uv <- u1_est; u2_est_ss_uv <- u2_est
  v11_est_ss_uv <- v11_est; v12_est_ss_uv <- v12_est
  v21_est_ss_uv <- v21_est; v22_est_ss_uv <- v22_est

  # MSE
  mse_v11_ss_uv <- mse_L2(v11, v11_est)
  mse_v12_ss_uv <- mse_L2(v12, v12_est)
  mse_v21_ss_uv <- mse_L2(v21, v21_est)
  mse_v22_ss_uv <- mse_L2(v22, v22_est)
  mse_V_ss_uv <- (mse_v11_ss_uv + mse_v12_ss_uv + mse_v21_ss_uv + mse_v22_ss_uv) / 4

  mse_u1_ss_uv <- mse_L2(u1, u1_est)
  mse_u2_ss_uv <- mse_L2(u2, u2_est)
  mse_U_ss_uv <- (mse_u1_ss_uv + mse_u1_ss_uv) / 2

  # ISE
  ISE_u1_ss_uv <- ise(u1, u1_est)
  ISE_u2_ss_uv <- ise(u2, u2_est)
  ISE_v11_ss_uv <- ise(v11, v11_est)
  ISE_v12_ss_uv <- ise(v12, v12_est)
  ISE_v21_ss_uv <- ise(v21, v21_est)
  ISE_v22_ss_uv <- ise(v22, v22_est)


  ###### SVD ######
  X1obj <- fdClass(X1, Smoothing_parameter = 0,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
  X2obj <- fdClass(X2, Smoothing_parameter = 0,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                  Sparsity_parameter = 0)
  w1 <- 1 / mean( apply(X1, 2, var) )
  w2 <- 1 / mean( apply(X2, 2, var) )

  print("SVD")
  SimulTest_SVD <- ReMPCA(hd = Xobj,
                          centerhds = FALSE,
                          num_pcs = 2,
                          nfolds_u = 5,
                          nfolds_v = NULL,
                          thresh = 1e-10,
                          maxit = 100,
                          tuning_iter = 1,
                          parallel = FALSE,
                          weights = 0,
                          smoothness_type = "Second_order",
                          sparse_tuning_type = "SCAD",
                          tuning_order = "Sparsity",
                          cv.pick = "min",
                          sparse_tuning_u = NULL,
                          sparse_tuning_v = NULL,
                          smooth_tuning_u = NULL,
                          smooth_tuning_v = NULL)

  # Results
  algo <- SimulTest_SVD


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_SVD <- u1_est; u2_est_SVD  <- u2_est
  v11_est_SVD  <- v11_est; v12_est_SVD  <- v12_est
  v21_est_SVD  <- v21_est; v22_est_SVD <- v22_est

  # MSE
  mse_v11_SVD  <- mse_L2(v11, v11_est)
  mse_v12_SVD  <- mse_L2(v12, v12_est)
  mse_v21_SVD  <- mse_L2(v21, v21_est)
  mse_v22_SVD  <- mse_L2(v22, v22_est)
  mse_V_SVD  <- (mse_v11_SVD  + mse_v12_SVD  + mse_v21_SVD  + mse_v22_SVD ) / 4

  mse_u1_SVD  <- mse_L2(u1, u1_est)
  mse_u2_SVD  <- mse_L2(u2, u2_est)
  mse_U_SVD  <- (mse_u1_SVD  + mse_u2_SVD) / 2

  # ISE
  ISE_u1_SVD  <- ise(u1, u1_est)
  ISE_u2_SVD  <- ise(u2, u2_est)
  ISE_v11_SVD  <- ise(v11, v11_est)
  ISE_v12_SVD <- ise(v12, v12_est)
  ISE_v21_SVD <- ise(v21, v21_est)
  ISE_v22_SVD <- ise(v22, v22_est)

  # Ratio Relative to two-way smooth + sparse
  R_ISE_u1_SVD <- ISE_u1_SVD/ISE_u1_ss_uv
  R_ISE_u2_SVD <- ISE_u2_SVD/ISE_u2_ss_uv
  R_ISE_v11_SVD <- ISE_v11_SVD/ISE_v11_ss_uv
  R_ISE_v12_SVD <- ISE_v12_SVD/ISE_v12_ss_uv
  R_ISE_v21_SVD <- ISE_v21_SVD/ISE_v21_ss_uv
  R_ISE_v22_SVD <- ISE_v22_SVD/ISE_v22_ss_uv


  ###### Sparsity and Smoothness on u only ######
  X1obj <- fdClass(X1, Smoothing_parameter = 0,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
  X2obj <- fdClass(X2, Smoothing_parameter = 0,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
                  Sparsity_parameter = c(0,15,30,35,40,42,45,50,55))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  print("Smoothness & Sparsity on u")
  SimulTest_ss_u <- ReMPCA(hd = Xobj,
                           centerhds = FALSE,
                           num_pcs = 2,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "SCAD",
                           tuning_order = "Sparsity",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = NULL)


  # Results
  algo <- SimulTest_ss_u


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_ss_u <- u1_est; u2_est_ss_u <- u2_est
  v11_est_ss_u <- v11_est; v12_est_ss_u <- v12_est
  v21_est_ss_u <- v21_est; v22_est_ss_u <- v22_est

  # MSE
  mse_v11_ss_u <- mse_L2(v11, v11_est)
  mse_v12_ss_u <- mse_L2(v12, v12_est)
  mse_v21_ss_u <- mse_L2(v21, v21_est)
  mse_v22_ss_u <- mse_L2(v22, v22_est)
  mse_V_ss_u <- (mse_v11_ss_u + mse_v11_ss_u + mse_v21_ss_u + mse_v22_ss_u) / 4

  mse_u1_ss_u <- mse_L2(u1, u1_est)
  mse_u2_ss_u <- mse_L2(u2, u2_est)
  mse_U_ss_u <- (mse_u1_ss_u + mse_u2_ss_u) / 2

  # ISE
  ISE_u1_ss_u <- ise(u1, u1_est)
  ISE_u2_ss_u <- ise(u2, u2_est)
  ISE_v11_ss_u <- ise(v11, v11_est)
  ISE_v12_ss_u <- ise(v12, v12_est)
  ISE_v21_ss_u <- ise(v21, v21_est)
  ISE_v22_ss_u <- ise(v22, v22_est)

  # Ratio Relative to two-way smooth + sparse
  R_ISE_u1_ss_u <- ISE_u1_ss_u/ISE_u1_ss_uv
  R_ISE_u2_ss_u <- ISE_u2_ss_u/ISE_u2_ss_uv
  R_ISE_v11_ss_u <- ISE_v11_ss_u/ISE_v11_ss_uv
  R_ISE_v12_ss_u <- ISE_v12_ss_u/ISE_v12_ss_uv
  R_ISE_v21_ss_u <- ISE_v21_ss_u/ISE_v21_ss_uv
  R_ISE_v22_ss_u <- ISE_v22_ss_u/ISE_v22_ss_uv



  ###### Two way Smoothness ######
  X1obj <- fdClass(X1, Smoothing_parameter = NULL,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
  X2obj <- fdClass(X2, Smoothing_parameter = NULL,
                   Sparsity_parameter = 0) #round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = NULL,
                  Sparsity_parameter = 0)

  print("Two way Smoothness")
  SimulTest_sm_uv <- ReMPCA(hd = Xobj,
                            centerhds = FALSE,
                            num_pcs = 2,
                            nfolds_u = 5,
                            nfolds_v = NULL,
                            thresh = 1e-10,
                            maxit = 100,
                            tuning_iter = 1,
                            parallel = FALSE,
                            weights = 0,
                            smoothness_type = "Second_order",
                            sparse_tuning_type = "SCAD",
                            tuning_order = "Sparsity",
                            cv.pick = "min",
                            sparse_tuning_u = NULL,
                            sparse_tuning_v = NULL,
                            smooth_tuning_u = NULL,
                            smooth_tuning_v = NULL)


  # Results
  algo <- SimulTest_sm_uv


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_sm_uv <- u1_est; u2_est_sm_uv <- u2_est
  v11_est_sm_uv <- v11_est; v12_est_sm_uv <- v12_est
  v21_est_sm_uv <- v21_est; v22_est_sm_uv <- v22_est

  # MSE
  mse_v11_sm_uv <- mse_L2(v11, v11_est)
  mse_v12_sm_uv <- mse_L2(v12, v12_est)
  mse_v21_sm_uv <- mse_L2(v21, v21_est)
  mse_v22_sm_uv <- mse_L2(v22, v22_est)
  mse_V_sm_uv <- (mse_v11_sm_uv + mse_v12_sm_uv + mse_v21_sm_uv + mse_v22_sm_uv) / 4

  mse_u1_sm_uv <- mse_L2(u1, u1_est)
  mse_u2_sm_uv <- mse_L2(u2, u2_est)
  mse_U_sm_uv <- (mse_u1_sm_uv + mse_u2_sm_uv) / 2

  # ISE
  ISE_u1_sm_uv <- ise(u1, u1_est)
  ISE_u2_sm_uv <- ise(u2, u2_est)
  ISE_v11_sm_uv <- ise(v11, v11_est)
  ISE_v12_sm_uv <- ise(v12, v12_est)
  ISE_v21_sm_uv <- ise(v21, v21_est)
  ISE_v22_sm_uv <- ise(v22, v22_est)


  # Ratio Relative to two-way smooth + sparse
  R_ISE_u1_sm_uv <- ISE_u1_sm_uv/ISE_u1_ss_uv
  R_ISE_u2_sm_uv <- ISE_u2_sm_uv/ISE_u2_ss_uv
  R_ISE_v11_sm_uv <- ISE_v11_sm_uv/ISE_v11_ss_uv
  R_ISE_v12_sm_uv <- ISE_v12_sm_uv/ISE_v12_ss_uv
  R_ISE_v21_sm_uv <- ISE_v21_sm_uv/ISE_v21_ss_uv
  R_ISE_v22_sm_uv <- ISE_v22_sm_uv/ISE_v22_ss_uv


  ###### Two-way Sparsity ######
  X1obj <- fdClass(X1, Smoothing_parameter = 0,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
  X2obj <- fdClass(X2, Smoothing_parameter = 0,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                  Sparsity_parameter = c(0,15,30,35,40,42,45,50,55))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  print("Two-way Sparsity")
  SimulTest_sp_uv <- ReMPCA(hd = Xobj,
                            centerhds = FALSE,
                            num_pcs = 2,
                            nfolds_u = 5,
                            nfolds_v = NULL,
                            thresh = 1e-10,
                            maxit = 100,
                            tuning_iter = 1,
                            parallel = FALSE,
                            weights = 0,
                            smoothness_type = "Second_order",
                            sparse_tuning_type = "SCAD",
                            tuning_order = "Sparsity",
                            cv.pick = "min",
                            sparse_tuning_u = NULL,
                            sparse_tuning_v = NULL,
                            smooth_tuning_u = NULL,
                            smooth_tuning_v = NULL)


  # Results
  algo <- SimulTest_sp_uv


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_sp_uv <- u1_est; u2_est_sp_uv <- u2_est
  v11_est_sp_uv <- v11_est; v12_est_sp_uv <- v12_est
  v21_est_sp_uv <- v21_est; v22_est_sp_uv <- v22_est

  # MSE
  mse_v11_sp_uv <- mse_L2(v11, v11_est)
  mse_v12_sp_uv <- mse_L2(v12, v12_est)
  mse_v21_sp_uv <- mse_L2(v21, v21_est)
  mse_v22_sp_uv <- mse_L2(v22, v22_est)
  mse_V_sp_uv <- (mse_v11_sp_uv + mse_v12_sp_uv + mse_v21_sp_uv + mse_v22_sp_uv) / 4

  mse_u1_sp_uv <- mse_L2(u1, u1_est)
  mse_u2_sp_uv <- mse_L2(u2, u2_est)
  mse_U_sp_uv <- (mse_u1_sp_uv + mse_u2_sp_uv) / 2

  # ISE
  ISE_u1_sp_uv <- ise(u1, u1_est)
  ISE_u2_sp_uv <- ise(u2, u2_est)
  ISE_v11_sp_uv <- ise(v11, v11_est)
  ISE_v12_sp_uv <- ise(v12, v12_est)
  ISE_v21_sp_uv <- ise(v21, v21_est)
  ISE_v22_sp_uv <- ise(v22, v22_est)

  # Ratio Relative to two-way smooth + sparse
  R_ISE_u1_sp_uv <- ISE_u1_sp_uv/ISE_u1_ss_uv
  R_ISE_u2_sp_uv <- ISE_u2_sp_uv/ISE_u2_ss_uv
  R_ISE_v11_sp_uv <- ISE_v11_sp_uv/ISE_v11_ss_uv
  R_ISE_v12_sp_uv <- ISE_v12_sp_uv/ISE_v12_ss_uv
  R_ISE_v21_sp_uv <- ISE_v21_sp_uv/ISE_v21_ss_uv
  R_ISE_v22_sp_uv <- ISE_v22_sp_uv/ISE_v22_ss_uv



  ###### Sparsity and Smoothness on v only ######
  X1obj <- fdClass(X1, Smoothing_parameter = NULL,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))
  X2obj <- fdClass(X2, Smoothing_parameter = NULL,
                   Sparsity_parameter = c(0,15,30,35,50,60,69))#round(seq(0,N-1, length.out = round(N/4, digits = 0)), digits = 0))

  Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                  Sparsity_parameter = 0)

  print("Smoothness & Sparsity on v")
  SimulTest_ss_v <- ReMPCA(hd = Xobj,
                           centerhds = FALSE,
                           num_pcs = 2,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "SCAD",
                           tuning_order = "Sparsity",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = NULL)


  # Results
  algo <- SimulTest_ss_v


  # Don't change!
  u2_est <- algo$PCScores[,2]
  u1_est <- algo$PCScores[,1]

  v12_est <- algo$PCFunctions[[2]][[1]]
  v12_est <- v12_est/norm_vec(v12_est)
  v11_est <- algo$PCFunctions[[1]][[1]]
  v11_est <- v11_est/norm_vec(v11_est)

  v22_est <- algo$PCFunctions[[2]][[2]]
  v22_est <- v22_est/norm_vec(v22_est)
  v21_est <- algo$PCFunctions[[1]][[2]]
  v21_est <- v21_est/norm_vec(v21_est)

  # Align signs
  if (mean((u1 - u1_est)^2) > mean((u1 + u1_est)^2)) u1_est <- -u1_est
  if (mean((u2 - u2_est)^2) > mean((u2 + u2_est)^2)) u2_est <- -u2_est
  if (mean((v11 - v11_est)^2) > mean((v11 + v11_est)^2)) v11_est <- -v11_est
  if (mean((v12 - v12_est)^2) > mean((v12 + v12_est)^2)) v12_est <- -v12_est
  if (mean((v21 - v21_est)^2) > mean((v21 + v21_est)^2)) v21_est <- -v21_est
  if (mean((v22 - v22_est)^2) > mean((v22 + v22_est)^2)) v22_est <- -v22_est

  # Change
  u1_est_ss_v <- u1_est; u2_est_ss_v <- u2_est
  v11_est_ss_v <- v11_est; v12_est_ss_v <- v12_est
  v21_est_ss_v <- v21_est; v22_est_ss_v <- v22_est

  # MSE
  mse_v11_ss_v <- mse_L2(v11, v11_est)
  mse_v12_ss_v <- mse_L2(v12, v12_est)
  mse_v21_ss_v <- mse_L2(v21, v21_est)
  mse_v22_ss_v <- mse_L2(v22, v22_est)
  mse_V_ss_v <- (mse_v11_ss_v + mse_v12_ss_v + mse_v21_ss_v + mse_v22_ss_v) / 4

  mse_u1_ss_v <- mse_L2(u1, u1_est)
  mse_u2_ss_v <- mse_L2(u2, u2_est)
  mse_U_ss_v <- (mse_u1_ss_v + mse_u2_ss_v) / 2

  # ISE
  ISE_u1_ss_v <- ise(u1, u1_est)
  ISE_u2_ss_v <- ise(u2, u2_est)
  ISE_v11_ss_v <- ise(v11, v11_est)
  ISE_v12_ss_v <- ise(v12, v12_est)
  ISE_v21_ss_v <- ise(v21, v21_est)
  ISE_v22_ss_v <- ise(v22, v22_est)

  # Ratio Relative to two-way smooth + sparse
  R_ISE_u1_ss_v <- ISE_u1_ss_v/ISE_u1_ss_uv
  R_ISE_u2_ss_v <- ISE_u2_ss_v/ISE_u2_ss_uv
  R_ISE_v11_ss_v <- ISE_v11_ss_v/ISE_v11_ss_uv
  R_ISE_v12_ss_v <- ISE_v12_ss_v/ISE_v12_ss_uv
  R_ISE_v21_ss_v <- ISE_v21_ss_v/ISE_v21_ss_uv
  R_ISE_v22_ss_v <- ISE_v22_ss_v/ISE_v22_ss_uv


  results.u1.table <- data.frame(param = rep("u1", 6),
                                 method = c("SVD",
                                            "Smooth & Sparse u",
                                            "Two-way Smoothness",
                                            "Two-way Sparsity",
                                            "Smooth & Sparse v",
                                            "Smooth & Sparse u & v"),
                                 MSE = c(mse_u1_SVD,
                                         mse_u1_ss_u,
                                         mse_u1_sm_uv,
                                         mse_u1_sp_uv,
                                         mse_u1_ss_v,
                                         mse_u1_ss_uv),
                                 ISE = c(ISE_u1_SVD,
                                         ISE_u1_ss_u,
                                         ISE_u1_sm_uv,
                                         ISE_u1_sp_uv,
                                         ISE_u1_ss_v,
                                         ISE_u1_ss_uv),
                                 R_ISE = c(R_ISE_u1_SVD,
                                           R_ISE_u1_ss_u,
                                           R_ISE_u1_sm_uv,
                                           R_ISE_u1_sp_uv,
                                           R_ISE_u1_ss_v,
                                           1))

  results.u2.table <- data.frame(param = rep("u2", 6),
                                 method = c("SVD",
                                            "Smooth & Sparse u",
                                            "Two-way Smoothness",
                                            "Two-way Sparsity",
                                            "Smooth & Sparse v",
                                            "Smooth & Sparse u & v"),
                                 MSE = c(mse_u2_SVD,
                                         mse_u2_ss_u,
                                         mse_u2_sm_uv,
                                         mse_u2_sp_uv,
                                         mse_u2_ss_v,
                                         mse_u2_ss_uv),
                                 ISE = c(ISE_u2_SVD,
                                         ISE_u2_ss_u,
                                         ISE_u2_sm_uv,
                                         ISE_u2_sp_uv,
                                         ISE_u2_ss_v,
                                         ISE_u2_ss_uv),
                                 R_ISE = c(R_ISE_u2_SVD,
                                           R_ISE_u2_ss_u,
                                           R_ISE_u2_sm_uv,
                                           R_ISE_u2_sp_uv,
                                           R_ISE_u2_ss_v,
                                           1))

  results.v11.table <- data.frame(param = rep("v11", 6),
                                  method = c("SVD",
                                             "Smooth & Sparse u",
                                             "Two-way Smoothness",
                                             "Two-way Sparsity",
                                             "Smooth & Sparse v",
                                             "Smooth & Sparse u & v"),
                                  MSE = c(mse_v11_SVD,
                                          mse_v11_ss_u,
                                          mse_v11_sm_uv,
                                          mse_v11_sp_uv,
                                          mse_v11_ss_v,
                                          mse_v11_ss_uv),
                                  ISE = c(ISE_v11_SVD,
                                          ISE_v11_ss_u,
                                          ISE_v11_sm_uv,
                                          ISE_v11_sp_uv,
                                          ISE_v11_ss_v,
                                          ISE_v11_ss_uv),
                                  R_ISE = c(R_ISE_v11_SVD,
                                            R_ISE_v11_ss_u,
                                            R_ISE_v11_sm_uv,
                                            R_ISE_v11_sp_uv,
                                            R_ISE_v11_ss_v,
                                            1))

  results.v12.table <- data.frame(param = rep("v12", 6),
                                  method = c("SVD",
                                             "Smooth & Sparse u",
                                             "Two-way Smoothness",
                                             "Two-way Sparsity",
                                             "Smooth & Sparse v",
                                             "Smooth & Sparse u & v"),
                                  MSE = c(mse_v12_SVD,
                                          mse_v12_ss_u,
                                          mse_v12_sm_uv,
                                          mse_v12_sp_uv,
                                          mse_v12_ss_v,
                                          mse_v12_ss_uv),
                                  ISE = c(ISE_v12_SVD,
                                          ISE_v12_ss_u,
                                          ISE_v12_sm_uv,
                                          ISE_v12_sp_uv,
                                          ISE_v12_ss_v,
                                          ISE_v12_ss_uv),
                                  R_ISE = c(R_ISE_v12_SVD,
                                            R_ISE_v12_ss_u,
                                            R_ISE_v12_sm_uv,
                                            R_ISE_v12_sp_uv,
                                            R_ISE_v12_ss_v,
                                            1))

  results.v21.table <- data.frame(param = rep("v21", 6),
                                  method = c("SVD",
                                             "Smooth & Sparse u",
                                             "Two-way Smoothness",
                                             "Two-way Sparsity",
                                             "Smooth & Sparse v",
                                             "Smooth & Sparse u & v"),
                                  MSE = c(mse_v21_SVD,
                                          mse_v21_ss_u,
                                          mse_v21_sm_uv,
                                          mse_v21_sp_uv,
                                          mse_v21_ss_v,
                                          mse_v21_ss_uv),
                                  ISE = c(ISE_v21_SVD,
                                          ISE_v21_ss_u,
                                          ISE_v21_sm_uv,
                                          ISE_v21_sp_uv,
                                          ISE_v21_ss_v,
                                          ISE_v21_ss_uv),
                                  R_ISE = c(R_ISE_v21_SVD,
                                            R_ISE_v21_ss_u,
                                            R_ISE_v21_sm_uv,
                                            R_ISE_v21_sp_uv,
                                            R_ISE_v21_ss_v,
                                            1))

  results.v22.table <- data.frame(param = rep("v22", 6),
                                  method = c("SVD",
                                             "Smooth & Sparse u",
                                             "Two-way Smoothness",
                                             "Two-way Sparsity",
                                             "Smooth & Sparse v",
                                             "Smooth & Sparse u & v"),
                                  MSE = c(mse_v22_SVD,
                                          mse_v22_ss_u,
                                          mse_v22_sm_uv,
                                          mse_v22_sp_uv,
                                          mse_v22_ss_v,
                                          mse_v22_ss_uv),
                                  ISE = c(ISE_v22_SVD,
                                          ISE_v22_ss_u,
                                          ISE_v22_sm_uv,
                                          ISE_v22_sp_uv,
                                          ISE_v22_ss_v,
                                          ISE_v22_ss_uv),
                                  R_ISE = c(R_ISE_v22_SVD,
                                            R_ISE_v22_ss_u,
                                            R_ISE_v22_sm_uv,
                                            R_ISE_v22_sp_uv,
                                            R_ISE_v22_ss_v,
                                            1))


  return(list(N = N, sigma = sigma, v11 = v11, v12 = v12,
              v21 = v21, v22 = v22, u1 = u1, u2 = u2,

              ResultsTabel = rbind(results.u1.table,
                                   results.u2.table,
                                   results.v11.table,
                                   results.v12.table,
                                   results.v21.table,
                                   results.v22.table),
              # SVD
              u1_est_SVD = u1_est_SVD,
              u2_est_SVD = u2_est_SVD,
              v11_est_SVD = v11_est_SVD,
              v12_est_SVD = v12_est_SVD,
              v21_est_SVD = v21_est_SVD,
              v22_est_SVD = v22_est_SVD,

              # Smoothness & Sparsity on u:
              u1_est_ss_u = u1_est_ss_u,
              u2_est_ss_u = u2_est_ss_u,
              v11_est_ss_u = v11_est_ss_u,
              v12_est_ss_u = v12_est_ss_u,
              v21_est_ss_u = v21_est_ss_u,
              v22_est_ss_u = v22_est_ss_u,

              # Two- Way Smoothness:
              u1_est_sm_uv = u1_est_sm_uv,
              u2_est_sm_uv = u2_est_sm_uv,
              v11_est_sm_uv = v11_est_sm_uv,
              v12_est_sm_uv = v12_est_sm_uv,
              v21_est_sm_uv = v21_est_sm_uv,
              v22_est_sm_uv = v22_est_sm_uv,

              # # Two- Way Sparsity:
              u1_est_sp_uv = u1_est_sp_uv,
              u2_est_sp_uv = u2_est_sp_uv,
              v11_est_sp_uv = v11_est_sp_uv,
              v12_est_sp_uv = v12_est_sp_uv,
              v21_est_sp_uv = v21_est_sp_uv,
              v22_est_sp_uv = v22_est_sp_uv,

              # Smoothness & Sparsity on v:
              u1_est_ss_v = u1_est_ss_v,
              u2_est_ss_v = u2_est_ss_v,
              v11_est_ss_v = v11_est_ss_v,
              v12_est_ss_v = v12_est_ss_v,
              v21_est_ss_v = v21_est_ss_v,
              v22_est_ss_v = v22_est_ss_v,

              # Two-way Sparsity and smoothness:
              u1_est_ss_uv = u1_est_ss_uv,
              u2_est_ss_uv = u2_est_ss_uv,
              v11_est_ss_uv = v11_est_ss_uv,
              v12_est_ss_uv = v12_est_ss_uv,
              v21_est_ss_uv = v21_est_ss_uv,
              v22_est_ss_uv = v22_est_ss_uv))
}

result1 <- TwoWaySimulation(N = 101, sigma = 0.05, random_seed =50)

