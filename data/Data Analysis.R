# library(ReMPCA)
#
# load("~/ReMPCA/data/electrical_power_data.rda")
# load("~/ReMPCA/data/motion_sense_data.rda")
#
# normfunc <- function(x) sqrt(sum(x^2))
#
# ######### Electrical Power Data #########
# raw1 <- as.matrix(electrical_power_data$power[ , -(1:3)])
# raw2 <- as.matrix(electrical_power_data$voltage[ , -(1:3)])
#
# # rescale to unit integrated variance
# w1 <- 1 / mean( apply(raw1, 2, var) )
# w2 <- 1 / mean( apply(raw2, 2, var) )
# X1_tilde <- sqrt(w1) * raw1
# X2_tilde <- sqrt(w2) * raw2
# X_tilde <- cbind(X1_tilde, X2_tilde)
#
# library(grDevices)   # for adjustcolor()
# cols <- adjustcolor("grey40", alpha.f = 0.05)
#
# X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# # # make a time‑of‑day grid (0 to 24 hours, 288 points)
# # timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
# # matplot(timeGrid, t(X1obj),
# #         type="l", lwd=1, col=cols,     # ← semi‑transparent lines
# #         xlab="Time (hours)", ylab="Active Power")
# #
# # matplot(timeGrid, t(X2obj),
# #         type="l", lwd=1, col=cols,     # ← semi‑transparent lines
# #         xlab="Time (hours)", ylab="Voltage")
#
#
# # SVD #
# Electrical_power_SVD <- ReMPCA(hd = Xobj,
#                                centerhds = FALSE,
#                                num_pcs = 3,
#                                nfolds_u = 5,
#                                nfolds_v = NULL,
#                                thresh = 1e-10,
#                                maxit = 100,
#                                tuning_iter = 1,
#                                parallel = FALSE,
#                                weights = 0,
#                                smoothness_type = "Second_order",
#                                sparse_tuning_type = "soft",
#                                tuning_order = "Sparsity",
#                                cv.pick = "min",
#                                sparse_tuning_u = NULL,
#                                sparse_tuning_v = NULL,
#                                smooth_tuning_u = NULL,
#                                smooth_tuning_v = NULL)
#
#
# Electrical_SVD_v11 <- Electrical_power_SVD$PCFunctions[[1]][[1]]
# Electrical_SVD_v12 <- Electrical_power_SVD$PCFunctions[[1]][[2]]
# Electrical_SVD_v21 <- Electrical_power_SVD$PCFunctions[[2]][[1]]
# Electrical_SVD_v22 <- Electrical_power_SVD$PCFunctions[[2]][[2]]
# Electrical_SVD_v31 <- Electrical_power_SVD$PCFunctions[[3]][[1]]
# Electrical_SVD_v32 <- Electrical_power_SVD$PCFunctions[[3]][[2]]
#
#
# # 2) pull out the corresponding score‐vectors and compute their singular‐values
# u1 <- Electrical_power_SVD$PCScores[,1]
# u2 <- Electrical_power_SVD$PCScores[,2]
# u3 <- Electrical_power_SVD$PCScores[,3]
#
# sigma1 <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v11,Electrical_SVD_v12))) )^2)  #sqrt( sum((X1_tilde %*% Electrical_SVD_v11)^2) )  # or sd(u1)*norm(u1) if u1 norm=1
# sigma2 <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v21,Electrical_SVD_v22))) )^2)  #sqrt( sum((X1_tilde %*% Electrical_SVD_v21)^2) )
# sigma3 <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2)  #sqrt( sum((X1_tilde %*% Electrical_SVD_v31)^2) )
#
# # 3) build “data‐scaled” PC‐curves for each variable
# pc1_pwr <- sigma1 * Electrical_SVD_v11
# pc2_pwr <- sigma2 * Electrical_SVD_v21
# pc3_pwr <- sigma3 * Electrical_SVD_v31
#
# pc1_vol <- sigma1 * Electrical_SVD_v12
# pc2_vol <- sigma2 * Electrical_SVD_v22
# pc3_vol <- sigma3 * Electrical_SVD_v32
#
# #  build “data‐scaled” PC‐curves for each variable
# pc1_pwr <- - sqrt( sum(as.matrix(raw1) %*% as.matrix(data.frame(c(Electrical_SVD_v11))) )^2) * Electrical_SVD_v11 # sigma*v
# pc2_pwr <- -sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v21,Electrical_SVD_v22))) )^2) * Electrical_SVD_v21 # sigma*v
# pc3_pwr <- -sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2) * Electrical_SVD_v31 # sigma*v
#
# pc1_vol <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2) * Electrical_SVD_v12 # sigma*v
# pc2_vol <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v21,Electrical_SVD_v22))) )^2) * Electrical_SVD_v22 # sigma*v
# pc3_vol <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2) * Electrical_SVD_v32 # sigma*v
#
# # 4) compute means
# mu_pwr <- colMeans(raw1)
# mu_vol <- colMeans(raw2)
#
# # 5) plot
# par(mfrow = c(2,3), mar = c(4,4,2,1))
#
# # ---- First row: Active Power ----
# # PC1 in red
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main='PC 1')
# #lines(timeGrid, pc1_pwr, col='black', lwd=2)
# #lines(timeGrid, +mu_pwr - pc1_pwr, col='seagreen4', lwd=2, type = "b", pch  = "+")
# lines(timeGrid, mu_pwr + pc1_pwr, col='black', lwd=2)#, type = "b", pch  = "+")
#
# # PC2 in blue
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main ='PC 2')
# #lines(timeGrid, pc2_pwr, col='black', lwd=2)
# lines(timeGrid, mu_pwr + pc2_pwr, col='black', lwd=2)#, type = "b", pch  = "+")
# #lines(timeGrid, (mu_pwr - pc2_pwr), col='firebrick', lty=2, type = "b", pch  = "+")
#
# # PC3 in black
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main = 'PC 3')
# #lines(timeGrid, -pc3_pwr, col='black', lwd=2)
# lines(timeGrid, mu_pwr + pc3_pwr, col='black', lwd=2)#, type = "b", pch  = "+")
# #lines(timeGrid, (mu_pwr - pc3_pwr), col='firebrick', lty=2, type = "b", pch  = "+")
#
# # ---- Second row: Voltage ----
# # PC1 in red
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# lines(timeGrid, mu_vol+pc1_vol, col='black', lwd=2)
# #lines(timeGrid, mu_vol + pc1_vol, col='seagreen4', lwd=2, type = "b", pch  = "+")
# #lines(timeGrid, (mu_vol - pc1_vol), col='firebrick', lty=2, type = "b", pch  = "+")
#
# # PC2 in blue
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# #lines(timeGrid, pc2_vol, col='black', lwd=2)
# #lines(timeGrid, mu_vol + pc2_vol, col='seagreen4', lwd=2, type = "b", pch  = "+")
# lines(timeGrid, (mu_vol - pc2_vol), col='black', lwd=2)#, type = "b", pch  = "+")
#
# # PC3 in black
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# #lines(timeGrid, pc3_vol, col='black', lwd=2)
# lines(timeGrid, mu_vol + pc3_vol, col='black', lwd=2)#, type = "b", pch  = "+")
# #lines(timeGrid, (mu_vol - pc3_vol), col='firebrick', lty=2, type = "b", pch  = "+")
#
#
#
#
#
#
#
#
#
# # Sparsity and Smoothness on v #
# X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# Electrical_power_SSV <- ReMPCA(hd = Xobj,
#                                centerhds = FALSE,
#                                num_pcs = 3,
#                                nfolds_u = 5,
#                                nfolds_v = NULL,
#                                thresh = 1e-10,
#                                maxit = 100,
#                                tuning_iter = 1,
#                                parallel = FALSE,
#                                weights = 0,
#                                smoothness_type = "Second_order",
#                                sparse_tuning_type = "soft",
#                                tuning_order = "Sparsity",
#                                cv.pick = "min",
#                                sparse_tuning_u = NULL,
#                                sparse_tuning_v = NULL,
#                                smooth_tuning_u = NULL,
#                                smooth_tuning_v = list(4.485274e-05, 4.485274e-05))
#
# Electrical_SSV_v11 <- Electrical_power_SSV$PCFunctions[[1]][[1]]
# Electrical_SSV_v12 <- Electrical_power_SSV$PCFunctions[[1]][[2]]
# Electrical_SSV_v21 <- Electrical_power_SSV$PCFunctions[[2]][[1]]
# Electrical_SSV_v22 <- Electrical_power_SSV$PCFunctions[[2]][[2]]
# Electrical_SSV_v31 <- Electrical_power_SSV$PCFunctions[[3]][[1]]
# Electrical_SSV_v32 <- Electrical_power_SSV$PCFunctions[[3]][[2]]
#
#
#
# par(mfrow=c(2,3))
# matplot(Electrical_SSV_v11, type='l', col='red',lwd=2, main = "v11")
# points(Electrical_SVD_v11, type = 'l',lwd=2)
# matplot(Electrical_SSV_v21, type='l', col='red',lwd=2, main = "v21")
# points(Electrical_SVD_v21, type = 'l',lwd=2)
# matplot(Electrical_SSV_v31, type='l', col='red',lwd=2, main = "v31")
# points(Electrical_SVD_v31, type = 'l',lwd=2)
# matplot(Electrical_SSV_v12, type='l', col='red',lwd=2, main = "v12")
# points(Electrical_SVD_v12, type = 'l',lwd=2)
# matplot(Electrical_SSV_v22, type='l', col='red',lwd=2, main = "v22")
# points(Electrical_SVD_v22, type = 'l',lwd=2)
# matplot(Electrical_SSV_v32, type='l', col='red',lwd=2, main = "v32")
# points(Electrical_SVD_v32, type = 'l',lwd=2)
#
#
# #  build “data‐scaled” PC‐curves for each variable
# pc1_pwr_SSV <- - sqrt( sum(as.matrix(raw1) %*% as.matrix(data.frame(c(Electrical_SSV_v11))) )^2) * Electrical_SSV_v11 # sigma*v
# pc2_pwr_SSV <- -sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v21,Electrical_SVD_v22))) )^2) * Electrical_SSV_v21 # sigma*v
# pc3_pwr_SSV <- -sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2) * Electrical_SSV_v31 # sigma*v
#
# pc1_vol_SSV <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2)  * Electrical_SSV_v12 # sigma*v
# pc2_vol_SSV <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v21,Electrical_SVD_v22))) )^2) * Electrical_SSV_v22 # sigma*v
# pc3_vol_SSV <- sqrt( sum(as.matrix(Xobj$matrix) %*% as.matrix(data.frame(c(Electrical_SVD_v31,Electrical_SVD_v32))) )^2) * Electrical_SSV_v32 # sigma*v
#
# # 4) compute means
# mu_pwr <- colMeans(raw1)
# mu_vol <- colMeans(raw2)
#
# # 5) plot
# par(mfrow = c(2,3), mar = c(4,4,2,1))
#
# # ---- First row: Active Power ----
# # PC1 in red
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main='PC 1')
# lines(timeGrid, mu_pwr + pc1_pwr_SSV, col='red', lwd=3)
# #lines(timeGrid, mu_pwr + pc1_pwr, col='black', lwd=1 , lty  = 2)
#
# # PC2 in blue
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main ='PC 2')
# lines(timeGrid, mu_pwr + pc2_pwr_SSV, col='red', lwd=3)
# lines(timeGrid, mu_pwr + pc2_pwr, col='black', lwd=1 , lty  = 2)
#
#
# # PC3 in black
# matplot(timeGrid, t(raw1), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Active Power', main = 'PC 3')
# lines(timeGrid, mu_pwr + pc3_pwr_SSV, col='red', lwd=3)
# lines(timeGrid, mu_pwr + pc3_pwr, col='black', lwd=1 , lty  = 2)
#
#
# # ---- Second row: Voltage ----
# # PC1 in red
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# lines(timeGrid, mu_vol + pc1_vol_SSV, col='red', lwd=3)
# lines(timeGrid, mu_vol+pc1_vol, col='black', lwd=1 , lty  = 2)
#
#
# # PC2 in blue
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# lines(timeGrid, mu_vol - pc2_vol_SSV, col='red', lwd=3)
# lines(timeGrid, (mu_vol - pc2_vol), col='black', lwd=1 , lty  = 2)
#
# # PC3 in black
# matplot(timeGrid, t(raw2), type='l', col=cols, lwd=1,
#         xlab='Time (hours)', ylab='Voltage')
# lines(timeGrid, mu_vol + pc3_vol_SSV, col='red', lwd=3)
# lines(timeGrid, mu_vol + pc3_vol, col='black', lwd=1 , lty  = 2)
#
#
#
#
#
#
#
#
#
#
#
# ######## Motion Sense Data #########
# X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
#
# w1 <- 1 / mean( apply(motion_sense_data$user_acceleration, 2, var) )
# w2 <- 1 / mean( apply(motion_sense_data$pitch_attitude, 2, var) )
#
# library(grDevices)   # for adjustcolor()
# cols <- adjustcolor("grey40", alpha.f = 0.5)
# timeGrid <- seq(0, 1, length.out = ncol(motion_sense_data$user_acceleration))
#
# matplot(
#   (as.matrix(motion_sense_data$user_acceleration)),
#   type = "l", col = cols, lwd = 1,
#   xlab = "Time",
#   ylab = "User acceleration")
#
# matplot(
#   (as.matrix(motion_sense_data$pitch_attitude)),
#   type = "l", col = cols, lwd = 1,
#   xlab = "Time",
#   ylab = "Pitch attitude")
#
#
# # SVD #
# Motion_Sense_SVD <- ReMPCA(hd = Xobj,
#                            centerhds = FALSE,
#                            num_pcs = 3,
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
# Motion_SVD_v11 <- Motion_Sense_SVD$PCFunctions[[1]][[1]]
# Motion_SVD_v12 <- Motion_Sense_SVD$PCFunctions[[1]][[2]]
# Motion_SVD_v21 <- Motion_Sense_SVD$PCFunctions[[2]][[1]]
# Motion_SVD_v22 <- Motion_Sense_SVD$PCFunctions[[2]][[2]]
# Motion_SVD_v31 <- Motion_Sense_SVD$PCFunctions[[3]][[1]]
# Motion_SVD_v32 <- Motion_Sense_SVD$PCFunctions[[3]][[2]]
#
#
# par(mfrow = c(2,3))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32')
#
# ##Smoothness on v #
# Motion_Sense_SmV <- ReMPCA(hd = Xobj,
#                            centerhds = FALSE,
#                            num_pcs = 3,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Smoothness",
#                            cv.pick = "min",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = list(4.485274e-05, 4.485274e-05))
#
# Motion_SmV_v11 <- Motion_Sense_SmV$PCFunctions[[1]][[1]]
# Motion_SmV_v12 <- Motion_Sense_SmV$PCFunctions[[1]][[2]]
# Motion_SmV_v21 <- Motion_Sense_SmV$PCFunctions[[2]][[1]]
# Motion_SmV_v22 <- Motion_Sense_SmV$PCFunctions[[2]][[2]]
# Motion_SmV_v31 <- Motion_Sense_SmV$PCFunctions[[3]][[1]]
# Motion_SmV_v32 <- Motion_Sense_SmV$PCFunctions[[3]][[2]]
#
#
# par(mfrow = c(2,3))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
# points(Motion_SmV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
# points(Motion_SmV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
# points(Motion_SmV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
# points(Motion_SmV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
# points(Motion_SmV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
# points(Motion_SmV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')
#
#
#
#
# # Sparsity and Smoothness on v #
# X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
# X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# Motion_Sense_SSV <- ReMPCA(hd = Xobj,
#                            centerhds = FALSE,
#                            num_pcs = 3,
#                            nfolds_u = 5,
#                            nfolds_v = NULL,
#                            thresh = 1e-10,
#                            maxit = 100,
#                            tuning_iter = 1,
#                            parallel = FALSE,
#                            weights = 0,
#                            smoothness_type = "Second_order",
#                            sparse_tuning_type = "soft",
#                            tuning_order = "Smoothness",
#                            cv.pick = "1se",
#                            sparse_tuning_u = NULL,
#                            sparse_tuning_v = NULL,
#                            smooth_tuning_u = NULL,
#                            smooth_tuning_v = list(4.485274e-05, 4.485274e-05))
#
# Motion_SSV_v11 <- Motion_Sense_SSV$PCFunctions[[1]][[1]]
# Motion_SSV_v12 <- Motion_Sense_SSV$PCFunctions[[1]][[2]]
# Motion_SSV_v21 <- Motion_Sense_SSV$PCFunctions[[2]][[1]]
# Motion_SSV_v22 <- Motion_Sense_SSV$PCFunctions[[2]][[2]]
# Motion_SSV_v31 <- Motion_Sense_SSV$PCFunctions[[3]][[1]]
# Motion_SSV_v32 <- Motion_Sense_SSV$PCFunctions[[3]][[2]]
#
# matplot(-Motion_SSV_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
# matplot(-Motion_SSV_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
# matplot(Motion_SSV_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
# matplot(-Motion_SSV_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
# matplot(-Motion_SSV_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
# matplot(Motion_SSV_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
#
# par(mfrow = c(2,3))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
# points(Motion_SSV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
# points(Motion_SSV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
# points(Motion_SSV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
# points(Motion_SSV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
# points(Motion_SSV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
# points(Motion_SSV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')
#
#
#
# # Sparsity and Smoothness on u and v #
# X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter =c(0,4,45,91,95)) #round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
# X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter =c(0,4,45,91,95)) #round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
#
# Xobj <- hdClass(list(X1obj,X2obj),
#                 Smoothing_parameter = 0,
#                 Sparsity_parameter = c(0,10,20,30,40,50,60,70,80,85))
#
# Motion_Sense_SS_UV <- ReMPCA(hd = Xobj,
#                              centerhds = FALSE,
#                              num_pcs = 3,
#                              nfolds_u = 5,
#                              nfolds_v = NULL,
#                              thresh = 1e-10,
#                              maxit = 100,
#                              tuning_iter = 1,
#                              parallel = FALSE,
#                              weights = 0,
#                              smoothness_type = "Second_order",
#                              sparse_tuning_type = "soft",
#                              tuning_order = "Smoothness",
#                              cv.pick = "1se",
#                              sparse_tuning_u = NULL,
#                              sparse_tuning_v = NULL,
#                              smooth_tuning_u = NULL,
#                              smooth_tuning_v = list(4.485274e-05,3.027727e-06))
#
#
#
# Motion_SS_UV_v11 <- Motion_Sense_SS_UV$PCFunctions[[1]][[1]]
# Motion_SS_UV_v12 <- Motion_Sense_SS_UV$PCFunctions[[1]][[2]]
# Motion_SS_UV_v21 <- Motion_Sense_SS_UV$PCFunctions[[2]][[1]]
# Motion_SS_UV_v22 <- Motion_Sense_SS_UV$PCFunctions[[2]][[2]]
# Motion_SS_UV_v31 <- Motion_Sense_SS_UV$PCFunctions[[3]][[1]]
# Motion_SS_UV_v32 <- Motion_Sense_SS_UV$PCFunctions[[3]][[2]]
#
# Motion_Sense_SS_UV0 <- ReMPCA(hd = Xobj,
#                              centerhds = FALSE,
#                              num_pcs = 3,
#                              nfolds_u = 5,
#                              nfolds_v = NULL,
#                              thresh = 1e-10,
#                              maxit = 100,
#                              tuning_iter = 1,
#                              parallel = FALSE,
#                              weights = 0,
#                              smoothness_type = "Second_order",
#                              sparse_tuning_type = "soft",
#                              tuning_order = "Smoothness",
#                              cv.pick = "1se",
#                              sparse_tuning_u = NULL,
#                              sparse_tuning_v = list(0,0),
#                              smooth_tuning_u = NULL,
#                              smooth_tuning_v = list(0,0))
#
# Motion_SS_UV_v110 <- Motion_Sense_SS_UV0$PCFunctions[[1]][[1]]
# Motion_SS_UV_v120 <- Motion_Sense_SS_UV0$PCFunctions[[1]][[2]]
# Motion_SS_UV_v210 <- Motion_Sense_SS_UV0$PCFunctions[[2]][[1]]
# Motion_SS_UV_v220 <- Motion_Sense_SS_UV0$PCFunctions[[2]][[2]]
# Motion_SS_UV_v310 <- Motion_Sense_SS_UV0$PCFunctions[[3]][[1]]
# Motion_SS_UV_v320 <- Motion_Sense_SS_UV0$PCFunctions[[3]][[2]]
#
# par(mfrow = c(2,3))
# matplot(Motion_SS_UV_v110, type = 'l', , col = 'gray', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
# points(Motion_SS_UV_v11, type = 'l', col = 'black', lwd=2,xlab='Time')
# matplot(Motion_SS_UV_v210, type = 'l', , col = 'gray', main = 'v21',xlab='Time',  lwd=2)
# points(Motion_SS_UV_v21, type = 'l', col = 'black', lwd=2,xlab='Time')
# matplot(Motion_SS_UV_v310, type = 'l', , col = 'gray', main = 'v31',xlab='Time',  lwd=2)
# points(Motion_SS_UV_v31, type = 'l', col = 'black', lwd=2,xlab='Time')
# matplot(Motion_SS_UV_v120, type = 'l', , col = 'gray', main = 'v12',xlab='Time', ylab='Pitch attitude', lwd=2)
# points(Motion_SS_UV_v12, type = 'l', col = 'black', lwd=2,xlab='Time')
# matplot(Motion_SS_UV_v220, type = 'l', , col = 'gray', main = 'v22',xlab='Time',  lwd=2)
# points(Motion_SS_UV_v22, type = 'l', col = 'black', lwd=2,xlab='Time')
#
#
#
# par(mfrow = c(2,3))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
# points(Motion_SS_UV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
# points(Motion_SS_UV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
# points(Motion_SS_UV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
# points(Motion_SS_UV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
# points(Motion_SS_UV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
# points(Motion_SS_UV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')
#
#
# plot(Motion_Sense_SVD$PCScores[,1], main = 'PC scores 1')
# plot(Motion_Sense_SVD$PCScores[,2], main = 'PC scores 2')
# plot(Motion_Sense_SVD$PCScores[,3], main = 'PC scores 3')
# plot(Motion_Sense_SS_UV$PCScores[,1])
# plot(Motion_Sense_SS_UV$PCScores[,2])
# abline(h=0,col='red')
# plot(Motion_Sense_SS_UV$PCScores[,3])
# abline(h=0,col='red')
#
