# library(ReMPCA)
#
# load("~/ReMPCA/data/electrical_power_data.rda")
# load("~/ReMPCA/data/motion_sense_data.rda")
#
# ######### Electrical Power Data #########
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
# # SVD #
# Electrical_power_SVD <- ReMPCA(hd = Xobj,
#                                centerhds = TRUE,
#                                num_pcs = 6,
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
# Electrical_SVD_v41 <- Electrical_power_SVD$PCFunctions[[4]][[1]]
# Electrical_SVD_v42 <- Electrical_power_SVD$PCFunctions[[4]][[2]]
# Electrical_SVD_v51 <- Electrical_power_SVD$PCFunctions[[5]][[1]]
# Electrical_SVD_v52 <- Electrical_power_SVD$PCFunctions[[5]][[2]]
# Electrical_SVD_v61 <- Electrical_power_SVD$PCFunctions[[6]][[1]]
# Electrical_SVD_v62 <- Electrical_power_SVD$PCFunctions[[6]][[2]]
#
# par(mfrow=c(2,3))
# matplot(Electrical_SVD_v11, type = 'l', main = "v11 - SVD")
# matplot(Electrical_SVD_v21, type = 'l', main = "v21 - SVD")
# matplot(Electrical_SVD_v31, type = 'l', main = "v31 - SVD")
# matplot(Electrical_SVD_v12, type = 'l', main = "v12 - SVD")
# matplot(Electrical_SVD_v22, type = 'l', main = "v22 - SVD")
# matplot(Electrical_SVD_v32, type = 'l', main = "v32 - SVD")
#
# matplot(Electrical_SVD_v41, type = 'l', main = "v41 - SVD")
# matplot(Electrical_SVD_v51, type = 'l', main = "v51 - SVD")
# matplot(Electrical_SVD_v61, type = 'l', main = "v61 - SVD")
# matplot(Electrical_SVD_v42, type = 'l', main = "v42 - SVD")
# matplot(Electrical_SVD_v52, type = 'l', main = "v52 - SVD")
# matplot(Electrical_SVD_v62, type = 'l', main = "v62 - SVD")
#
# # Sparsity and Smoothness on v #
# X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,287, length.out = round(288/10, digits = 0)), digits = 0))
#
# X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,287, length.out = round(288/10, digits = 0)), digits = 0))
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# Electrical_power_SSV <- ReMPCA(hd = Xobj,
#                                centerhds = TRUE,
#                                num_pcs = 6,
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
# Electrical_SSV_v11 <- Electrical_power_SSV$PCFunctions[[1]][[1]]
# Electrical_SSV_v12 <- Electrical_power_SSV$PCFunctions[[1]][[2]]
# Electrical_SSV_v21 <- Electrical_power_SSV$PCFunctions[[2]][[1]]
# Electrical_SSV_v22 <- Electrical_power_SSV$PCFunctions[[2]][[2]]
# Electrical_SSV_v31 <- Electrical_power_SSV$PCFunctions[[3]][[1]]
# Electrical_SSV_v32 <- Electrical_power_SSV$PCFunctions[[3]][[2]]
# Electrical_SSV_v41 <- Electrical_power_SSV$PCFunctions[[4]][[1]]
# Electrical_SSV_v42 <- Electrical_power_SSV$PCFunctions[[4]][[2]]
# Electrical_SSV_v51 <- Electrical_power_SSV$PCFunctions[[5]][[1]]
# Electrical_SSV_v52 <- Electrical_power_SSV$PCFunctions[[5]][[2]]
# Electrical_SSV_v61 <- Electrical_power_SSV$PCFunctions[[6]][[1]]
# Electrical_SSV_v62 <- Electrical_power_SSV$PCFunctions[[6]][[2]]
#
#
# par(mfrow=c(2,3))
# matplot(Electrical_SSV_v11, type='l')
# points(Electrical_SVD_v11, type = 'l', main = "v11")
# matplot(Electrical_SSV_v21, type='l')
# points(Electrical_SVD_v21, type = 'l', main = "v21")
# matplot(Electrical_SSV_v31, type='l')
# points(Electrical_SVD_v31, type = 'l', main = "v31")
# matplot(Electrical_SSV_v12, type='l')
# points(Electrical_SVD_v12, type = 'l', main = "v12")
# matplot(Electrical_SSV_v22, type='l')
# points(Electrical_SVD_v22, type = 'l', main = "v22")
# matplot(Electrical_SSV_v32, type='l')
# points(Electrical_SVD_v32, type = 'l', main = "v32")
#
#
# matplot(Electrical_SSV_v41, type='l')
# points(Electrical_SVD_v41, type = 'l', main = "v41")
# matplot(Electrical_SSV_v51, type='l')
# points(Electrical_SVD_v51, type = 'l', main = "v51")
# matplot(Electrical_SSV_v61, type='l')
# points(Electrical_SVD_v61, type = 'l', main = "v61")
# matplot(Electrical_SSV_v42, type='l')
# points(Electrical_SVD_v42, type = 'l', main = "v42")
# matplot(Electrical_SSV_v52, type='l')
# points(Electrical_SVD_v52, type = 'l', main = "v52")
# matplot(Electrical_SSV_v62, type='l')
# points(Electrical_SVD_v62, type = 'l', main = "v62")
#
#
#
#
#
#
#
#
#
# ######### Motion Sense Data #########
# X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# # SVD #
# Motion_Sense_SVD <- ReMPCA(hd = Xobj,
#                            centerhds = FALSE,
#                            num_pcs = 4,
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
#
# Motion_SVD_v11 <- Motion_Sense_SVD$PCFunctions[[1]][[1]]
# Motion_SVD_v12 <- Motion_Sense_SVD$PCFunctions[[1]][[2]]
# Motion_SVD_v21 <- Motion_Sense_SVD$PCFunctions[[2]][[1]]
# Motion_SVD_v22 <- Motion_Sense_SVD$PCFunctions[[2]][[2]]
# Motion_SVD_v31 <- Motion_Sense_SVD$PCFunctions[[3]][[1]]
# Motion_SVD_v32 <- Motion_Sense_SVD$PCFunctions[[3]][[2]]
# Motion_SVD_v41 <- Motion_Sense_SVD$PCFunctions[[4]][[1]]
# Motion_SVD_v42 <- Motion_Sense_SVD$PCFunctions[[4]][[2]]
#
#
#
# # Sparsity and Smoothness on v #
# X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
# X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# Motion_Sense_SSV <- ReMPCA(hd = Xobj,
#                            centerhds = TRUE,
#                            num_pcs = 4,
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
# Motion_SSV_v11 <- Motion_Sense_SSV$PCFunctions[[1]][[1]]
# Motion_SSV_v12 <- Motion_Sense_SSV$PCFunctions[[1]][[2]]
# Motion_SSV_v21 <- Motion_Sense_SSV$PCFunctions[[2]][[1]]
# Motion_SSV_v22 <- Motion_Sense_SSV$PCFunctions[[2]][[2]]
# Motion_SSV_v31 <- Motion_Sense_SSV$PCFunctions[[3]][[1]]
# Motion_SSV_v32 <- Motion_Sense_SSV$PCFunctions[[3]][[2]]
# Motion_SSV_v41 <- Motion_Sense_SSV$PCFunctions[[4]][[1]]
# Motion_SSV_v42 <- Motion_Sense_SSV$PCFunctions[[4]][[2]]
#
#
# par(mfrow = c(2,4))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11')
# points(Motion_SSV_v11, type = 'l', col = 'red')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12')
# points(Motion_SSV_v12, type = 'l', col = 'red')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21')
# points(Motion_SSV_v21, type = 'l', col = 'red')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22')
# points(Motion_SSV_v22, type = 'l', col = 'red')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31')
# points(Motion_SSV_v31, type = 'l', col = 'red')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32')
# points(Motion_SSV_v32, type = 'l', col = 'red')
# matplot(Motion_SVD_v41, type = 'l', main = 'v41')
# points(Motion_SSV_v41, type = 'l', col = 'red')
# matplot(Motion_SVD_v42, type = 'l', main = 'v42')
# points(Motion_SSV_v42, type = 'l', col = 'red')
#
#
#
# # Sparsity and Smoothness on u and v #
# X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
# X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
#                  Smoothing_parameter = NULL,
#                  Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))
#
#
# Xobj <- hdClass(list(X1obj,X2obj),
#                 Smoothing_parameter = NULL,
#                 Sparsity_parameter = round(seq(0,199, length.out = round(200/4, digits = 0)), digits = 0))
#
# Motion_Sense_SS_UV <- ReMPCA(hd = Xobj,
#                              centerhds = TRUE,
#                              num_pcs = 4,
#                              nfolds_u = 5,
#                              nfolds_v = NULL,
#                              thresh = 1e-10,
#                              maxit = 100,
#                              tuning_iter = 1,
#                              parallel = FALSE,
#                              weights = 0,
#                              smoothness_type = "Second_order",
#                              sparse_tuning_type = "soft",
#                              tuning_order = "Sparsity",
#                              cv.pick = "min",
#                              sparse_tuning_u = NULL,
#                              sparse_tuning_v = NULL,
#                              smooth_tuning_u = NULL,
#                              smooth_tuning_v = NULL)
#
#
#
# Motion_SS_UV_v11 <- Motion_Sense_SS_UV$PCFunctions[[1]][[1]]
# Motion_SS_UV_v12 <- Motion_Sense_SS_UV$PCFunctions[[1]][[2]]
# Motion_SS_UV_v21 <- Motion_Sense_SS_UV$PCFunctions[[2]][[1]]
# Motion_SS_UV_v22 <- Motion_Sense_SS_UV$PCFunctions[[2]][[2]]
# Motion_SS_UV_v31 <- Motion_Sense_SS_UV$PCFunctions[[3]][[1]]
# Motion_SS_UV_v32 <- Motion_Sense_SS_UV$PCFunctions[[3]][[2]]
# Motion_SS_UV_v41 <- Motion_Sense_SS_UV$PCFunctions[[4]][[1]]
# Motion_SS_UV_v42 <- Motion_Sense_SS_UV$PCFunctions[[4]][[2]]
#
#
# par(mfrow = c(2,4))
# matplot(Motion_SVD_v11, type = 'l', main = 'v11')
# points(Motion_SS_UV_v11, type = 'l', col = 'red')
# matplot(Motion_SVD_v12, type = 'l', main = 'v12')
# points(Motion_SS_UV_v12, type = 'l', col = 'red')
# matplot(Motion_SVD_v21, type = 'l', main = 'v21')
# points(Motion_SS_UV_v21, type = 'l', col = 'red')
# matplot(Motion_SVD_v22, type = 'l', main = 'v22')
# points(Motion_SS_UV_v22, type = 'l', col = 'red')
# matplot(Motion_SVD_v31, type = 'l', main = 'v31')
# points(Motion_SS_UV_v31, type = 'l', col = 'red')
# matplot(Motion_SVD_v32, type = 'l', main = 'v32')
# points(Motion_SS_UV_v32, type = 'l', col = 'red')
# matplot(Motion_SVD_v41, type = 'l', main = 'v41')
# points(Motion_SS_UV_v41, type = 'l', col = 'red')
# matplot(Motion_SVD_v42, type = 'l', main = 'v42')
# points(Motion_SS_UV_v42, type = 'l', col = 'red')
#
#
#
