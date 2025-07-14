library(ReMPCA)

load("~/ReMPCA/data/electrical_power_data.rda")
load("~/ReMPCA/data/motion_sense_data.rda")

######### Electrical Power Data #########
X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

# SVD #
Electrical_power_SVD <- ReMPCA(hd = Xobj,
                               centerhds = TRUE,
                               num_pcs = 6,
                               nfolds_u = 5,
                               nfolds_v = NULL,
                               thresh = 1e-10,
                               maxit = 100,
                               tuning_iter = 1,
                               parallel = FALSE,
                               weights = 0,
                               smoothness_type = "Second_order",
                               sparse_tuning_type = "soft",
                               tuning_order = "Sparsity",
                               cv.pick = "min",
                               sparse_tuning_u = NULL,
                               sparse_tuning_v = NULL,
                               smooth_tuning_u = NULL,
                               smooth_tuning_v = NULL)



# Sparsity and Smoothness on v #
X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,287, length.out = round(288/4, digits = 0)), digits = 0))

X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,287, length.out = round(288/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

Electrical_power_SSV <- ReMPCA(hd = Xobj,
                               centerhds = TRUE,
                               num_pcs = 6,
                               nfolds_u = 5,
                               nfolds_v = NULL,
                               thresh = 1e-10,
                               maxit = 100,
                               tuning_iter = 1,
                               parallel = FALSE,
                               weights = 0,
                               smoothness_type = "Second_order",
                               sparse_tuning_type = "soft",
                               tuning_order = "Sparsity",
                               cv.pick = "min",
                               sparse_tuning_u = NULL,
                               sparse_tuning_v = NULL,
                               smooth_tuning_u = NULL,
                               smooth_tuning_v = NULL)


######### Motion Sense Data #########
X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

# SVD #
Motion_Sense_SVD <- ReMPCA(hd = Xobj,
                           centerhds = TRUE,
                           num_pcs = 4,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Sparsity",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = NULL)






# Sparsity and Smoothness on v #
X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))

X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

Motion_Sense_SSV <- ReMPCA(hd = Xobj,
                           centerhds = TRUE,
                           num_pcs = 4,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Sparsity",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = NULL)




# Sparsity and Smoothness on u and v #
X1obj <- fdClass(as.matrix(motion_sense_data$user_acceleration),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))

X2obj <- fdClass(as.matrix(motion_sense_data$pitch_attitude),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj),
                Smoothing_parameter = NULL,
                Sparsity_parameter = round(seq(0,199, length.out = round(200/4, digits = 0)), digits = 0))

Motion_Sense_SS_UV <- ReMPCA(hd = Xobj,
                             centerhds = TRUE,
                             num_pcs = 4,
                             nfolds_u = 5,
                             nfolds_v = NULL,
                             thresh = 1e-10,
                             maxit = 100,
                             tuning_iter = 1,
                             parallel = FALSE,
                             weights = 0,
                             smoothness_type = "Second_order",
                             sparse_tuning_type = "soft",
                             tuning_order = "Sparsity",
                             cv.pick = "min",
                             sparse_tuning_u = NULL,
                             sparse_tuning_v = NULL,
                             smooth_tuning_u = NULL,
                             smooth_tuning_v = NULL)







