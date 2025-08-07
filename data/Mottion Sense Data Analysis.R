library(ReMPCA)

load("~/ReMPCA/data/motion_sense_data.rda")

normfunc <- function(x) sqrt(sum(x^2))

######## Motion Sense Data #########
X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)


w1 <- 1 / mean( apply(motion_sense_data$user_acceleration, 2, var) )
w2 <- 1 / mean( apply(motion_sense_data$pitch_attitude, 2, var) )

library(grDevices)   # for adjustcolor()
cols <- adjustcolor("grey40", alpha.f = 0.5)
timeGrid <- seq(0, 1, length.out = ncol(motion_sense_data$user_acceleration))

matplot(
  (as.matrix(motion_sense_data$user_acceleration)),
  type = "l", col = cols, lwd = 1,
  xlab = "Time",
  ylab = "User acceleration")

matplot(
  (as.matrix(motion_sense_data$pitch_attitude)),
  type = "l", col = cols, lwd = 1,
  xlab = "Time",
  ylab = "Pitch attitude")


# SVD #
Motion_Sense_SVD <- ReMPCA(hd = Xobj,
                           centerhds = TRUE,
                           num_pcs = 3,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = c(w1,w2),
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Smoothness",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = NULL)

Motion_SVD_v11 <- Motion_Sense_SVD$PCFunctions[[1]][[1]]
Motion_SVD_v12 <- Motion_Sense_SVD$PCFunctions[[1]][[2]]
Motion_SVD_v21 <- Motion_Sense_SVD$PCFunctions[[2]][[1]]
Motion_SVD_v22 <- Motion_Sense_SVD$PCFunctions[[2]][[2]]
Motion_SVD_v31 <- Motion_Sense_SVD$PCFunctions[[3]][[1]]
Motion_SVD_v32 <- Motion_Sense_SVD$PCFunctions[[3]][[2]]


par(mfrow = c(2,3))
matplot(Motion_SVD_v11, type = 'l', main = 'v11')
matplot(Motion_SVD_v12, type = 'l', main = 'v12')
matplot(Motion_SVD_v21, type = 'l', main = 'v21')
matplot(Motion_SVD_v22, type = 'l', main = 'v22')
matplot(Motion_SVD_v31, type = 'l', main = 'v31')
matplot(Motion_SVD_v32, type = 'l', main = 'v32')

##Smoothness on v #
Motion_Sense_SmV <- ReMPCA(hd = Xobj,
                           centerhds = TRUE,
                           num_pcs = 3,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = c(w1,w2),
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Smoothness",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = list(4.485274e-05, 4.485274e-05))

Motion_SmV_v11 <- Motion_Sense_SmV$PCFunctions[[1]][[1]]
Motion_SmV_v12 <- Motion_Sense_SmV$PCFunctions[[1]][[2]]
Motion_SmV_v21 <- Motion_Sense_SmV$PCFunctions[[2]][[1]]
Motion_SmV_v22 <- Motion_Sense_SmV$PCFunctions[[2]][[2]]
Motion_SmV_v31 <- Motion_Sense_SmV$PCFunctions[[3]][[1]]
Motion_SmV_v32 <- Motion_Sense_SmV$PCFunctions[[3]][[2]]


par(mfrow = c(2,3))
matplot(-Motion_SmV_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
matplot(-Motion_SmV_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
matplot(-Motion_SmV_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
matplot(-Motion_SmV_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
matplot(-Motion_SmV_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
matplot(-Motion_SmV_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')


par(mfrow = c(2,3))
matplot(-Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(-Motion_SmV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
points(-Motion_SmV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
points(-Motion_SmV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
points(-Motion_SmV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
points(-Motion_SmV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
points(-Motion_SmV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')


# Sparsity and Smoothness on v #
X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))

X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
                 Smoothing_parameter = NULL,
                 Sparsity_parameter = round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

Motion_Sense_SSV <- ReMPCA(hd = Xobj,
                           centerhds = TRUE,
                           num_pcs = 3,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = c(w1,w2),
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Smoothness",
                           cv.pick = "min",
                           sparse_tuning_u = NULL,
                           sparse_tuning_v = NULL,
                           smooth_tuning_u = NULL,
                           smooth_tuning_v = list(4.485274e-05, 4.485274e-05))

Motion_SSV_v11 <- Motion_Sense_SSV$PCFunctions[[1]][[1]]
Motion_SSV_v12 <- Motion_Sense_SSV$PCFunctions[[1]][[2]]
Motion_SSV_v21 <- Motion_Sense_SSV$PCFunctions[[2]][[1]]
Motion_SSV_v22 <- Motion_Sense_SSV$PCFunctions[[2]][[2]]
Motion_SSV_v31 <- Motion_Sense_SSV$PCFunctions[[3]][[1]]
Motion_SSV_v32 <- Motion_Sense_SSV$PCFunctions[[3]][[2]]

matplot(-Motion_SSV_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
matplot(-Motion_SSV_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
matplot(-Motion_SSV_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
matplot(-Motion_SSV_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
matplot(-Motion_SSV_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
matplot(-Motion_SSV_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')

par(mfrow = c(2,3))
matplot(-Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(-Motion_SSV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
points(-Motion_SSV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
points(-Motion_SSV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
points(-Motion_SSV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
points(-Motion_SSV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(-Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
points(-Motion_SSV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')




# Sparsity and Smoothness on u and v #
X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter =round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))

X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter =round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj),
                Smoothing_parameter = 0,
                Sparsity_parameter = c(0,10,20,30,40,50,60,70,80,85))

Motion_Sense_SS_UV <- ReMPCA(hd = Xobj,
                             centerhds = TRUE,
                             num_pcs = 3,
                             nfolds_u = 5,
                             nfolds_v = NULL,
                             thresh = 1e-10,
                             maxit = 100,
                             tuning_iter = 1,
                             parallel = FALSE,
                             weights = c(w1,w2),
                             smoothness_type = "Second_order",
                             sparse_tuning_type = "hard",
                             tuning_order = "Sparsity",
                             cv.pick = "min",
                             sparse_tuning_u = NULL,
                             sparse_tuning_v = NULL,
                             smooth_tuning_u = NULL,
                             smooth_tuning_v = list(4.485274e-05, 4.485274e-05))



Motion_SS_UV_v11 <- Motion_Sense_SS_UV$PCFunctions[[1]][[1]]
Motion_SS_UV_v12 <- Motion_Sense_SS_UV$PCFunctions[[1]][[2]]
Motion_SS_UV_v21 <- Motion_Sense_SS_UV$PCFunctions[[2]][[1]]
Motion_SS_UV_v22 <- Motion_Sense_SS_UV$PCFunctions[[2]][[2]]
Motion_SS_UV_v31 <- Motion_Sense_SS_UV$PCFunctions[[3]][[1]]
Motion_SS_UV_v32 <- Motion_Sense_SS_UV$PCFunctions[[3]][[2]]

pdf("Motion Sense PCS.pdf", width=10, height=5)
par(
  mfrow  = c(2, 3),
  # bottom, left, top, right margins in lines:
  #   - very small bottom on the top row (we'll draw only bottom-row axes)
  #   - just enough left for col‑1 y‑labels
  #   - tiny residues on col‑2/3
  mar    = c(3, 3.5, 1.2, 1),
  # outer margins: only bottom for the “Time” label
  oma    = c(2.2, 0, 0, 0),
  mgp    = c(1.8, 0.4, 0),   # pull axis titles and labels in
  tcl    = -0.2,            # shorter ticks
  xaxs   = "i",             # no 4% padding on X
  yaxs   = "i"              # no 4% padding on Y
)
matplot(-Motion_SVD_v11, type = 'l', main = 'PC1', ylab='User acceleration', lwd=1.5, col = 'gray',xaxt = "n")
points(-Motion_SS_UV_v11, type = 'l', col = 'black', lwd=2)
matplot(-Motion_SVD_v21, type = 'l', main = 'PC2', lwd=1.5,ylab = "", col = 'gray',xaxt = "n")
points(-Motion_SS_UV_v21, type = 'l', col = 'black', lwd=2)
matplot(-Motion_SVD_v31, type = 'l', main = 'PC3', lwd=1.5,ylab = "", col = 'gray',xaxt = "n")
points(-Motion_SS_UV_v31, type = 'l', col = 'black', lwd=2)
matplot(-Motion_SVD_v12, type = 'l', lwd=1.5,xlab='Time', ylab='Pitch attitude', col = 'gray')
points(-Motion_SS_UV_v12, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(-Motion_SVD_v22, type = 'l', lwd=1.5,xlab='Time', col = 'gray',ylab = "")
points(-Motion_SS_UV_v22, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(-Motion_SVD_v32, type = 'l', lwd=1.5,xlab='Time', col = 'gray',ylab = "")
points(-Motion_SS_UV_v32, type = 'l', col = 'black', lwd=2,xlab='Time')
dev.off()


plot(Motion_Sense_SVD$PCScores[,1], main = 'PC scores 1')
plot(Motion_Sense_SVD$PCScores[,2], main = 'PC scores 2')
plot(Motion_Sense_SVD$PCScores[,3], main = 'PC scores 3')
plot(Motion_Sense_SS_UV$PCScores[,1])
plot(Motion_Sense_SS_UV$PCScores[,2])
abline(h=0,col='red')
plot(Motion_Sense_SS_UV$PCScores[,3])
abline(h=0,col='red')

