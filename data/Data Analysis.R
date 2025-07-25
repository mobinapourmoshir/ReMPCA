library(ReMPCA)

load("~/ReMPCA/data/electrical_power_data.rda")
load("~/ReMPCA/data/motion_sense_data.rda")

normfunc <- function(x) sqrt(sum(x^2))

######### Electrical Power Data #########
raw1 <- as.matrix(electrical_power_data$power[ , -(1:3)])
raw2 <- as.matrix(electrical_power_data$voltage[ , -(1:3)])

# rescale to unit integrated variance
w1 <- 1 / mean( apply(raw1, 2, var) )
w2 <- 1 / mean( apply(raw2, 2, var) )
X1_tilde <- sqrt(w1) * raw1
X2_tilde <- sqrt(w2) * raw2
X_tilde <- cbind(X1_tilde, X2_tilde)

library(grDevices)   # for adjustcolor()
cols <- adjustcolor("grey40", alpha.f = 0.05)

X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)

X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = 0)


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

# # make a time‑of‑day grid (0 to 24 hours, 288 points)
# timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
# matplot(timeGrid, t(X1obj),
#         type="l", lwd=1, col=cols,     # ← semi‑transparent lines
#         xlab="Time (hours)", ylab="Active Power")
#
# matplot(timeGrid, t(X2obj),
#         type="l", lwd=1, col=cols,     # ← semi‑transparent lines
#         xlab="Time (hours)", ylab="Voltage")

w1 <- 1 / mean( apply(raw1, 2, var) )
w2 <- 1 / mean( apply(raw2, 2, var) )

# SVD #
Electrical_power_SVD <- ReMPCA(hd = Xobj,
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
                               sparse_tuning_type = "SCAD",
                               tuning_order = "Sparsity",
                               cv.pick = "min",
                               sparse_tuning_u = NULL,
                               sparse_tuning_v = NULL,
                               smooth_tuning_u = NULL,
                               smooth_tuning_v = NULL)


Electrical_SVD_v11 <- Electrical_power_SVD$PCFunctions[[1]][[1]]
Electrical_SVD_v12 <- Electrical_power_SVD$PCFunctions[[1]][[2]]
Electrical_SVD_v21 <- -Electrical_power_SVD$PCFunctions[[2]][[1]]
Electrical_SVD_v22 <- -Electrical_power_SVD$PCFunctions[[2]][[2]]
Electrical_SVD_v31 <- Electrical_power_SVD$PCFunctions[[3]][[1]]
Electrical_SVD_v32 <- Electrical_power_SVD$PCFunctions[[3]][[2]]


# plot
par(mfrow = c(2,3), mar = c(4,4,2,1))
timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))

matplot(timeGrid, Electrical_SVD_v11, type ='l', xlab = "Active Power")
matplot(timeGrid, Electrical_SVD_v21, type ='l')
matplot(timeGrid, Electrical_SVD_v31, type ='l')
matplot(timeGrid, Electrical_SVD_v12, type ='l', xlab ="Voltage")
matplot(timeGrid, Electrical_SVD_v22, type ='l')
matplot(timeGrid, Electrical_SVD_v32, type ='l')




# Sparsity and Smoothness on v #
X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = NULL)

X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
                 Smoothing_parameter = 0,
                 Sparsity_parameter = NULL)


Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
                Sparsity_parameter = 0)

Electrical_power_SSV <- ReMPCA(hd = Xobj,
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
                               sparse_tuning_type = "SCAD",
                               tuning_order = "Sparsity",
                               cv.pick = "min",
                               sparse_tuning_u = NULL,
                               sparse_tuning_v = NULL,
                               smooth_tuning_u = NULL,
                               smooth_tuning_v = list(4.485274e-05,3.027727e-06))

Electrical_SSV_v11 <- Electrical_power_SSV$PCFunctions[[1]][[1]]
Electrical_SSV_v12 <- Electrical_power_SSV$PCFunctions[[1]][[2]]
Electrical_SSV_v21 <- Electrical_power_SSV$PCFunctions[[2]][[1]]
Electrical_SSV_v22 <- Electrical_power_SSV$PCFunctions[[2]][[2]]
Electrical_SSV_v31 <- Electrical_power_SSV$PCFunctions[[3]][[1]]
Electrical_SSV_v32 <- Electrical_power_SSV$PCFunctions[[3]][[2]]


# Plot PCs
pdf("Electrical_PCS.pdf", width = 9, height = 6)

par(
  mfrow  = c(2, 3),
  # bottom, left, top, right margins in lines:
  #   - very small bottom on the top row (we'll draw only bottom-row axes)
  #   - just enough left for col‑1 y‑labels
  #   - tiny residues on col‑2/3
  mar    = c(2, 3.5, 1.2, 1),
  # outer margins: only bottom for the “Time” label
  oma    = c(2.2, 0, 0, 0),
  mgp    = c(1.8, 0.4, 0),   # pull axis titles and labels in
  tcl    = -0.2,            # shorter ticks
  xaxs   = "i",             # no 4% padding on X
  yaxs   = "i"              # no 4% padding on Y
)

# build a 0–24 time axis to match your data length
G        <- length(Electrical_SVD_v11)
time_seq <- seq(0, 24, length.out = G)
pct      <- c(34.1, 19.25, 7.19)

## Top row (no x‑axes)
plot(time_seq, Electrical_SVD_v11, type='l', col='red',
     main=bquote(PC~1~"(" * .(pct[1]) * "%)"),
     ylab="Active power", xaxt='n', xaxs="i", yaxs="i")
lines(time_seq, Electrical_SSV_v11)

plot(time_seq, Electrical_SVD_v21, type='l', col='red',
     main=bquote(PC~2~"(" * .(pct[2]) * "%)"),
     ylab="", xaxt='n')
lines(time_seq, -Electrical_SSV_v21)

plot(time_seq, Electrical_SVD_v31, type='l', col='red',
     main=bquote(PC~3~"(" * .(pct[3]) * "%)"),
     ylab="", xaxt='n')
lines(time_seq, Electrical_SSV_v31)

## Bottom row (draw x‑axes, only col‑1 y‑label)
plot(time_seq, Electrical_SVD_v12, type='l', col='red',
     ylab="Voltage", xaxt='n')
lines(time_seq, Electrical_SSV_v12)
axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))

plot(time_seq, Electrical_SVD_v22, type='l', col='red',
     ylab="", xaxt='n')
lines(time_seq, -Electrical_SSV_v22)
axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))

plot(time_seq, Electrical_SVD_v32, type='l', col='red',
     ylab="", xaxt='n')
lines(time_seq, Electrical_SSV_v32)
axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))

# single common x‑label
mtext("Time", side=1, outer=TRUE, line=1, cex=1.0)

dev.off()



#### Plot u
Electrical_SVD_u1 <- Electrical_power_SVD$PCScores[,1]
Electrical_SVD_u2 <- Electrical_power_SVD$PCScores[,2]
Electrical_SVD_u3 <- Electrical_power_SVD$PCScores[,3]

Electrical_SSV_u1 <- Electrical_power_SSV$PCScores[,1]
Electrical_SSV_u2 <- Electrical_power_SSV$PCScores[,2]
Electrical_SSV_u3 <- Electrical_power_SSV$PCScores[,3]

plot(data.frame(Electrical_SVD_u1,Electrical_SVD_u2,Electrical_SVD_u3), col= c('red', 'navyblue','darkgreen'))





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
                           centerhds = FALSE,
                           num_pcs = 3,
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
                           centerhds = FALSE,
                           num_pcs = 3,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
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
matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(Motion_SmV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
points(Motion_SmV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
points(Motion_SmV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
points(Motion_SmV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
points(Motion_SmV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
points(Motion_SmV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')




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
                           centerhds = FALSE,
                           num_pcs = 3,
                           nfolds_u = 5,
                           nfolds_v = NULL,
                           thresh = 1e-10,
                           maxit = 100,
                           tuning_iter = 1,
                           parallel = FALSE,
                           weights = 0,
                           smoothness_type = "Second_order",
                           sparse_tuning_type = "soft",
                           tuning_order = "Smoothness",
                           cv.pick = "1se",
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
matplot(Motion_SSV_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
matplot(-Motion_SSV_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
matplot(-Motion_SSV_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
matplot(Motion_SSV_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')

par(mfrow = c(2,3))
matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(Motion_SSV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
points(Motion_SSV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
points(Motion_SSV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
points(Motion_SSV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
points(Motion_SSV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
points(Motion_SSV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')



# Sparsity and Smoothness on u and v #
X1obj <- fdClass(t(as.matrix(motion_sense_data$user_acceleration)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter =c(0,4,45,91,95)) #round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))

X2obj <- fdClass(t(as.matrix(motion_sense_data$pitch_attitude)),
                 Smoothing_parameter = 0,
                 Sparsity_parameter =c(0,4,45,91,95)) #round(seq(0,95, length.out = round(96/4, digits = 0)), digits = 0))


Xobj <- hdClass(list(X1obj,X2obj),
                Smoothing_parameter = 0,
                Sparsity_parameter = c(0,10,20,30,40,50,60,70,80,85))

Motion_Sense_SS_UV <- ReMPCA(hd = Xobj,
                             centerhds = FALSE,
                             num_pcs = 3,
                             nfolds_u = 5,
                             nfolds_v = NULL,
                             thresh = 1e-10,
                             maxit = 100,
                             tuning_iter = 1,
                             parallel = FALSE,
                             weights = 0,
                             smoothness_type = "Second_order",
                             sparse_tuning_type = "soft",
                             tuning_order = "Smoothness",
                             cv.pick = "1se",
                             sparse_tuning_u = NULL,
                             sparse_tuning_v = NULL,
                             smooth_tuning_u = NULL,
                             smooth_tuning_v = list(4.485274e-05,3.027727e-06))



Motion_SS_UV_v11 <- Motion_Sense_SS_UV$PCFunctions[[1]][[1]]
Motion_SS_UV_v12 <- Motion_Sense_SS_UV$PCFunctions[[1]][[2]]
Motion_SS_UV_v21 <- Motion_Sense_SS_UV$PCFunctions[[2]][[1]]
Motion_SS_UV_v22 <- Motion_Sense_SS_UV$PCFunctions[[2]][[2]]
Motion_SS_UV_v31 <- Motion_Sense_SS_UV$PCFunctions[[3]][[1]]
Motion_SS_UV_v32 <- Motion_Sense_SS_UV$PCFunctions[[3]][[2]]

Motion_Sense_SS_UV0 <- ReMPCA(hd = Xobj,
                             centerhds = FALSE,
                             num_pcs = 3,
                             nfolds_u = 5,
                             nfolds_v = NULL,
                             thresh = 1e-10,
                             maxit = 100,
                             tuning_iter = 1,
                             parallel = FALSE,
                             weights = 0,
                             smoothness_type = "Second_order",
                             sparse_tuning_type = "soft",
                             tuning_order = "Smoothness",
                             cv.pick = "1se",
                             sparse_tuning_u = NULL,
                             sparse_tuning_v = list(0,0),
                             smooth_tuning_u = NULL,
                             smooth_tuning_v = list(0,0))

Motion_SS_UV_v110 <- Motion_Sense_SS_UV0$PCFunctions[[1]][[1]]
Motion_SS_UV_v120 <- Motion_Sense_SS_UV0$PCFunctions[[1]][[2]]
Motion_SS_UV_v210 <- Motion_Sense_SS_UV0$PCFunctions[[2]][[1]]
Motion_SS_UV_v220 <- Motion_Sense_SS_UV0$PCFunctions[[2]][[2]]
Motion_SS_UV_v310 <- Motion_Sense_SS_UV0$PCFunctions[[3]][[1]]
Motion_SS_UV_v320 <- Motion_Sense_SS_UV0$PCFunctions[[3]][[2]]

par(mfrow = c(2,3))
matplot(Motion_SS_UV_v110, type = 'l', , col = 'gray', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(Motion_SS_UV_v11, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(Motion_SS_UV_v210, type = 'l', , col = 'gray', main = 'v21',xlab='Time',  lwd=2)
points(Motion_SS_UV_v21, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(Motion_SS_UV_v310, type = 'l', , col = 'gray', main = 'v31',xlab='Time',  lwd=2)
points(Motion_SS_UV_v31, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(Motion_SS_UV_v120, type = 'l', , col = 'gray', main = 'v12',xlab='Time', ylab='Pitch attitude', lwd=2)
points(Motion_SS_UV_v12, type = 'l', col = 'black', lwd=2,xlab='Time')
matplot(Motion_SS_UV_v220, type = 'l', , col = 'gray', main = 'v22',xlab='Time',  lwd=2)
points(Motion_SS_UV_v22, type = 'l', col = 'black', lwd=2,xlab='Time')



par(mfrow = c(2,3))
matplot(Motion_SVD_v11, type = 'l', main = 'v11',xlab='Time', ylab='User acceleration', lwd=2)
points(Motion_SS_UV_v11, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v21, type = 'l', main = 'v21', lwd=2,xlab='Time')
points(Motion_SS_UV_v21, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v31, type = 'l', main = 'v31', lwd=2,xlab='Time')
points(Motion_SS_UV_v31, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v12, type = 'l', main = 'v12', lwd=2,xlab='Time', ylab='Pitch attitude')
points(Motion_SS_UV_v12, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v22, type = 'l', main = 'v22', lwd=2,xlab='Time')
points(Motion_SS_UV_v22, type = 'l', col = 'red', lwd=2,xlab='Time')
matplot(Motion_SVD_v32, type = 'l', main = 'v32', lwd=2,xlab='Time')
points(Motion_SS_UV_v32, type = 'l', col = 'red', lwd=2,xlab='Time')


plot(Motion_Sense_SVD$PCScores[,1], main = 'PC scores 1')
plot(Motion_Sense_SVD$PCScores[,2], main = 'PC scores 2')
plot(Motion_Sense_SVD$PCScores[,3], main = 'PC scores 3')
plot(Motion_Sense_SS_UV$PCScores[,1])
plot(Motion_Sense_SS_UV$PCScores[,2])
abline(h=0,col='red')
plot(Motion_Sense_SS_UV$PCScores[,3])
abline(h=0,col='red')

