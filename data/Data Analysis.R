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
# w1 <- 1 / mean( apply(raw1, 2, var) )
# w2 <- 1 / mean( apply(raw2, 2, var) )
#
# # SVD #
# Electrical_power_SVD <- ReMPCA(hd = Xobj,
#                                centerhds = TRUE,
#                                num_pcs = 3,
#                                nfolds_u = 5,
#                                nfolds_v = NULL,
#                                thresh = 1e-10,
#                                maxit = 100,
#                                tuning_iter = 1,
#                                parallel = FALSE,
#                                weights = c(w1,w2),
#                                smoothness_type = "Second_order",
#                                sparse_tuning_type = "SCAD",
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
# Electrical_SVD_v21 <- -Electrical_power_SVD$PCFunctions[[2]][[1]]
# Electrical_SVD_v22 <- -Electrical_power_SVD$PCFunctions[[2]][[2]]
# Electrical_SVD_v31 <- Electrical_power_SVD$PCFunctions[[3]][[1]]
# Electrical_SVD_v32 <- Electrical_power_SVD$PCFunctions[[3]][[2]]
#
#
# # plot
# par(mfrow = c(2,3), mar = c(4,4,2,1))
# timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
#
# matplot(timeGrid, Electrical_SVD_v11, type ='l', xlab = "Active Power")
# matplot(timeGrid, Electrical_SVD_v21, type ='l')
# matplot(timeGrid, Electrical_SVD_v31, type ='l')
# matplot(timeGrid, Electrical_SVD_v12, type ='l', xlab ="Voltage")
# matplot(timeGrid, Electrical_SVD_v22, type ='l')
# matplot(timeGrid, Electrical_SVD_v32, type ='l')
#
#
#
#
# # Sparsity and Smoothness on v #
# X1obj <- fdClass(as.matrix(electrical_power_data$power[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = NULL)
#
# X2obj <- fdClass(as.matrix(electrical_power_data$voltage[,-c(1,2,3)]),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = NULL)
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# Electrical_power_SSV <- ReMPCA(hd = Xobj,
#                                centerhds = TRUE,
#                                num_pcs = 3,
#                                nfolds_u = 5,
#                                nfolds_v = NULL,
#                                thresh = 1e-10,
#                                maxit = 100,
#                                tuning_iter = 1,
#                                parallel = FALSE,
#                                weights = c(w1,w2),
#                                smoothness_type = "Second_order",
#                                sparse_tuning_type = "SCAD",
#                                tuning_order = "Sparsity",
#                                cv.pick = "min",
#                                sparse_tuning_u = NULL,
#                                sparse_tuning_v = NULL,
#                                smooth_tuning_u = NULL,
#                                smooth_tuning_v = list(4.485274e-05,3.027727e-06))
#
# Electrical_SSV_v11 <- Electrical_power_SSV$PCFunctions[[1]][[1]]
# Electrical_SSV_v12 <- Electrical_power_SSV$PCFunctions[[1]][[2]]
# Electrical_SSV_v21 <- Electrical_power_SSV$PCFunctions[[2]][[1]]
# Electrical_SSV_v22 <- Electrical_power_SSV$PCFunctions[[2]][[2]]
# Electrical_SSV_v31 <- Electrical_power_SSV$PCFunctions[[3]][[1]]
# Electrical_SSV_v32 <- Electrical_power_SSV$PCFunctions[[3]][[2]]
#
#
# # Plot PCs
# pdf("Electrical_PCS.pdf", width = 9, height = 6)
#
# par(
#   mfrow  = c(2, 3),
#   # bottom, left, top, right margins in lines:
#   #   - very small bottom on the top row (we'll draw only bottom-row axes)
#   #   - just enough left for col‑1 y‑labels
#   #   - tiny residues on col‑2/3
#   mar    = c(2, 3.5, 1.2, 1),
#   # outer margins: only bottom for the “Time” label
#   oma    = c(2.2, 0, 0, 0),
#   mgp    = c(1.8, 0.4, 0),   # pull axis titles and labels in
#   tcl    = -0.2,            # shorter ticks
#   xaxs   = "i",             # no 4% padding on X
#   yaxs   = "i"              # no 4% padding on Y
# )
#
# # build a 0–24 time axis to match your data length
# G        <- length(Electrical_SVD_v11)
# time_seq <- seq(0, 24, length.out = G)
# pct      <- c(34.1, 19.25, 7.19)
#
# ## Top row (no x‑axes)
# plot(time_seq, Electrical_SVD_v11, type='l', col='red',
#      main=bquote(PC~1~"(" * .(pct[1]) * "%)"),
#      ylab="Active power", xaxt='n', xaxs="i", yaxs="i")
# lines(time_seq, Electrical_SSV_v11)
#
# plot(time_seq, Electrical_SVD_v21, type='l', col='red',
#      main=bquote(PC~2~"(" * .(pct[2]) * "%)"),
#      ylab="", xaxt='n')
# lines(time_seq, -Electrical_SSV_v21)
#
# plot(time_seq, Electrical_SVD_v31, type='l', col='red',
#      main=bquote(PC~3~"(" * .(pct[3]) * "%)"),
#      ylab="", xaxt='n')
# lines(time_seq, Electrical_SSV_v31)
#
# ## Bottom row (draw x‑axes, only col‑1 y‑label)
# plot(time_seq, Electrical_SVD_v12, type='l', col='red',
#      ylab="Voltage", xaxt='n')
# lines(time_seq, Electrical_SSV_v12)
# axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))
#
# plot(time_seq, Electrical_SVD_v22, type='l', col='red',
#      ylab="", xaxt='n')
# lines(time_seq, -Electrical_SSV_v22)
# axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))
#
# plot(time_seq, Electrical_SVD_v32, type='l', col='red',
#      ylab="", xaxt='n')
# lines(time_seq, Electrical_SSV_v32)
# axis(1, at=seq(0,24,by=6), mgp=c(1.8,0.4,0))
#
# # single common x‑label
# mtext("Time", side=1, outer=TRUE, line=1, cex=1.0)
#
# dev.off()
#
#
#
# #### Plot u
# Electrical_SVD_u1 <- Electrical_power_SVD$PCScores[,1]
# Electrical_SVD_u2 <- Electrical_power_SVD$PCScores[,2]
# Electrical_SVD_u3 <- Electrical_power_SVD$PCScores[,3]
#
# Electrical_SSV_u1 <- Electrical_power_SSV$PCScores[,1]
# Electrical_SSV_u2 <- Electrical_power_SSV$PCScores[,2]
# Electrical_SSV_u3 <- Electrical_power_SSV$PCScores[,3]
#
# plot(data.frame(Electrical_SVD_u1,Electrical_SVD_u2,Electrical_SVD_u3), col= c('red', 'navyblue','darkgreen'))
#
#
#
# # Holiday and Workday
# Powerdf <- electrical_power_data$power
# Voltagrdf <- electrical_power_data$voltage
# Powerdf$date <- as.Date(with(Powerdf, paste(year, month, day, sep = "-")), format = "%Y-%m-%d")
# Voltagrdf$date <- as.Date(with(Voltagrdf, paste(year, month, day, sep = "-")), format = "%Y-%m-%d")
#
# Powerdf$label <- ifelse(weekdays(Powerdf$date) %in% c("Saturday", "Sunday"), "Holiday", "Weekday")
# Voltagrdf$label <- ifelse(weekdays(Voltagrdf$date) %in% c("Saturday", "Sunday"), "Holiday", "Weekday")
# us_holidays <- as.Date(c("2006-12-25", "2007-01-01", "2007-07-04", "2007-12-25", "2008-01-01", "2008-12-25",
#                          "2009-01-01", "2009-07-04", "2009-12-25", "2010-01-01", "2010-07-04", "2010-12-25"))
#
# Powerdf$label[Powerdf$date %in% us_holidays] <- "Holiday"
# Voltagrdf$label[Voltagrdf$date %in% us_holidays] <- "Holiday"
#
# Powerdf_Holiday <- Powerdf[which(Powerdf$label=="Holiday"),]
# Powerdf_Holiday <-  Powerdf_Holiday[,-c(1,2,3,292,293)]
# Powerdf_Weekday <- Powerdf[which(Powerdf$label=="Weekday"),]
# Powerdf_Weekday <-  Powerdf_Weekday[,-c(1,2,3,292,293)]
#
# Voltagrdf_Holiday <- Voltagrdf[which(Voltagrdf$label=="Holiday"),]
# Voltagrdf_Holiday <-  Voltagrdf_Holiday[,-c(1,2,3,292,293)]
# Voltagrdf_Weekday <- Voltagrdf[which(Voltagrdf$label=="Weekday"),]
# Voltagrdf_Weekday <-  Voltagrdf_Weekday[,-c(1,2,3,292,293)]
#
#
# ##################  Holiday - Sparse and Smooth ##################
# X1obj <- fdClass(as.matrix(Powerdf_Holiday),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(as.matrix(Voltagrdf_Holiday),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# w1_Holiday <- 1 / mean( apply(as.matrix(Powerdf_Holiday), 2, var) )
# w2_Holiday <- 1 / mean( apply(as.matrix(Voltagrdf_Holiday), 2, var) )
#
# Electrical_Holiday_SSV <- ReMPCA(hd = Xobj,
#                                centerhds = FALSE,
#                                num_pcs = 1,
#                                nfolds_u = 5,
#                                nfolds_v = NULL,
#                                thresh = 1e-10,
#                                maxit = 100,
#                                tuning_iter = 1,
#                                parallel = FALSE,
#                                weights = c(w1_Holiday,w2_Holiday),
#                                smoothness_type = "Second_order",
#                                sparse_tuning_type = "SCAD",
#                                tuning_order = "Sparsity",
#                                cv.pick = "min",
#                                sparse_tuning_u = NULL,
#                                sparse_tuning_v = NULL,
#                                smooth_tuning_u = NULL,
#                                smooth_tuning_v = list(4.485274e-05,3.027727e-06))
#
# Electrical_SSV_v11_Holiday <- Electrical_Holiday_SSV$PCFunctions[[1]][[1]]
# Electrical_SSV_v12_Holiday <- Electrical_Holiday_SSV$PCFunctions[[1]][[2]]
# # Electrical_SSV_v21_Holiday <- Electrical_Holiday_SSV$PCFunctions[[2]][[1]]
# # Electrical_SSV_v22_Holiday <- Electrical_Holiday_SSV$PCFunctions[[2]][[2]]
# # Electrical_SSV_v31_Holiday <- Electrical_Holiday_SSV$PCFunctions[[3]][[1]]
# # Electrical_SSV_v32_Holiday <- Electrical_Holiday_SSV$PCFunctions[[3]][[2]]
#
# # plot
# #par(mfrow = c(2,3), mar = c(4,4,2,1))
# par(mfrow = c(1,2), mar = c(4,4,2,1))
# timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
# matplot(timeGrid, Electrical_SSV_v11_Holiday, type ='l', xlab = "Active Power")
# # matplot(timeGrid, Electrical_SSV_v21_Holiday, type ='l')
# # matplot(timeGrid, Electrical_SSV_v31_Holiday, type ='l')
# matplot(timeGrid, Electrical_SSV_v12_Holiday, type ='l', xlab ="Voltage")
# # matplot(timeGrid, Electrical_SSV_v22_Holiday, type ='l')
# # matplot(timeGrid, Electrical_SSV_v32_Holiday, type ='l')
#
#
#
# ##################  Workday - Sparse and Smooth ##################
# X1obj <- fdClass(as.matrix(Powerdf_Weekday),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
# X2obj <- fdClass(as.matrix(Voltagrdf_Weekday),
#                  Smoothing_parameter = 0,
#                  Sparsity_parameter = 0)
#
#
# Xobj <- hdClass(list(X1obj,X2obj), Smoothing_parameter = 0,
#                 Sparsity_parameter = 0)
#
# w1_Workday <- 1 / mean( apply(as.matrix(Powerdf_Weekday), 2, var) )
# w2_Workday <- 1 / mean( apply(as.matrix(Voltagrdf_Weekday), 2, var) )
#
# Electrical_Workday_SSV <- ReMPCA(hd = Xobj,
#                                  centerhds = FALSE,
#                                  num_pcs = 1,
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
#                                  smooth_tuning_v = list(4.485274e-05,3.027727e-06))
#
# Electrical_SSV_v11_Workday <- Electrical_Workday_SSV$PCFunctions[[1]][[1]]
# Electrical_SSV_v12_Workday <- Electrical_Workday_SSV$PCFunctions[[1]][[2]]
# # Electrical_SSV_v21_Workday <- Electrical_Workday_SSV$PCFunctions[[2]][[1]]
# # Electrical_SSV_v22_Workday <- Electrical_Workday_SSV$PCFunctions[[2]][[2]]
# # Electrical_SSV_v31_Workday <- Electrical_Workday_SSV$PCFunctions[[3]][[1]]
# # Electrical_SSV_v32_Workday <- Electrical_Workday_SSV$PCFunctions[[3]][[2]]
#
# # plot
# #par(mfrow = c(2,3), mar = c(4,4,2,1))
# par(mfrow = c(1,2), mar = c(4,4,2,1))
# timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
# matplot(timeGrid, -Electrical_SSV_v11_Workday, type ='l', xlab = "Active Power")
# # matplot(timeGrid, Electrical_SSV_v21_Holiday, type ='l')
# # matplot(timeGrid, Electrical_SSV_v31_Holiday, type ='l')
# matplot(timeGrid, -Electrical_SSV_v12_Workday, type ='l', xlab ="Voltage")
# # matplot(timeGrid, Electrical_SSV_v22_Holiday, type ='l')
# # matplot(timeGrid, Electrical_SSV_v32_Holiday, type ='l')
#
#
#
# ######################## Plots of Holiday/Workday ########################
# Electrical_SSV_v11_Workday <- Electrical_SSV_v11_Workday/ normfunc(Electrical_SSV_v11_Workday)
# Electrical_SSV_v11_Holiday <- Electrical_SSV_v11_Holiday/ normfunc(Electrical_SSV_v11_Holiday)
# Electrical_SSV_v12_Workday <- Electrical_SSV_v12_Workday/ normfunc(Electrical_SSV_v12_Workday)
# Electrical_SSV_v12_Holiday <- Electrical_SSV_v12_Holiday/ normfunc(Electrical_SSV_v12_Holiday)
#
# pdf("Workday-Holiday PCs.pdf", width=10, height=4)
# par(
#   mfrow = c(2,3),
#   mar   = c(1, 1, 2, 1),    # bottom, left, top, right – just enough for ticks + titles
#   oma   = c(4, 4, 0, 0),    # bottom, left, top, right – space for outer labels
#   mgp   = c(2, 0.5, 0),
#   tcl   = -0.3
# )
# timeGrid <- seq(0, 24, length.out = ncol(X1_tilde))
#
# # — compute common y–limits for each panel —
# ylim_power <- range(
#   -Electrical_SSV_v11_Workday,
#   Electrical_SSV_v11_Holiday,
#   na.rm = TRUE
# )
#
# ylim_voltage <- range(
#   -Electrical_SSV_v12_Workday,
#   Electrical_SSV_v12_Holiday,
#   na.rm = TRUE
# )
#
# ylim_power   <- grDevices::extendrange(ylim_power,   f = 0.05)
# ylim_voltage <- grDevices::extendrange(ylim_voltage, f = 0.05)
#
#
# matplot(timeGrid, -Electrical_SSV_v11_Workday, type ='l',col = "firebrick4",
#         ylab = "Active Power", ylim = ylim_power, main = "PCs",xaxt  = "n", xlab = "")
# points(timeGrid, Electrical_SSV_v11_Holiday, type ='l',col = "seagreen4")
#
# matplot(timeGrid, Electrical_SSV_v11_Holiday, type ='l',col = "seagreen4",
#         ylab ="", ylim = ylim_power, main = "Holiday",xaxt  = "n", xlab = "",yaxt  = "n")
# matplot(timeGrid, -Electrical_SSV_v11_Workday, type ='l',col = "firebrick4",
#         ylab ="", ylim = ylim_power, main = "Workday", xaxt  = "n",xlab = "",yaxt  = "n")
#
# matplot(timeGrid, -Electrical_SSV_v12_Workday, type ='l',col = "firebrick4",
#         ylab = "Voltage", ylim = ylim_voltage, xlab = "Time")
# points(timeGrid, Electrical_SSV_v12_Holiday, type ='l',col = "seagreen4", ylab ="")
#
# matplot(timeGrid, Electrical_SSV_v12_Holiday, type ='l',col = "seagreen4",
#         ylab ="", ylim = ylim_voltage, xlab = "Time",yaxt  = "n")
# matplot(timeGrid, -Electrical_SSV_v12_Workday, type ='l',col = "firebrick4",
#         ylab ="", ylim = ylim_voltage, xlab = "Time",yaxt  = "n")
#
# mtext("Time (hours)",        side=1, outer=TRUE, line=2.0, cex=1.1)
# mtext("Active Power", side=2, outer=TRUE, line=2, at=0.75, cex=1.1)
# mtext("Voltage",      side=2, outer=TRUE, line=2, at=0.25, cex=1.1)
#
# dev.off()
#
#
# ######################## Holiday/Workday PC scores ########################
# pdf("Workday-Holiday PC scores.pdf", width=7, height=5)
# u_Holiday <- Electrical_Holiday_SSV$PCScores[,1]
# u_Workday <- Electrical_Workday_SSV$PCScores[,1]
# par(mfrow=c(1,2), mar   = c(1, 2, 1, 1))
# boxplot(
#   u_Holiday,
#   main = "Holiday",
#   col  = "seagreen4"
# )
# boxplot(
#   u_Workday,
#   main = "Workday",
#   col  = "firebrick4"
# )
# dev.off()
#
#
