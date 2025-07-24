# ###################################
# #           Box Plot              #
# ###################################
# pdf("MSE_boxplots.pdf", width = 13, height = 8)
# # 1) pull out your PC1/PC2 MSE vectors (as before)
# svd_pc1    <- sapply(combined_list, `[[`, "mse_svd")["PC1", ]
# smooth_pc1 <- sapply(combined_list, `[[`, "mse_smooth")["PC1", ]
# sparse_pc1 <- sapply(combined_list, `[[`, "mse_sparse")["PC1", ]
# ss_pc1     <- sapply(combined_list, `[[`, "mse_smooth_sparse")["PC1", ]
#
# svd_pc2    <- sapply(combined_list, `[[`, "mse_svd")["PC2", ]
# smooth_pc2 <- sapply(combined_list, `[[`, "mse_smooth")["PC2", ]
# sparse_pc2 <- sapply(combined_list, `[[`, "mse_sparse")["PC2", ]
# ss_pc2     <- sapply(combined_list, `[[`, "mse_smooth_sparse")["PC2", ]
#
# methods <- c("SVD","Smooth","Sparse","Smooth+\nSparse")
# cols    <- c("#1b9e77","darkgoldenrod1","#7570b3","brown1")
#
# # 2) common y‑limits so panels align
# ymin <- min(c(svd_pc1, smooth_pc1, sparse_pc1, ss_pc1,
#               svd_pc2, smooth_pc2, sparse_pc2, ss_pc2))
# ymax <- max(c(svd_pc1, smooth_pc1, sparse_pc1, ss_pc1,
#               svd_pc2, smooth_pc2, sparse_pc2, ss_pc2))
#
# # 3) Set up a 1×2 layout *and* reserve extra bottom space in the outer margins
# par(
#   mfrow = c(1,2),
#   mar   = c(5, 6, 4, 1),   # more space on the left
#   oma   = c(0, 0, 0, 0),
#   mgp   = c(4, 1, 0)       # axis titles one line further out
# )
#
# # PC1
# boxplot(svd_pc1, smooth_pc1, sparse_pc1, ss_pc1,
#         names = methods, col = cols, notch = TRUE,
#         outline = FALSE,
#         main = "PC 1: Distribution of MSE",
#         ylab = "MSE", las = 1)
# grid(nx=NA, ny=NULL, col="lightgray", lty="dotted")
#
# # PC2
# boxplot(svd_pc2, smooth_pc2, sparse_pc2, ss_pc2,
#         names = methods, col = cols, notch = TRUE,
#         outline = FALSE,
#         main = "PC 2: Distribution of MSE",
#         ylab = "", las = 1)
# grid(nx=NA, ny=NULL, col="lightgray", lty="dotted")
#
# # mtext("Simulation settings: N = 101, σ = 4; 100 replicates",
# #       side   = 1,
# #       line   = 1,
# #       outer  = TRUE,
# #       adj    = 0.5,
# #       cex    = 0.8,
# #       col    = "gray40")
#
# dev.off()
#
# ###################################
# #        Box Plot - MSE           #
# ###################################
# pdf("MSE_methods_boxplot.pdf", width = 5, height = 4.5)
#
# # 1×1 layout, small bottom margin + small outer bottom margin
# par(
#   mar = c(2, 6, 2, 1) + 0.1,  # bottom, left, top, right
#   oma = c(0, 0, 0, 0),        # just 1 line outer bottom
#   mgp = c(4, 1, 0)
# )
#
# methods <- c("SVD", "Smooth", "Sparse", "Smooth+\nSparse")
# cols    <- c("#1b9e77", "darkgoldenrod1", "#7570b3", "brown1")
#
# boxplot(
#   mse_svd,
#   mse_smooth,
#   mse_sparse,
#   mse_smooth_sparse,
#   names   = methods,
#   col     = cols,
#   notch   = TRUE,
#   outline = FALSE,
#   main    = "Distribution of MSE by Method",
#   ylab    = "MSE",
#   las     = 1,       # horizontal x-labels
#   cex.axis= 0.9      # slightly smaller axis text
# )
#
# grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
#
# # caption in the *outer* margin at line=0 (tight under axis)
# # mtext(
# #   expression(
# #     paste(
# #       "Simulation settings: N = 101, ",
# #       sigma, " = 4; 100 replicates"
# #     )
# #   ),
# #   side   = 1,
# #   line   = 0,
# #   outer  = TRUE,
# #   adj    = 0.5,
# #   cex    = 0.8,
# #   col    = "gray40"
# # )
#
# dev.off()
#
#
#
# ###################################
# #        Summary table            #
# ###################################
# library(kableExtra)
# mse_list <- list(
#   SVD            = mse_svd,
#   Smooth         = mse_smooth,
#   Sparse         = mse_sparse,
#   `Smooth+Sparse`= mse_smooth_sparse
# )
# methods <- names(mse_list)
#
# Q1     <- sapply(mse_list, function(x) quantile(x, .25))
# Median <- sapply(mse_list, median)
# Mean   <- sapply(mse_list, mean)
# Q3     <- sapply(mse_list, function(x) quantile(x, .75))
#
# summary_tbl <- data.frame(
#   Method = methods,
#   Q1     = Q1,
#   Mean   = Mean,
#   Q3     = Q3,
#   row.names = NULL,
#   check.names = FALSE)
#
# library(knitr)
# kable(
#   summary_tbl,
#   format    = "latex",
#   booktabs  = TRUE,
#   digits    = 7,
#   caption   = "Simulation summary of multivariate MSE by method",
#   row.names = FALSE
# )
#
#
#
# ###################################
# #    Summary table - PC 1,2       #
# ###################################
# library(dplyr)
# library(tidyr)
# library(kableExtra)
#
# # --- 1) build your long summaries as before ---
# pc1_df <- data.frame(
#   Method = rep(c("SVD","Smooth","Sparse","Smooth+Sparse"), each = length(svd_pc1)),
#   PC     = "FPC1",
#   Values = c(svd_pc1, smooth_pc1, sparse_pc1, ss_pc1)
# )
# pc2_df <- data.frame(
#   Method = rep(c("SVD","Smooth","Sparse","Smooth+Sparse"), each = length(svd_pc2)),
#   PC     = "FPC2",
#   Values = c(svd_pc2, smooth_pc2, sparse_pc2, ss_pc2)
# )
# df_all <- bind_rows(pc1_df, pc2_df)
#
# summaries <- df_all %>%
#   group_by(Method, PC) %>%
#   summarize(
#     Q1     = quantile(Values, .25),
#     Median = median(Values),
#     Mean   = mean(Values),
#     Q3     = quantile(Values, .75),
#     .groups = "drop"
#   )
#
# # --- 2) split into two tables ---
# table_pc1 <- summaries %>%
#   filter(PC == "FPC1") %>%
#   select(-PC)
#
# table_pc2 <- summaries %>%
#   filter(PC == "FPC2") %>%
#   select(-PC)
#
# # PC 1 table → LaTeX
# table_pc1 %>%
#   kbl(
#     format    = "latex",
#     booktabs  = TRUE,
#     digits    = 5,                        # show 7 significant digits
#     caption   = "PC 1: Quartiles and Mean by Method",
#     label     = "tab:pc1",
#     col.names = c("Method", "Q1", "Median", "Mean", "Q3")
#   ) %>%
#   kable_styling(latex_options = "hold_position")
#
#
#
# table_pc2 %>%
#   kbl(
#     format    = "latex",
#     booktabs  = TRUE,
#     digits    = 4,
#     caption   = "PC 2: Quartiles and Mean by Method",
#     label     = "tab:pc2",
#     col.names = c("Method", "Q1", "Median", "Mean", "Q3")
#   ) %>%
#   kable_styling(latex_options = c("hold_position"))
#
#
#
# ###################################
# #               PCs               #
# ###################################
# pdf("MFPCA-PC-Simulation1.pdf", width = 5, height = 6)
# # 1) compute a common y‑range
# allY <- c(result1$v11_est_SVD, result1$v12_est_SVD,
#           result1$v21_est_SVD, result1$v22_est_SVD)
# ylims <- range(allY)
#
# # 2) set up 2×2, small inner margins and space for outer labels/caption
# par(
#   mfrow = c(2,2),
#   mar   = c(0, 0, 1.5, 0),    # no x‐axis labels on top row, no y on right cols
#   oma   = c(7, 4, 2, 1)     # space for outer x, y, and caption
# )
#
# # Top‐left: Var 1, PC 1
# plot(result1$v11_est_SVD, type="l", col="red", ylim=ylims,
#      xaxt="n", xlab="", ylab="", yaxt="s", lwd = 1.5)
# lines(result1$v11, col="black", lty=2, lwd = 1.5)
# title("PC 1", line=0.5, cex.main=1)
# mtext("Var 1", side=2, line=2.5, outer=FALSE)
#
# # Top‐right: Var 1, PC 2
# plot(result1$v12_est_SVD, type="l", col="red", ylim=ylims,
#      xaxt="n", yaxt="n", xlab="", ylab="", lwd = 1.5)
# lines(result1$v12, col="black", lty=2, lwd = 1.5)
# title("PC 2", line=0.5, cex.main=1)
#
# # Bottom‐left: Var 2, PC 1
# plot(result1$v21_est_SVD, type="l", col="red", ylim=ylims,
#      yaxt="s", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v21, col="black", lty=2, lwd = 1.5)
# mtext("Var 2", side=2, line=2.5, outer=FALSE)
# mtext("Time",  side=1, line=2.5, outer=FALSE)
#
# # Bottom‐right: Var 2, PC 2
# plot(result1$v22_est_SVD, type="l", col="red", ylim=ylims,
#      yaxt="n", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v22, col="black", lty=2, lwd = 1.5)
# mtext("Time", side=1, line=2.5, outer=FALSE)
#
# # 3) Add the overall caption spanning the bottom
# mtext(
#   "(a) MFPCA",
#   side   = 1,
#   line   = 4,
#   outer  = TRUE,
#   cex    = 0.9
# )
# dev.off()
#
# #############
# pdf("SmoothMFPCA-PC-Simulation1.pdf", width = 5, height = 6)
# # 1) compute a common y‑range
# allY <- c(result1$v11_est_smooth_multi, result1$v12_est_smooth_multi,
#           result1$v21_est_smooth_multi, result1$v22_est_smooth_multi)
# ylims <- range(allY)
#
# # 2) set up 2×2, small inner margins and space for outer labels/caption
# par(
#   mfrow = c(2,2),
#   mar   = c(0, 0, 1.5, 0),    # no x‐axis labels on top row, no y on right cols
#   oma   = c(7, 4, 2, 1)     # space for outer x, y, and caption
# )
#
# # Top‐left: Var 1, PC 1
# plot(result1$v11_est_smooth_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", xlab="", ylab="", yaxt="s", lwd = 1.5)
# lines(result1$v11, col="black", lty=2, lwd = 1.5)
# title("PC 1", line=0.5, cex.main=1)
# mtext("Var 1", side=2, line=2.5, outer=FALSE)
#
# # Top‐right: Var 1, PC 2
# plot(result1$v12_est_smooth_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", yaxt="n", xlab="", ylab="", lwd = 1.5)
# lines(result1$v12, col="black", lty=2, lwd = 1.5)
# title("PC 2", line=0.5, cex.main=1)
#
# # Bottom‐left: Var 2, PC 1
# plot(result1$v21_est_smooth_multi, type="l", col="red", ylim=ylims,
#      yaxt="s", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v21, col="black", lty=2, lwd = 1.5)
# mtext("Var 2", side=2, line=2.5, outer=FALSE)
# mtext("Time",  side=1, line=2.5, outer=FALSE)
#
# # Bottom‐right: Var 2, PC 2
# plot(result1$v22_est_smooth_multi, type="l", col="red", ylim=ylims,
#      yaxt="n", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v22, col="black", lty=2, lwd = 1.5)
# mtext("Time", side=1, line=2.5, outer=FALSE)
#
# # 3) Add the overall caption spanning the bottom
# mtext(
#   "(b) Smoothed MFPCA",
#   side   = 1,
#   line   = 4,
#   outer  = TRUE,
#   cex    = 0.9
# )
# dev.off()
#
# #############
# pdf("SparseMFPCA-PC-Simulation1.pdf", width = 5, height = 6)
# # 1) compute a common y‑range
# allY <- c(result1$v11_est_sparse_multi, result1$v12_est_sparse_multi,
#           result1$v21_est_sparse_multi, result1$v22_est_sparse_multi)
# ylims <- range(allY)
#
# # 2) set up 2×2, small inner margins and space for outer labels/caption
# par(
#   mfrow = c(2,2),
#   mar   = c(0, 0, 1.5, 0),    # no x‐axis labels on top row, no y on right cols
#   oma   = c(7, 4, 2, 1)     # space for outer x, y, and caption
# )
#
# # Top‐left: Var 1, PC 1
# plot(result1$v11_est_sparse_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", xlab="", ylab="", yaxt="s", lwd = 1.5)
# lines(result1$v11, col="black", lty=2, lwd = 1.5)
# title("PC 1", line=0.5, cex.main=1)
# mtext("Var 1", side=2, line=2.5, outer=FALSE)
#
# # Top‐right: Var 1, PC 2
# plot(result1$v12_est_sparse_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", yaxt="n", xlab="", ylab="", lwd = 1.5)
# lines(result1$v12, col="black", lty=2, lwd = 1.5)
# title("PC 2", line=0.5, cex.main=1)
#
# # Bottom‐left: Var 2, PC 1
# plot(result1$v21_est_sparse_multi, type="l", col="red", ylim=ylims,
#      yaxt="s", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v21, col="black", lty=2, lwd = 1.5)
# mtext("Var 2", side=2, line=2.5, outer=FALSE)
# mtext("Time",  side=1, line=2.5, outer=FALSE)
#
# # Bottom‐right: Var 2, PC 2
# plot(result1$v22_est_sparse_multi, type="l", col="red", ylim=ylims,
#      yaxt="n", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v22, col="black", lty=2, lwd = 1.5)
# mtext("Time", side=1, line=2.5, outer=FALSE)
#
# # 3) Add the overall caption spanning the bottom
# mtext(
#   "(c) Sparse MFPCA",
#   side   = 1,
#   line   = 4,
#   outer  = TRUE,
#   cex    = 0.9
# )
# dev.off()
#
# #############
# pdf("SmSpMFPCA-PC-Simulation1.pdf", width = 5, height = 6)
# # 1) compute a common y‑range
# allY <- c(result1$v11_est_smooth_sparse_multi, result1$v12_est_smooth_sparse_multi,
#           result1$v21_est_smooth_sparse_multi, result1$v22_est_smooth_sparse_multi)
# ylims <- range(allY)
#
# # 2) set up 2×2, small inner margins and space for outer labels/caption
# par(
#   mfrow = c(2,2),
#   mar   = c(0, 0, 1.5, 0),    # no x‐axis labels on top row, no y on right cols
#   oma   = c(7, 4, 2, 1)     # space for outer x, y, and caption
# )
#
# # Top‐left: Var 1, PC 1
# plot(result1$v11_est_smooth_sparse_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", xlab="", ylab="", yaxt="s", lwd = 1.5)
# lines(result1$v11, col="black", lty=2, lwd = 1.5)
# title("PC 1", line=0.5, cex.main=1)
# mtext("Var 1", side=2, line=2.5, outer=FALSE)
#
# # Top‐right: Var 1, PC 2
# plot(result1$v12_est_smooth_sparse_multi, type="l", col="red", ylim=ylims,
#      xaxt="n", yaxt="n", xlab="", ylab="", lwd = 1.5)
# lines(result1$v12, col="black", lty=2, lwd = 1.5)
# title("PC 2", line=0.5, cex.main=1)
#
# # Bottom‐left: Var 2, PC 1
# plot(result1$v21_est_smooth_sparse_multi, type="l", col="red", ylim=ylims,
#      yaxt="s", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v21, col="black", lty=2, lwd = 1.5)
# mtext("Var 2", side=2, line=2.5, outer=FALSE)
# mtext("Time",  side=1, line=2.5, outer=FALSE)
#
# # Bottom‐right: Var 2, PC 2
# plot(result1$v22_est_smooth_sparse_multi, type="l", col="red", ylim=ylims,
#      yaxt="n", xlab="", ylab="", xaxt="s", lwd = 1.5)
# lines(result1$v22, col="black", lty=2, lwd = 1.5)
# mtext("Time", side=1, line=2.5, outer=FALSE)
#
# # 3) Add the overall caption spanning the bottom
# mtext(
#   "(d) Smoothed and Sparse MFPCA",
#   side   = 1,
#   line   = 4,
#   outer  = TRUE,
#   cex    = 0.9
# )
# dev.off()
#
#
#
#
