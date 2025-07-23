# pdf("MSE_boxplots.pdf", width = 10, height = 10)
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
#   oma   = c(3, 0, 0, 0),
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
# mtext("Simulation settings: N = 101, σ = 4; 100 replicates",
#       side   = 1,
#       line   = 1,
#       outer  = TRUE,
#       adj    = 0.5,
#       cex    = 0.8,
#       col    = "gray40")
#
# dev.off()
