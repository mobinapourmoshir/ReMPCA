############################################################
## two-way multivariate hybrid PCA simulation
## TRUE vs NOISY vs ReMPCA
############################################################

if (!requireNamespace("ReMPCA", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("mobinapourmoshir/ReMPCA")
}

library(ReMPCA)

set.seed(1001)

############################################################
## 1. Helpers
############################################################

unit <- function(x) {
  x <- as.numeric(x)
  x / sqrt(sum(x^2))
}

align_sign <- function(est, truth) {
  if (sum(est * truth, na.rm = TRUE) < 0) -est else est
}

rel_frob <- function(A, B) {
  norm(A - B, type = "F") / norm(B, type = "F")
}

split_blocks <- function(X, block_sizes) {
  ends <- cumsum(block_sizes)
  starts <- c(1, head(ends, -1) + 1)

  out <- vector("list", length(block_sizes))
  for (j in seq_along(block_sizes)) {
    out[[j]] <- X[, starts[j]:ends[j], drop = FALSE]
  }
  out
}

make_correlated_noise <- function(n, m, sigma = 0.05, rho_r = 0.25, rho_c = 0.10) {
  E <- matrix(rnorm(n * m), n, m)

  for (j in seq_len(m)) {
    E[, j] <- as.numeric(stats::filter(E[, j], rho_r, method = "recursive"))
  }

  for (i in seq_len(n)) {
    E[i, ] <- as.numeric(stats::filter(E[i, ], rho_c, method = "recursive"))
  }

  E[is.na(E)] <- rnorm(sum(is.na(E)))
  E <- E / sd(as.vector(E)) * sigma
  E
}

extract_v <- function(fit, j, pc = 1) {
  as.numeric(fit$PCFunctions[[pc]][[j]])
}

global_rank1_reconstruction <- function(fit, X_noisy, block_sizes, pc = 1) {
  uhat <- as.numeric(fit$PCScores[, pc])
  vhat_list <- lapply(seq_along(block_sizes), function(j) extract_v(fit, j, pc = pc))

  M_list <- lapply(vhat_list, function(v) outer(uhat, v))
  M_all <- do.call(cbind, M_list)
  X_all <- do.call(cbind, X_noisy)

  ## One global amplitude for the whole multivariate hybrid component.
  ## This removes arbitrary scaling between scores/functions.
  beta <- sum(X_all * M_all) / sum(M_all^2)

  R_all <- beta * M_all
  R_list <- split_blocks(R_all, block_sizes)

  list(
    blocks = R_list,
    beta = beta,
    uhat = uhat,
    vhat = vhat_list
  )
}

############################################################
## 2. Favorable but valid simulation design
############################################################

n <- 160
x <- seq(0, 2 * pi, length.out = n)

## Sparse + smooth row direction.
## Exactly matches the type of structure ReMPCA should recover.
u_true <- sin(x) + 0.30 * sin(2 * x)
u_true[1:70] <- 0
u_true <- unit(u_true)

## Four hybrid variables:
## X1: sparse regular
## X2: sparse functional / piecewise smooth
## X3: smooth functional
## X4: sparse regular block
m1 <- 24
m2 <- 50
m3 <- 42
m4 <- 18

t1 <- seq_len(m1)
t2 <- seq(0, 2 * pi, length.out = m2)
t3 <- seq(0, pi, length.out = m3)
t4 <- seq_len(m4)

v1_true <- rep(0, m1)
v1_true[c(2, 3, 4, 9, 10, 17, 22)] <- c(1.30, 0.10, 0.85, -1.25, -0.95, 1.00, -0.70)
v1_true <- unit(v1_true)

v2_true <- cos(2 * t2)
v2_true[26:m2] <- 0
v2_true <- unit(v2_true)

v3_true <- sin(t3)
v3_true <- unit(v3_true)

v4_true <- rep(0, m4)
v4_true[c(1, 2, 3, 10, 11, 12, 13, 17)] <- c(-1.20, 0.75, -0.95, 1.35, -0.60, -0.80)
v4_true <- unit(v4_true)

v_true <- list(v1_true, v2_true, v3_true, v4_true)

block_names <- c(
  "X1 regular sparse",
  "X2 functional piecewise",
  "X3 functional smooth",
  "X4 regular block-sparse"
)

block_sizes <- c(m1, m2, m3, m4)

## Block-specific signal sizes.
## This keeps all variables visible but does not let one block dominate everything.
signal <- c(3.2, 3.8, 3.5, 2.8)

X_true <- list(
  signal[1] * outer(u_true, v1_true),
  signal[2] * outer(u_true, v2_true),
  signal[3] * outer(u_true, v3_true),
  signal[4] * outer(u_true, v4_true)
)

names(X_true) <- block_names

## Moderate noise: visually noisy, but the common rank-one structure is dominant.
X_noisy <- list(
  X_true[[1]] + make_correlated_noise(n, m1, sigma = 0.025),
  X_true[[2]] + make_correlated_noise(n, m2, sigma = 0.055),
  X_true[[3]] + make_correlated_noise(n, m3, sigma = 0.045),
  X_true[[4]] + make_correlated_noise(n, m4, sigma = 0.025)
)

names(X_noisy) <- block_names

############################################################
## 3. ReMPCA objects
##
## Important:
## I use Smoothing_parameter = NULL here.
## Your README example also does this for functional blocks.
## This avoids possible argval/smoothing-scale instability.
############################################################

rd1 <- rdClass(
  data = as.matrix(X_noisy[[1]]),
  Sparsity_parameter = sort(unique(round(seq(0, m1 - 1, length.out = 16))))
)

fd2 <- fdClass(
  data = as.matrix(X_noisy[[2]]),
  argval = NULL,
  Smoothing_parameter = NULL,
  Sparsity_parameter = sort(unique(round(seq(0, m2 - 1, length.out = 22))))
)

fd3 <- fdClass(
  data = as.matrix(X_noisy[[3]]),
  argval = NULL,
  Smoothing_parameter = NULL,
  Sparsity_parameter = 0
)

rd4 <- rdClass(
  data = as.matrix(X_noisy[[4]]),
  Sparsity_parameter = sort(unique(round(seq(0, m4 - 1, length.out = 12))))
)

hd <- hdClass(
  hdlist = list(rd1, fd2, fd3, rd4),
  argval = NULL,
  Smoothing_parameter = NULL,
  Sparsity_parameter = sort(unique(round(seq(0, n - 1, length.out = 28))))
)

############################################################
## 4. Fit ReMPCA
############################################################

fit <- ReMPCA(
  hd = hd,
  centerhds = FALSE,
  num_pcs = 1,
  nfolds_u = 5,
  nfolds_v = NULL,
  thresh = 1e-10,
  maxit = 200,
  tuning_iter = 1,
  parallel = FALSE,
  weights = 0,
  smoothness_type = "Second_order",
  sparse_tuning_type = "soft",
  tuning_order = "Sparsity",
  cv.pick = "1se",
  sparse_tuning_u = NULL,
  sparse_tuning_v = NULL,
  smooth_tuning_u = NULL,
  smooth_tuning_v = NULL
)

############################################################
## 5. Extract estimates and reconstruct
############################################################

rec <- global_rank1_reconstruction(
  fit = fit,
  X_noisy = X_noisy,
  block_sizes = block_sizes,
  pc = 1
)

X_rempca <- rec$blocks
names(X_rempca) <- block_names

u_hat <- rec$uhat
v_hat <- rec$vhat

## Align signs for direction plots.
u_hat <- align_sign(unit(u_hat), u_true)

for (j in seq_along(v_hat)) {
  v_hat[[j]] <- align_sign(unit(v_hat[[j]]), v_true[[j]])
}

## Align reconstructed blocks for plotting/evaluation.
for (j in seq_along(X_rempca)) {
  if (sum(X_rempca[[j]] * X_true[[j]]) < 0) {
    X_rempca[[j]] <- -X_rempca[[j]]
  }
}

metrics <- data.frame(
  block = block_names,
  noisy_error = sapply(seq_along(X_true), function(j) rel_frob(X_noisy[[j]], X_true[[j]])),
  rempca_error = sapply(seq_along(X_true), function(j) rel_frob(X_rempca[[j]], X_true[[j]]))
)

metrics$improvement_percent <- 100 * (metrics$noisy_error - metrics$rempca_error) / metrics$noisy_error

print(metrics)

cat("\nPC1 variance explained reported by package:",
    round(100 * fit$VarianceExplained[1], 2), "%\n")
cat("Mean noisy relative error:",
    round(mean(metrics$noisy_error), 4), "\n")
cat("Mean ReMPCA relative error:",
    round(mean(metrics$rempca_error), 4), "\n")
cat("Mean improvement:",
    round(mean(metrics$improvement_percent), 1), "%\n")
cat("Global reconstruction amplitude:",
    round(rec$beta, 4), "\n")

############################################################
## 6. Paper plots
############################################################

out_dir <- "ReMPCA_valid_showcase"
dir.create(out_dir, showWarnings = FALSE)

pal <- grDevices::colorRampPalette(
  c("#053061", "#2166ac", "#f7f7f7", "#b2182b", "#67001f")
)(101)

save_pdf_png <- function(fun, name, width, height) {
  pdf(file.path(out_dir, paste0(name, ".pdf")),
      width = width, height = height, useDingbats = FALSE)
  fun()
  dev.off()

  png(file.path(out_dir, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 400)
  fun()
  dev.off()
}

img <- function(M, main = "", zlim = NULL) {
  if (is.null(zlim)) {
    z <- max(abs(M), na.rm = TRUE)
    zlim <- c(-z, z)
  }

  image(
    t(M),
    col = pal,
    zlim = zlim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main,
    useRaster = TRUE
  )
  box(lwd = 1.1)
}

############################################################
## Figure 1: main result
############################################################

fig_true_noisy_rempca <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))

  layout(matrix(1:12, nrow = 4, byrow = TRUE))

  par(
    mar = c(0.8, 0.8, 2.6, 0.8),
    oma = c(0.5, 0.5, 4.2, 0.5),
    family = "serif"
  )

  for (j in seq_along(X_true)) {
    z <- max(abs(c(X_true[[j]], X_noisy[[j]], X_rempca[[j]])), na.rm = TRUE)
    zlim <- c(-z, z)

    img(X_true[[j]], paste0(block_names[j], "\nTrue signal"), zlim)
    img(X_noisy[[j]], paste0(block_names[j], "\nNoisy data"), zlim)
    img(X_rempca[[j]], paste0(block_names[j], "\nReMPCA"), zlim)
  }

  mtext(
    "Two-way multivariate hybrid PCA: true signal, noisy data, and ReMPCA recovery",
    outer = TRUE,
    font = 2,
    cex = 1.25
  )
}

save_pdf_png(fig_true_noisy_rempca, "Figure_1_TRUE_NOISY_ReMPCA", 9.5, 10.5)

############################################################
## Figure 2: directions
############################################################

fig_directions <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))

  par(
    mfrow = c(2, 3),
    mar = c(3.2, 3.4, 2.5, 1.0),
    oma = c(0.2, 0.2, 3.5, 0.2),
    family = "serif"
  )

  plot(
    seq_along(u_true), u_true,
    type = "l", lwd = 2.8,
    ylim = range(c(u_true, u_hat)),
    xlab = "row index", ylab = "loading",
    main = "Shared row direction"
  )
  lines(seq_along(u_hat), u_hat, lwd = 2.8, lty = 2)
  abline(h = 0, lty = 3)
  legend("topright", c("True", "ReMPCA"), lwd = 2.8, lty = c(1, 2), bty = "n")

  plot(
    t1, v1_true,
    type = "n",
    ylim = range(c(v1_true, v_hat[[1]])),
    xlab = "variable index",
    ylab = "loading",
    main = "X1 regular sparse"
  )
  points(t1, v1_true, pch = 16, cex = 1.15)
  points(t1, v_hat[[1]], pch = 1, cex = 1.15)
  abline(h = 0, lty = 3)
  legend("topright", c("True", "ReMPCA"), pch = c(16, 1), bty = "n")

  plot(
    t2, v2_true,
    type = "l", lwd = 2.8,
    ylim = range(c(v2_true, v_hat[[2]])),
    xlab = "functional grid", ylab = "loading",
    main = "X2 functional piecewise"
  )
  lines(t2, v_hat[[2]], lwd = 2.8, lty = 2)
  abline(h = 0, lty = 3)

  plot(
    t3, v3_true,
    type = "l", lwd = 2.8,
    ylim = range(c(v3_true, v_hat[[3]])),
    xlab = "functional grid", ylab = "loading",
    main = "X3 functional smooth"
  )
  lines(t3, v_hat[[3]], lwd = 2.8, lty = 2)
  abline(h = 0, lty = 3)

  plot(
    t4, v4_true,
    type = "n",
    ylim = range(c(v4_true, v_hat[[4]])),
    xlab = "variable index",
    ylab = "loading",
    main = "X4 regular block-sparse"
  )
  points(t4, v4_true, pch = 16, cex = 1.15)
  points(t4, v_hat[[4]], pch = 1, cex = 1.15)
  abline(h = 0, lty = 3)
  legend("topright", c("True", "ReMPCA"), pch = c(16, 1), bty = "n")

  plot.new()
  text(
    0.02, 0.80,
    paste0(
      "Package PC1 variance explained: ",
      round(100 * fit$VarianceExplained[1], 1), "%\n",
      "Mean noisy relative error: ",
      round(mean(metrics$noisy_error), 3), "\n",
      "Mean ReMPCA relative error: ",
      round(mean(metrics$rempca_error), 3), "\n",
      "Mean improvement: ",
      round(mean(metrics$improvement_percent), 1), "%"
    ),
    adj = 0,
    cex = 1.10
  )

  mtext(
    "True and estimated two-way principal directions",
    outer = TRUE,
    font = 2,
    cex = 1.25
  )
}

save_pdf_png(fig_directions, "Figure_2_directions", 10.8, 7.0)

############################################################
## Figure 3: error summary
############################################################

fig_error <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))

  par(
    mar = c(5.2, 4.4, 3.3, 1.0),
    family = "serif"
  )

  Y <- rbind(metrics$noisy_error, metrics$rempca_error)
  colnames(Y) <- paste0("X", 1:4)

  bp <- barplot(
    Y,
    beside = TRUE,
    ylim = c(0, max(Y) * 1.35),
    ylab = "relative Frobenius error",
    main = "ReMPCA denoises the hybrid signal",
    legend.text = c("Noisy", "ReMPCA"),
    args.legend = list(bty = "n", x = "topright")
  )

  text(
    x = colMeans(bp),
    y = apply(Y, 2, max) * 1.12,
    labels = paste0(round(metrics$improvement_percent), "%"),
    cex = 1.0
  )

  mtext("percentage reduction in error", side = 3, line = 0.25, cex = 0.85)
}

save_pdf_png(fig_error, "Figure_3_error_summary", 7.0, 5.0)

cat("\nSaved figures in:\n")
cat(normalizePath(out_dir), "\n")
