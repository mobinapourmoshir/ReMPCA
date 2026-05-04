###############################################################################
# Regularized Multivariate Two-Way Hybrid PCA using ReMPCA
# Author: Mobina-style simulation workflow
#
# Goal:
#   1. Simulate multivariate hybrid data:
#        - functional variable 1
#        - functional variable 2
#        - regular sparse variable
#        - regular dense variable
#   2. Put them into ReMPCA objects:
#        fdClass(), rdClass(), hdClass()
#   3. Fit smooth + sparse two-way hybrid PCA:
#        ReMPCA()
#   4. Make nice plots:
#        - true vs estimated row score u
#        - true vs estimated block loadings v
#        - true/noisy/reconstructed heatmaps
#        - built-in ReMPCA plots
#   5. Save useful tables:
#        - block summary table
#        - tuning parameter table
#        - performance metrics table
###############################################################################

###############################################################################
# 0. Packages
###############################################################################

needed_packages <- c(
  "remotes",
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "knitr",
  "gt",
  "patchwork"
)

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

install_if_missing(needed_packages)

# Install your package from GitHub if it is not already installed
if (!requireNamespace("ReMPCA", quietly = TRUE)) {
  remotes::install_github("mobinapourmoshir/ReMPCA", force = TRUE)
}

library(ReMPCA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(knitr)
library(gt)
library(patchwork)

###############################################################################
# 1. Output folders
###############################################################################

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

###############################################################################
# 2. Helper functions
###############################################################################

unit_norm <- function(x) {
  x <- as.numeric(x)
  den <- sqrt(sum(x^2))
  if (den == 0) return(x)
  x / den
}

center_vec <- function(x) {
  as.numeric(x - mean(x))
}

rmse <- function(a, b) {
  sqrt(mean((as.numeric(a) - as.numeric(b))^2))
}

abs_cor <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
  abs(cor(a, b))
}

safe_extract <- function(object, name) {
  if (!is.null(object[[name]])) object[[name]] else NA
}

collapse_value <- function(x, digits = 4) {
  if (is.null(x)) return(NA_character_)
  x <- unlist(x)
  if (length(x) == 0) return(NA_character_)
  paste(signif(as.numeric(x), digits), collapse = ", ")
}

matrix_to_heat_df <- function(M, block, version) {
  as.data.frame(as.table(M)) |>
    as_tibble() |>
    transmute(
      row = as.integer(Var1),
      feature = as.integer(Var2),
      value = as.numeric(Freq),
      block = block,
      version = version
    )
}

safe_pkg_plot <- function(function_name, fit_object) {
  # Some plotting functions may be exported, and some may only exist in namespace.
  # This wrapper tries both safely.
  fn <- NULL

  if (exists(function_name, where = asNamespace("ReMPCA"), inherits = FALSE)) {
    fn <- get(function_name, envir = asNamespace("ReMPCA"))
  } else if (exists(function_name, mode = "function")) {
    fn <- get(function_name, mode = "function")
  }

  if (is.function(fn)) {
    try(fn(fit_object), silent = TRUE)
  } else {
    message("Plot function not found: ", function_name)
  }
}

###############################################################################
# 3. Simulate multivariate two-way hybrid data
###############################################################################

set.seed(2026)

# Number of observations / row domain points
n <- 150
s_grid <- seq(0, 1, length.out = n)

# True left singular vector u:
# smooth + sparse over the row direction.
u_true <- sin(2 * pi * s_grid) + 0.35 * cos(6 * pi * s_grid)
u_true[s_grid < 0.18] <- 0
u_true[s_grid > 0.82] <- 0
u_true <- center_vec(u_true)
u_true <- unit_norm(u_true)

# Block dimensions
m_fd1 <- 70
m_fd2 <- 55
q_rd1 <- 25
q_rd2 <- 14

t_fd1 <- seq(0, 1, length.out = m_fd1)
t_fd2 <- seq(0, 1, length.out = m_fd2)

# Functional loading 1: smooth and sparse on the domain
v_fd1_true <- sin(2 * pi * t_fd1) * exp(-2 * (t_fd1 - 0.45)^2)
v_fd1_true[t_fd1 < 0.12 | t_fd1 > 0.90] <- 0
v_fd1_true <- unit_norm(v_fd1_true)

# Functional loading 2: piecewise smooth
v_fd2_true <- cos(3 * pi * t_fd2)
v_fd2_true[t_fd2 > 0.55 & t_fd2 < 0.78] <- 0
v_fd2_true <- unit_norm(v_fd2_true)

# Regular loading 1: sparse scalar/vector covariates
v_rd1_true <- rep(0, q_rd1)
active_rd1 <- c(2, 5, 9, 13, 18, 24)
v_rd1_true[active_rd1] <- c(1.6, -1.2, 1.1, -0.9, 1.4, -1.0)
v_rd1_true <- unit_norm(v_rd1_true)

# Regular loading 2: dense regular covariates
v_rd2_true <- seq(-1, 1, length.out = q_rd2) + rnorm(q_rd2, sd = 0.08)
v_rd2_true <- unit_norm(v_rd2_true)

# Put true loadings into a list
v_true_list <- list(
  FD1_smooth_sparse = v_fd1_true,
  FD2_piecewise = v_fd2_true,
  RD1_sparse = v_rd1_true,
  RD2_dense = v_rd2_true
)

block_names <- names(v_true_list)
block_types <- c("Functional", "Functional", "Regular", "Regular")
block_dims <- c(m_fd1, m_fd2, q_rd1, q_rd2)

# Block-specific signal strengths and noise levels
signal_strength <- c(8.0, 7.0, 6.5, 5.5)
noise_sd <- c(0.25, 0.35, 0.20, 0.20)

# Generate true rank-one signal and noisy observed matrices
X_true_list <- map2(v_true_list, signal_strength, \(v, strength) {
  strength * outer(u_true, v)
})

X_noisy_list <- map2(X_true_list, noise_sd, \(X_signal, sd_noise) {
  X_signal + matrix(rnorm(length(X_signal), sd = sd_noise), nrow = nrow(X_signal))
})

names(X_true_list) <- block_names
names(X_noisy_list) <- block_names

X_fd1 <- X_noisy_list$FD1_smooth_sparse
X_fd2 <- X_noisy_list$FD2_piecewise
X_rd1 <- X_noisy_list$RD1_sparse
X_rd2 <- X_noisy_list$RD2_dense

X_full_noisy <- do.call(cbind, X_noisy_list)
X_full_true <- do.call(cbind, X_true_list)

###############################################################################
# 4. Define tuning grids
###############################################################################

# Row-direction tuning: this is what makes it two-way.
# Smoothing_parameter in hdClass controls row smoothness for u.
# Sparsity_parameter in hdClass controls row sparsity for u.
smooth_grid_u <- 2^seq(-18, 4, length.out = 12)
sparse_grid_u <- unique(round(seq(0, n - 1, length.out = 16)))

# Column-direction tuning for functional variables
smooth_grid_fd1 <- 2^seq(-18, 4, length.out = 12)
smooth_grid_fd2 <- 2^seq(-18, 4, length.out = 12)

sparse_grid_fd1 <- unique(round(seq(0, m_fd1 - 1, length.out = 12)))
sparse_grid_fd2 <- unique(round(seq(0, m_fd2 - 1, length.out = 12)))

# Column-direction tuning for regular variables:
# regular variables have sparsity only, not smoothness.
sparse_grid_rd1 <- unique(round(seq(0, q_rd1 - 1, length.out = 10)))
sparse_grid_rd2 <- unique(round(seq(0, q_rd2 - 1, length.out = 8)))

###############################################################################
# 5. Create ReMPCA hybrid data objects
###############################################################################

fd_object1 <- fdClass(
  data = as.matrix(X_fd1),
  argval = t_fd1,
  Smoothing_parameter = smooth_grid_fd1,
  Sparsity_parameter = sparse_grid_fd1
)

fd_object2 <- fdClass(
  data = as.matrix(X_fd2),
  argval = t_fd2,
  Smoothing_parameter = smooth_grid_fd2,
  Sparsity_parameter = sparse_grid_fd2
)

rd_object1 <- rdClass(
  data = as.matrix(X_rd1),
  Sparsity_parameter = sparse_grid_rd1
)

rd_object2 <- rdClass(
  data = as.matrix(X_rd2),
  Sparsity_parameter = sparse_grid_rd2
)

hd_object <- hdClass(
  hdlist = list(fd_object1, fd_object2, rd_object1, rd_object2),
  argval = s_grid,
  Smoothing_parameter = smooth_grid_u,
  Sparsity_parameter = sparse_grid_u
)

###############################################################################
# 6. Fit ReMPCA
###############################################################################

# For simulation, centerhds = FALSE is useful because the true data are generated
# from a rank-one signal around zero. For real data, centerhds = TRUE is usually
# more appropriate.
fit_rempca <- ReMPCA(
  hd = hd_object,
  centerhds = FALSE,
  num_pcs = 1,
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
  cv.pick = "1se",
  sparse_tuning_u = NULL,
  sparse_tuning_v = NULL,
  smooth_tuning_u = NULL,
  smooth_tuning_v = NULL
)

print(fit_rempca)

###############################################################################
# 7. Extract estimated scores and loadings
###############################################################################

u_hat_raw <- as.numeric(fit_rempca$PCScores[, 1])

v_hat_raw_list <- list(
  FD1_smooth_sparse = as.numeric(fit_rempca$PCFunctions[[1]][[1]]),
  FD2_piecewise = as.numeric(fit_rempca$PCFunctions[[1]][[2]]),
  RD1_sparse = as.numeric(fit_rempca$PCFunctions[[1]][[3]]),
  RD2_dense = as.numeric(fit_rempca$PCFunctions[[1]][[4]])
)

# Align sign for comparison with truth.
# PCA signs are arbitrary, so this does not change the fitted model.
sign_factor <- ifelse(cor(u_true, u_hat_raw) < 0, -1, 1)

u_hat <- sign_factor * u_hat_raw
v_hat_list <- map(v_hat_raw_list, \(v) sign_factor * v)

# Reconstructed blocks from estimated component
X_recon_list <- map(v_hat_list, \(vhat) {
  outer(u_hat, vhat)
})

names(X_recon_list) <- block_names
X_full_recon <- do.call(cbind, X_recon_list)

###############################################################################
# 8. Tables
###############################################################################

block_table <- tibble(
  Block = block_names,
  Type = block_types,
  Columns = block_dims,
  Noise_SD = noise_sd,
  Signal_Strength = signal_strength,
  Smoothing_Tuned = c(TRUE, TRUE, FALSE, FALSE),
  Sparsity_Tuned = c(TRUE, TRUE, TRUE, TRUE)
)

performance_table <- tibble(
  Block = block_names,
  Type = block_types,
  RMSE_Noisy_vs_True = map2_dbl(X_noisy_list, X_true_list, rmse),
  RMSE_Recon_vs_True = map2_dbl(X_recon_list, X_true_list, rmse),
  Relative_RMSE = RMSE_Recon_vs_True / RMSE_Noisy_vs_True,
  Loading_Correlation = map2_dbl(v_true_list, v_hat_list, abs_cor)
)

score_table <- tibble(
  Metric = c("Score correlation", "Score RMSE after unit normalization"),
  Value = c(
    abs_cor(u_true, u_hat),
    rmse(unit_norm(u_true), unit_norm(u_hat))
  )
)

tuning_table <- tibble(
  Parameter = c(
    "OptimalAlphaU",
    "OptimalGammaU",
    "OptimalAlphaV",
    "OptimalGammaV",
    "VarianceExplained"
  ),
  Value = c(
    collapse_value(safe_extract(fit_rempca, "OptimalAlphaU")),
    collapse_value(safe_extract(fit_rempca, "OptimalGammaU")),
    collapse_value(safe_extract(fit_rempca, "OptimalAlphaV")),
    collapse_value(safe_extract(fit_rempca, "OptimalGammaV")),
    collapse_value(safe_extract(fit_rempca, "VarianceExplained"))
  )
)

write.csv(block_table, "tables/block_summary_table.csv", row.names = FALSE)
write.csv(performance_table, "tables/performance_metrics_table.csv", row.names = FALSE)
write.csv(score_table, "tables/score_metrics_table.csv", row.names = FALSE)
write.csv(tuning_table, "tables/tuning_parameter_table.csv", row.names = FALSE)

# Pretty HTML tables
block_table |>
  gt() |>
  tab_header(
    title = "Simulated Hybrid Data Blocks",
    subtitle = "Functional and regular variables used in the multivariate two-way hybrid PCA simulation"
  ) |>
  gtsave("tables/block_summary_table.html")

performance_table |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  gt() |>
  tab_header(
    title = "ReMPCA Simulation Performance",
    subtitle = "Reconstruction accuracy and loading recovery by hybrid block"
  ) |>
  gtsave("tables/performance_metrics_table.html")

score_table |>
  mutate(Value = round(Value, 4)) |>
  gt() |>
  tab_header(
    title = "Recovered Row Score Performance",
    subtitle = "Comparison of estimated and true row-direction score"
  ) |>
  gtsave("tables/score_metrics_table.html")

tuning_table |>
  gt() |>
  tab_header(
    title = "Selected ReMPCA Tuning Parameters",
    subtitle = "Smoothness, sparsity, and variance explained"
  ) |>
  gtsave("tables/tuning_parameter_table.html")

cat("\nBlock summary table:\n")
print(kable(block_table, digits = 4))

cat("\nPerformance table:\n")
print(kable(performance_table, digits = 4))

cat("\nScore table:\n")
print(kable(score_table, digits = 4))

cat("\nTuning table:\n")
print(kable(tuning_table))

###############################################################################
# 9. Plot 1: true vs estimated row score u
###############################################################################

u_plot_df <- tibble(
  s = rep(s_grid, 2),
  value = c(unit_norm(u_true), unit_norm(u_hat)),
  Curve = rep(c("True u", "Estimated u"), each = n)
)

p_u <- ggplot(u_plot_df, aes(x = s, y = value, linetype = Curve)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "True vs Estimated Row Score",
    subtitle = "Two-way regularization acts on this row-direction component",
    x = "Row domain",
    y = "Unit-normalized score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "figures/01_true_vs_estimated_u.pdf",
  plot = p_u,
  width = 8,
  height = 4.5
)

ggsave(
  filename = "figures/01_true_vs_estimated_u.png",
  plot = p_u,
  width = 8,
  height = 4.5,
  dpi = 300
)

###############################################################################
# 10. Plot 2: true vs estimated loadings v for each hybrid block
###############################################################################

loading_df <- map2_dfr(
  names(v_true_list),
  seq_along(v_true_list),
  \(block_name, i) {
    true_v <- unit_norm(v_true_list[[i]])
    est_v <- unit_norm(v_hat_list[[i]])

    tibble(
      block = block_name,
      index = seq_along(true_v),
      True = true_v,
      Estimated = est_v
    ) |>
      pivot_longer(
        cols = c(True, Estimated),
        names_to = "Curve",
        values_to = "value"
      )
  }
)

p_loading <- ggplot(loading_df, aes(x = index, y = value, linetype = Curve)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.0, alpha = 0.75) +
  facet_wrap(~ block, scales = "free_x", ncol = 2) +
  labs(
    title = "True vs Estimated Hybrid Loadings",
    subtitle = "Functional and regular blocks are estimated jointly",
    x = "Feature index / grid point",
    y = "Unit-normalized loading"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "figures/02_true_vs_estimated_loadings.pdf",
  plot = p_loading,
  width = 10,
  height = 7
)

ggsave(
  filename = "figures/02_true_vs_estimated_loadings.png",
  plot = p_loading,
  width = 10,
  height = 7,
  dpi = 300
)

###############################################################################
# 11. Plot 3: true, noisy, and reconstructed heatmaps
###############################################################################

heat_df <- bind_rows(
  map2_dfr(X_true_list, names(X_true_list), \(M, nm) {
    matrix_to_heat_df(M, nm, "True signal")
  }),
  map2_dfr(X_noisy_list, names(X_noisy_list), \(M, nm) {
    matrix_to_heat_df(M, nm, "Noisy data")
  }),
  map2_dfr(X_recon_list, names(X_recon_list), \(M, nm) {
    matrix_to_heat_df(M, nm, "ReMPCA reconstruction")
  })
)

heat_df$version <- factor(
  heat_df$version,
  levels = c("True signal", "Noisy data", "ReMPCA reconstruction")
)

p_heat <- ggplot(heat_df, aes(x = feature, y = row, fill = value)) +
  geom_raster() +
  scale_y_reverse() +
  facet_grid(block ~ version, scales = "free", space = "free_x") +
  labs(
    title = "Hybrid Data Blocks: True, Noisy, and Reconstructed",
    x = "Feature index / grid point",
    y = "Observation index",
    fill = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = "figures/03_true_noisy_reconstructed_heatmaps.pdf",
  plot = p_heat,
  width = 13,
  height = 9
)

ggsave(
  filename = "figures/03_true_noisy_reconstructed_heatmaps.png",
  plot = p_heat,
  width = 13,
  height = 9,
  dpi = 300
)

###############################################################################
# 12. Plot 4: RMSE improvement table as a bar plot
###############################################################################

rmse_plot_df <- performance_table |>
  select(Block, RMSE_Noisy_vs_True, RMSE_Recon_vs_True) |>
  pivot_longer(
    cols = c(RMSE_Noisy_vs_True, RMSE_Recon_vs_True),
    names_to = "Version",
    values_to = "RMSE"
  ) |>
  mutate(
    Version = recode(
      Version,
      RMSE_Noisy_vs_True = "Noisy vs true",
      RMSE_Recon_vs_True = "Reconstructed vs true"
    )
  )

p_rmse <- ggplot(rmse_plot_df, aes(x = Block, y = RMSE, fill = Version)) +
  geom_col(position = "dodge") +
  labs(
    title = "Noise Reduction by ReMPCA Reconstruction",
    x = "Hybrid block",
    y = "RMSE"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "figures/04_rmse_comparison.pdf",
  plot = p_rmse,
  width = 8,
  height = 5
)

ggsave(
  filename = "figures/04_rmse_comparison.png",
  plot = p_rmse,
  width = 8,
  height = 5,
  dpi = 300
)

###############################################################################
# 13. Plot 5: Perspective plots for functional blocks
###############################################################################

pdf("figures/05_functional_blocks_perspective_plots.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(2, 2, 3, 1))

persp(
  X_true_list$FD1_smooth_sparse,
  theta = 45,
  phi = 25,
  main = "FD1 True",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

persp(
  X_noisy_list$FD1_smooth_sparse,
  theta = 45,
  phi = 25,
  main = "FD1 Noisy",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

persp(
  X_recon_list$FD1_smooth_sparse,
  theta = 45,
  phi = 25,
  main = "FD1 Reconstructed",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

persp(
  X_true_list$FD2_piecewise,
  theta = 45,
  phi = 25,
  main = "FD2 True",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

persp(
  X_noisy_list$FD2_piecewise,
  theta = 45,
  phi = 25,
  main = "FD2 Noisy",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

persp(
  X_recon_list$FD2_piecewise,
  theta = 45,
  phi = 25,
  main = "FD2 Reconstructed",
  xlab = "Rows",
  ylab = "Grid",
  zlab = "Value"
)

dev.off()

###############################################################################
# 14. Built-in ReMPCA plots
###############################################################################

# These are the plots already included in your package.
# They are useful for quick package demonstration and debugging.

pdf("figures/06_builtin_ReMPCA_plots.pdf", width = 11, height = 8.5)

safe_pkg_plot("plot_pc_functions", fit_rempca)
safe_pkg_plot("plot_pc_scores", fit_rempca)
safe_pkg_plot("plot_cv_u", fit_rempca)
safe_pkg_plot("plot_cv_v", fit_rempca)
safe_pkg_plot("plot_gcv_u", fit_rempca)
safe_pkg_plot("plot_gcv_v", fit_rempca)

dev.off()

###############################################################################
# 15. Optional: compare against ordinary SVD baseline
###############################################################################

ordinary_svd <- svd(X_full_noisy)
u_svd <- ordinary_svd$u[, 1] * ordinary_svd$d[1]
v_svd <- ordinary_svd$v[, 1]

# Align sign
svd_sign <- ifelse(cor(u_true, u_svd) < 0, -1, 1)
u_svd <- svd_sign * u_svd
v_svd <- svd_sign * v_svd

X_svd_recon <- outer(u_svd, v_svd)

baseline_table <- tibble(
  Method = c("Noisy data", "Ordinary rank-one SVD", "ReMPCA"),
  RMSE_vs_True = c(
    rmse(X_full_noisy, X_full_true),
    rmse(X_svd_recon, X_full_true),
    rmse(X_full_recon, X_full_true)
  ),
  Score_Correlation = c(
    NA_real_,
    abs_cor(u_true, u_svd),
    abs_cor(u_true, u_hat)
  )
)

write.csv(
  baseline_table,
  "tables/baseline_comparison_table.csv",
  row.names = FALSE
)

baseline_table |>
  mutate(across(where(is.numeric), \(x) round(x, 4))) |>
  gt() |>
  tab_header(
    title = "Baseline Comparison",
    subtitle = "ReMPCA versus ordinary rank-one SVD"
  ) |>
  gtsave("tables/baseline_comparison_table.html")

p_baseline <- ggplot(baseline_table, aes(x = Method, y = RMSE_vs_True)) +
  geom_col() +
  labs(
    title = "Overall Reconstruction Error",
    subtitle = "Lower RMSE indicates better recovery of the true signal",
    x = NULL,
    y = "RMSE vs true signal"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = "figures/07_baseline_comparison.pdf",
  plot = p_baseline,
  width = 7,
  height = 5
)

ggsave(
  filename = "figures/07_baseline_comparison.png",
  plot = p_baseline,
  width = 7,
  height = 5,
  dpi = 300
)

###############################################################################
# 16. Final message
###############################################################################

cat("\nDone!\n")
cat("Figures saved in: figures/\n")
cat("Tables saved in: tables/\n")
cat("\nMain objects in memory:\n")
cat("  hd_object      : hybrid data object\n")
cat("  fit_rempca     : fitted ReMPCA model\n")
cat("  block_table    : simulated block information\n")
cat("  performance_table : per-block recovery metrics\n")
cat("  tuning_table   : selected tuning parameters\n")
cat("  baseline_table : ReMPCA vs ordinary SVD\n")
