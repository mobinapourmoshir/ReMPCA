## -----------------------------------------------------------------------------
## 0) Settings
## -----------------------------------------------------------------------------

N_RUNS    <- 100
BASE_SEED <- 20260609          # dataset seed = BASE_SEED + 10000 * scenario + run

# Power-algorithm stopping rule: stop when ||v_new - v_old||^2 < 1e-10.
# Every converging fit in this study needs at most about 2,800 iterations.
# 115 of the 7,200 fits never reach the tolerance (the iterates cycle between
# nearby sparsity patterns, even after 200,000 iterations); they are stopped
# at MAX_ITER, recorded as not converged (NA), and excluded from Table 2 and
# Figure 2. The number of converged runs per cell is n_success in the summary.
THRESH   <- 1e-10
MAX_ITER <- 10000L

N_CORES <- as.integer(Sys.getenv("REMPCA_CORES", max(1, parallel::detectCores() - 1)))
OUT_DIR <- Sys.getenv("REMPCA_OUT", file.path("inst", "reproducibility", "output_simulation"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# One BLAS thread per worker; the workers already run in parallel.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

## -----------------------------------------------------------------------------
## 1) Packages
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ReMPCA)      # get.pen(): second-order roughness smoother S_alpha
  library(PMA)         # PMD
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

get.pen <- ReMPCA:::get.pen

## -----------------------------------------------------------------------------
## 2) Rank-one regularized power algorithm (Algorithm 1 of the paper)
## -----------------------------------------------------------------------------

# Thresholding rule: tuning_parameter = number of entries set to zero.
sparse_pen <- function(y, tuning_parameter, type = "soft") {
  if (tuning_parameter == 0 || tuning_parameter > length(y)) return(y)
  lambda <- sort(abs(y))[tuning_parameter]
  if (type == "soft") return(sign(y) * pmax(abs(y) - lambda, 0))
  ifelse(abs(y) > lambda, y, 0)                        # hard
}

norm_vec      <- function(x) sqrt(sum(x^2))
normalize_vec <- function(x) {
  x <- as.vector(x); nx <- norm_vec(x)
  if (!is.finite(nx) || nx == 0) return(x)
  x / nx
}
fro_norm <- function(X) sqrt(sum(as.matrix(X)^2))

# X = [X_1 ... X_p]; gamma_u, gamma_v = numbers of zeros in u and in each v_j;
# S_u, S_v = smoothing matrices (identity = no smoothing).
power_algo <- function(X, p_blocks, gamma_u, gamma_v, S_u, S_v,
                       type = "soft", thresh = THRESH, max_iter = MAX_ITER) {
  X <- as.matrix(X)
  idx <- split(seq_len(ncol(X)), rep(seq_along(p_blocks), p_blocks))
  v_old <- svd(X)$v[, 1]
  err <- Inf; iter <- 0L

  while (err > thresh && iter < max_iter) {
    iter <- iter + 1L
    # u-step
    Xv <- X %*% v_old
    u_raw <- sparse_pen(Xv, gamma_u, type)
    if (all(u_raw == 0)) u_raw <- Xv
    u_new <- S_u %*% u_raw
    # v-step, block by block
    v_new <- c()
    for (j in seq_along(idx)) {
      Xu_j <- t(X[, idx[[j]]]) %*% u_new
      v_raw <- sparse_pen(Xu_j, as.integer(gamma_v[j]), type)
      if (all(v_raw == 0)) v_raw <- Xu_j
      v_new <- c(v_new, S_v[[j]] %*% as.vector(v_raw))
    }
    v_new <- v_new / norm_vec(v_new)
    err   <- sum((v_new - v_old)^2)
    v_old <- v_new
  }
  list(v = v_new, u = u_new / norm_vec(u_new), iterations = iter, converged = err <= thresh)
}

## -----------------------------------------------------------------------------
## 3) Simulation design: smooth and sparse truth in both u and v
## -----------------------------------------------------------------------------

n  <- 80                       # rows
p  <- c(35, 45, 40)            # grid points of the three functional blocks
signal_strength <- 4           # X_signal = 4 u v', ||u|| = ||v|| = 1

noise_levels <- data.frame(
  noise_level = c("Low noise", "Medium noise", "High noise"),
  noise_sd    = c(0.12, 0.40, 0.85),
  alpha_mult  = c(0.01, 0.50, 4.00),       # smoothing level grows with the noise
  stringsAsFactors = FALSE
)
sparsity_levels <- data.frame(
  sparsity_level = c("Low sparsity", "Medium sparsity", "High sparsity"),
  zero_prop      = c(0.10, 0.45, 0.72),    # fraction of exactly-zero entries
  stringsAsFactors = FALSE
)
base_alpha <- 2^(-7)

bump <- function(t, a, b) {
  out <- rep(0, length(t)); i <- which(t >= a & t <= b)
  out[i] <- sin(pi * (t[i] - a) / (b - a))^2
  out
}

# Two adjacent bumps of opposite sign covering a fraction 1 - zero_prop of [0, 1].
make_sparse_smooth_vec <- function(m, zero_prop, amp1, amp2, shift) {
  t <- seq(0, 1, length.out = m)
  w   <- max(0.25, min(1 - zero_prop, 0.90))
  gap <- min(0.06, max(0.02, 1 - w - 0.02))
  start1 <- (1 - w - gap) / 2
  end1   <- start1 + 0.55 * w
  start2 <- end1 + gap
  end2   <- start2 + 0.45 * w
  ends <- c(start1, end1, start2, end2) + shift
  if (min(ends) < 0.01) ends <- ends + (0.01 - min(ends))
  if (max(ends) > 0.99) ends <- ends - (max(ends) - 0.99)
  normalize_vec(amp1 * bump(t, ends[1], ends[2]) + amp2 * bump(t, ends[3], ends[4]))
}

make_truth <- function(zero_prop) {
  u  <- make_sparse_smooth_vec(n,    zero_prop,  1.10, -0.55,  0.00)
  v1 <- make_sparse_smooth_vec(p[1], zero_prop,  1.00, -0.65, -0.03)
  v2 <- make_sparse_smooth_vec(p[2], zero_prop, -0.85,  1.10,  0.02)
  v3 <- make_sparse_smooth_vec(p[3], zero_prop,  1.05, -0.75,  0.04)
  v  <- normalize_vec(c(v1, v2, v3))
  idx <- split(seq_along(v), rep(1:3, p))
  list(u = u, v = v, v_blocks = lapply(idx, function(i) v[i]), idx = idx,
       grid_u = seq(0, 1, length.out = n),
       grid_v = lapply(p, function(m) seq(0, 1, length.out = m)))
}

make_dataset <- function(seed, truth, noise_sd) {
  set.seed(seed)
  X_signal <- signal_strength * tcrossprod(truth$u, truth$v)
  X <- X_signal + matrix(rnorm(length(X_signal), sd = noise_sd), nrow(X_signal))
  list(X = X, X_signal = X_signal)
}

## -----------------------------------------------------------------------------
## 4) Methods and tuning parameters
## -----------------------------------------------------------------------------

methods <- data.frame(
  method_id = 1:8,
  method = c("No penalty", "Smooth+Sparse on u", "Smooth+Sparse on v",
             "Two-way Smoothness", "Two-way Sparsity",
             "Two-way Smoothness and Sparsity", "HFV-PCA", "PMD"),
  short  = c("No penalty", "Smooth+Sparse u", "Smooth+Sparse v", "Two-way Smooth",
             "Two-way Sparse", "Two-way Both", "HFV-PCA", "PMD"),
  smooth_u = c(FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE, FALSE),
  sparse_u = c(FALSE, TRUE,  FALSE, FALSE, TRUE,  TRUE,  FALSE, FALSE),
  smooth_v = c(FALSE, FALSE, TRUE,  TRUE,  FALSE, TRUE,  FALSE, FALSE),
  sparse_v = c(FALSE, FALSE, TRUE,  FALSE, TRUE,  TRUE,  FALSE, FALSE),
  stringsAsFactors = FALSE
)

# Tuning parameters are fixed in advance for every scenario (no CV/GCV):
#   smoothing  alpha = 2^-7 * alpha_mult(noise level) * method factor
#   sparsity   gamma = number of exactly-zero entries of the true u / v_j * method factor
# Methods with a single type of penalty use factor 1. For Two-way Smoothness
# and Sparsity the factors below balance the smoothing and sparsity penalties.
both_alpha_factor <- function(noise_level, zero_prop) {
  if (noise_level == "Low noise") {
    if (zero_prop <= 0.15) return(1.00)
    if (zero_prop >= 0.65) return(0.015)
    return(0.20)
  }
  if (noise_level == "Medium noise") 0.60 else 0.90
}
both_gamma_factor <- function(noise_level, zero_prop) {
  if (noise_level == "Low noise") return(if (zero_prop <= 0.15) 0.15 else 0.90)
  if (zero_prop <= 0.15) return(0.25)
  if (zero_prop <= 0.50) return(0.90)
  1.00
}

clip_int <- function(x, lo, hi) as.integer(pmax(lo, pmin(round(x), hi)))

smoother <- function(grid, alpha) {
  if (alpha <= 0) return(diag(length(grid)))
  Re(get.pen(grid, alpha = alpha, type = "Second_order")$S.alpha)
}

# Everything that depends only on the scenario (truth, tuning, smoothers).
build_scenario <- function(sc) {
  truth  <- make_truth(sc$zero_prop)
  zero_u <- sum(abs(truth$u) <= 1e-12)
  zero_v <- vapply(truth$v_blocks, function(z) sum(abs(z) <= 1e-12), integer(1))
  gamma_u_true <- clip_int(zero_u, 0, n - 2)
  gamma_v_true <- clip_int(zero_v, 0, p - 2)
  alpha_sc <- base_alpha * sc$alpha_mult

  specs <- lapply(seq_len(nrow(methods)), function(m) {
    mt <- methods[m, ]
    both <- mt$method == "Two-way Smoothness and Sparsity"
    fa <- if (both) both_alpha_factor(sc$noise_level, sc$zero_prop) else 1
    fg <- if (both) both_gamma_factor(sc$noise_level, sc$zero_prop) else 1
    alpha_u <- if (mt$smooth_u) fa * alpha_sc else 0
    alpha_v <- if (mt$smooth_v) rep(fa * alpha_sc, 3) else rep(0, 3)
    gamma_u <- if (mt$sparse_u) clip_int(fg * gamma_u_true, 0, n - 2) else 0L
    gamma_v <- if (mt$sparse_v) clip_int(fg * gamma_v_true, 0, p - 2) else rep(0L, 3)
    list(method_id = mt$method_id, method = mt$method,
         alpha_u = alpha_u, alpha_v = alpha_v, gamma_u = gamma_u, gamma_v = gamma_v,
         S_u = smoother(truth$grid_u, alpha_u),
         S_v = lapply(1:3, function(j) smoother(truth$grid_v[[j]], alpha_v[j])))
  })
  c(as.list(sc), list(truth = truth, specs = specs))
}

## -----------------------------------------------------------------------------
## 5) Fitting one method on one dataset
## -----------------------------------------------------------------------------

# HFV-PCA (Jang 2021), functional-only case: blockwise FPCA keeping 99% of the
# variance, then a joint PCA of the stacked scores.
fit_hfvpca <- function(X, truth) {
  bl <- lapply(truth$idx, function(idx) {
    s  <- sqrt(length(idx))
    pc <- stats::prcomp(X[, idx, drop = FALSE] / s, center = TRUE, scale. = FALSE)
    k  <- which(cumsum(pc$sdev^2) / sum(pc$sdev^2) >= 0.99)[1]
    list(idx = idx, s = s, pc = pc, k = k)
  })
  scores <- do.call(cbind, lapply(bl, function(b) b$pc$x[, seq_len(b$k), drop = FALSE]))
  w <- svd(scores, nu = 0, nv = 1)$v[, 1]
  v <- numeric(ncol(X)); off <- 0L
  for (b in bl) {
    j <- off + seq_len(b$k)
    v[b$idx] <- drop(b$pc$rotation[, seq_len(b$k), drop = FALSE] %*% w[j]) / b$s
    off <- off + b$k
  }
  v <- normalize_vec(v)
  list(u = drop(X %*% v), v = v, iterations = NA_integer_, converged = TRUE)
}

# PMD (Witten et al. 2009) with L1 penalties on u and v; the L1 bounds are set
# from the true numbers of non-zero entries, like the sparsity levels above.
fit_pmd <- function(X, truth) {
  fit <- PMA::PMD(X, type = "standard", sumabs = NULL,
                  sumabsu = max(1.01, sqrt(sum(abs(truth$u) > 1e-12))),
                  sumabsv = max(1.01, sqrt(sum(abs(truth$v) > 1e-12))),
                  K = 1, niter = 50, center = FALSE, trace = FALSE)
  list(u = fit$u, v = fit$v, iterations = NA_integer_, converged = TRUE)
}

# Rank-one reconstruction d u v' with d = u'Xv; error relative to the true signal.
fit_one <- function(X, X_signal, truth, spec) {
  fit <- switch(spec$method,
                "HFV-PCA" = fit_hfvpca(X, truth),
                "PMD"     = fit_pmd(X, truth),
                power_algo(X, p, spec$gamma_u, spec$gamma_v, spec$S_u, spec$S_v)
  )
  if (!fit$converged)
    return(data.frame(reconstruction_error = NA_real_, iterations = fit$iterations, converged = FALSE))
  u <- normalize_vec(fit$u); v <- normalize_vec(fit$v)
  d <- as.numeric(t(u) %*% X %*% v)
  data.frame(reconstruction_error = fro_norm(X_signal - d * tcrossprod(u, v)) / fro_norm(X_signal),
             iterations = fit$iterations, converged = TRUE)
}

run_job <- function(scenario_id, run) {
  sc   <- scenarios[[scenario_id]]
  seed <- BASE_SEED + 10000 * scenario_id + run
  dat  <- make_dataset(seed, sc$truth, sc$noise_sd)
  do.call(rbind, lapply(sc$specs, function(spec) {
    cbind(data.frame(scenario_id = scenario_id, noise_level = sc$noise_level,
                     sparsity_level = sc$sparsity_level, run = run, seed = seed,
                     method_id = spec$method_id, method = spec$method),
          fit_one(dat$X, dat$X_signal, sc$truth, spec))
  }))
}

## -----------------------------------------------------------------------------
## 6) Run all 9 x 100 datasets (8 methods each)
## -----------------------------------------------------------------------------

design <- expand.grid(noise_id = 1:3, sparsity_id = 1:3)
design <- design[order(design$noise_id, design$sparsity_id), ]
design <- cbind(noise_levels[design$noise_id, ], sparsity_levels[design$sparsity_id, ])
rownames(design) <- NULL
scenarios <- lapply(seq_len(nrow(design)), function(i) build_scenario(design[i, ]))

jobs <- expand.grid(run = seq_len(N_RUNS), scenario_id = seq_along(scenarios))
cat(sprintf("Running %d datasets x %d methods on %d cores ...\n", nrow(jobs), nrow(methods), N_CORES))
t0 <- Sys.time()

job_fun <- function(k) run_job(jobs$scenario_id[k], jobs$run[k])
if (N_CORES > 1) {
  cl <- makeCluster(N_CORES)
  invisible(clusterEvalQ(cl, Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")))
  clusterExport(cl, setdiff(ls(globalenv()), "cl"), envir = globalenv())
  res <- parLapplyLB(cl, seq_len(nrow(jobs)), job_fun)
  stopCluster(cl)
} else {
  res <- lapply(seq_len(nrow(jobs)), job_fun)
}
results <- do.call(rbind, res)
cat(sprintf("Done in %.1f min. Not converged within %d iterations: %d of %d fits.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), MAX_ITER,
            sum(!results$converged), nrow(results)))
write.csv(results, file.path(OUT_DIR, "reconstruction_errors_9scenarios_all_runs.csv"), row.names = FALSE)

tuning <- do.call(rbind, lapply(scenarios, function(sc) do.call(rbind, lapply(sc$specs, function(s)
  data.frame(noise_level = sc$noise_level, sparsity_level = sc$sparsity_level, method = s$method,
             alpha_u = s$alpha_u, alpha_v = s$alpha_v[1], gamma_u = s$gamma_u,
             gamma_v1 = s$gamma_v[1], gamma_v2 = s$gamma_v[2], gamma_v3 = s$gamma_v[3])))))
write.csv(tuning, file.path(OUT_DIR, "tuning_parameters_9scenarios.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## 7) Table 2
## -----------------------------------------------------------------------------

lev_noise <- noise_levels$noise_level
lev_spars <- sparsity_levels$sparsity_level

summary_tab <- results |>
  filter(converged) |>
  group_by(noise_level, sparsity_level, method_id, method) |>
  summarise(n_success = n(), mean = mean(reconstruction_error), sd = sd(reconstruction_error),
            median = median(reconstruction_error), .groups = "drop") |>
  mutate(noise_level = factor(noise_level, lev_noise),
         sparsity_level = factor(sparsity_level, lev_spars)) |>
  arrange(noise_level, method_id, sparsity_level)
write.csv(summary_tab, file.path(OUT_DIR, "reconstruction_error_summary_9scenarios.csv"), row.names = FALSE)

table_label <- c("No penalty", "Smooth+Sparse $\\bm{u}$", "Smooth+Sparse $\\bm{v}$",
                 "Two-way Smooth", "Two-way Sparse", "Two-way Smooth+Sparse", "HFV-PCA", "PMD")
cells <- summary_tab |>
  group_by(noise_level, sparsity_level) |>
  mutate(cell = sprintf("%.3f (%.3f)", mean, sd),
         cell = ifelse(mean == min(mean), paste0("\\textbf{", cell, "}"), cell)) |>
  ungroup()
tex <- c("\\begin{table}[!ht]", "\\centering", "\\scriptsize",
         "\\caption{Mean relative reconstruction error over simulation runs for the nine noise--sparsity scenarios. Standard deviations are shown in parentheses.}",
         "\\label{tab:reconstruction-error-9scenarios}",
         "\\resizebox{.8\\textwidth}{!}{%", "\\begin{tabular}{llccc}", "\\toprule",
         "Noise level & Method & Low sparsity & Medium sparsity & High sparsity \\\\", "\\midrule")
for (nl in lev_noise) {
  tex <- c(tex, sprintf("\\multirow{8}{*}{%s}", nl))
  for (m in 1:8) {
    x <- cells$cell[cells$noise_level == nl & cells$method_id == m]
    tex <- c(tex, sprintf("& %s\n& %s & %s & %s \\\\", table_label[m], x[1], x[2], x[3]))
  }
  tex <- c(tex, if (nl != tail(lev_noise, 1)) "\\midrule" else "\\bottomrule")
}
writeLines(c(tex, "\\end{tabular}}", "\\end{table}"),
           file.path(OUT_DIR, "reconstruction_error_latex_table_9scenarios.tex"))

## -----------------------------------------------------------------------------
## 8) Check against the published Table 2 (mean and sd, 3 decimals)
## -----------------------------------------------------------------------------

# Order: noise level (low, medium, high) > method 1-8 > sparsity (low, medium, high)
paper_mean <- c(
  0.449, 0.451, 0.449,  0.344, 0.352, 0.476,  0.314, 0.340, 0.504,  0.188, 0.275, 0.566,
  0.442, 0.391, 0.298,  0.186, 0.196, 0.195,  0.472, 0.468, 0.462,  0.449, 0.451, 0.388,
  2.081, 2.079, 2.066,  1.201, 1.310, 1.439,  1.194, 1.400, 1.456,  0.572, 0.831, 1.127,
  2.074, 1.948, 1.606,  0.553, 0.767, 1.052,  2.092, 2.072, 2.051,  2.081, 2.061, 1.819,
  4.251, 4.259, 4.251,  2.649, 2.664, 2.605,  2.551, 2.522, 2.416,  1.213, 1.348, 1.363,
  4.228, 3.913, 3.217,  1.213, 1.327, 1.275,  4.208, 4.213, 4.205,  4.251, 4.225, 3.779)
paper_sd <- c(
  0.023, 0.026, 0.021,  0.020, 0.021, 0.017,  0.021, 0.020, 0.016,  0.020, 0.018, 0.012,
  0.023, 0.028, 0.025,  0.020, 0.022, 0.026,  0.026, 0.025, 0.021,  0.023, 0.025, 0.018,
  0.106, 0.103, 0.098,  0.084, 0.114, 0.085,  0.110, 0.113, 0.043,  0.068, 0.053, 0.036,
  0.104, 0.093, 0.149,  0.068, 0.061, 0.053,  0.093, 0.099, 0.093,  0.106, 0.100, 0.142,
  0.070, 0.076, 0.071,  0.124, 0.133, 0.146,  0.120, 0.105, 0.122,  0.130, 0.090, 0.079,
  0.070, 0.073, 0.108,  0.132, 0.093, 0.071,  0.072, 0.077, 0.073,  0.070, 0.073, 0.071)

check <- summary_tab |>
  mutate(paper_mean = paper_mean, paper_sd = paper_sd,
         match = sprintf("%.3f", mean) == sprintf("%.3f", paper_mean) &
           sprintf("%.3f", sd)   == sprintf("%.3f", paper_sd))
write.csv(check, file.path(OUT_DIR, "table2_check.csv"), row.names = FALSE)
cat(sprintf("\nTable 2 check: %d of %d cells (mean and sd, 3 decimals) equal the published values.\n",
            sum(check$match), nrow(check)))

## -----------------------------------------------------------------------------
## 9) Figure 2
## -----------------------------------------------------------------------------

method_cols <- c("No penalty" = "#E64B35", "Smooth+Sparse u" = "#D9A300",
                 "Smooth+Sparse v" = "#4DBBD5", "Two-way Smooth" = "#00A087",
                 "Two-way Sparse" = "#7E57C2", "Two-way Both" = "#3C5488",
                 "HFV-PCA" = "#9C755F", "PMD" = "#E07BA8")

plot_data <- results |>
  filter(converged) |>
  mutate(method_short = factor(methods$short[method_id], levels = methods$short),
         noise_level = factor(noise_level, lev_noise),
         sparsity_level = factor(sparsity_level, lev_spars,
                                 labels = c("Low\nsparsity", "Medium\nsparsity", "High\nsparsity")))

fig2 <- ggplot(plot_data, aes(x = sparsity_level, y = reconstruction_error, fill = method_short)) +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dashed", colour = "grey75", linewidth = 0.45) +
  geom_boxplot(position = position_dodge2(width = 0.85, preserve = "single"),
               width = 0.70, alpha = 0.90, colour = "grey25", linewidth = 0.35,
               outlier.shape = 16, outlier.size = 0.45, outlier.alpha = 0.55) +
  facet_wrap(~ noise_level, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = method_cols) +
  scale_y_continuous(trans = "log10", expand = expansion(mult = c(0, 0)),
                     labels = function(x) format(x, digits = 3, scientific = FALSE)) +
  labs(x = NULL, y = "Relative reconstruction error, log scale", fill = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_grey(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
        legend.text = element_text(size = 17), legend.key.size = unit(0.65, "cm"),
        strip.text.x = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 17, face = "bold", colour = "black", hjust = 0.5,
                                   vjust = 0.5, lineheight = 0.9),
        axis.text.y = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 19, face = "bold"),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linewidth = 0.45),
        panel.grid.minor.y = element_line(colour = "white", linewidth = 0.25),
        panel.spacing.x = unit(0.8, "lines"), plot.margin = margin(8, 12, 8, 12))

ggsave(file.path(OUT_DIR, "reconstruction_error_boxplot_1x3new.pdf"), fig2,
       width = 15, height = 7.1, units = "in")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_simulation.txt"))
cat("Outputs written to", normalizePath(OUT_DIR), "\n")
