############################################################
## Simulation designed to show WHEN two-way smooth+sparse HPCA wins
## Copy-paste and run
############################################################

rm(list = ls())
set.seed(20260510)

############################################################
## 0. Packages
############################################################

need <- c("remotes", "ggplot2", "dplyr", "tidyr", "patchwork")
to_install <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

if (!requireNamespace("ReMPCA", quietly = TRUE)) {
  remotes::install_github(
    "mobinapourmoshir/ReMPCA",
    dependencies = TRUE,
    upgrade = "never"
  )
}

library(ReMPCA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

############################################################
## 1. User controls
############################################################

## For a quick check, use N <- 250 and NREP <- 1.
## For manuscript-style results, use N <- 600 or 1000 and NREP <- 10+.
N <- 350
NREP <- 1
M <- 3

## Run all four scenarios to show the contrast.
## To run only the key scenario, set:
## SCENARIOS_TO_RUN <- "D. both: smooth plus sparse"
SCENARIOS_TO_RUN <- c(
  "A. dense clean baseline",
  "B. smoothness needed",
  "C. sparsity needed",
  "D. both: smooth plus sparse"
)

## Main scenario to use for the Figure-6.3-style estimated-vs-true plot
PLOT_SCENARIO <- "D. both: smooth plus sparse"

## Important:
## hard thresholding has much less amplitude shrinkage than soft thresholding.
## This usually gives a fairer reconstruction comparison for score sparsity.
SPARSE_TYPE <- "hard"

## Use "min", not "1se", because "1se" can over-regularize in this package.
CV_PICK <- "min"

############################################################
## 2. Helper functions
############################################################

normalize_vec <- function(x) {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x^2))
  if (!is.finite(nrm) || nrm == 0) return(x)
  x / nrm
}

safe_center <- function(X) {
  X <- as.matrix(X)
  sweep(X, 2, colMeans(X), "-")
}

split_matrix_blocks <- function(X, block_sizes) {
  ends <- cumsum(block_sizes)
  starts <- c(1, head(ends, -1) + 1)
  out <- Map(function(a, b) X[, a:b, drop = FALSE], starts, ends)
  names(out) <- names(block_sizes)
  out
}

split_vector_blocks <- function(x, block_sizes) {
  ends <- cumsum(block_sizes)
  starts <- c(1, head(ends, -1) + 1)
  out <- Map(function(a, b) x[a:b], starts, ends)
  names(out) <- names(block_sizes)
  out
}

block_indices <- function(block_sizes) {
  ends <- cumsum(block_sizes)
  starts <- c(1, head(ends, -1) + 1)
  out <- Map(function(a, b) a:b, starts, ends)
  names(out) <- names(block_sizes)
  out
}

permn_base <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permn_base(x[-i]))
  }))
}

align_columns <- function(Vhat, Vtrue) {
  M <- ncol(Vtrue)
  perms <- permn_base(seq_len(M))

  best_score <- -Inf
  best_perm <- seq_len(M)
  best_signs <- rep(1, M)

  for (r in seq_len(nrow(perms))) {
    pp <- perms[r, ]
    signs <- rep(1, M)
    score <- 0
    for (m in seq_len(M)) {
      cc <- sum(Vtrue[, m] * Vhat[, pp[m]], na.rm = TRUE)
      signs[m] <- ifelse(cc < 0, -1, 1)
      score <- score + abs(cc)
    }
    if (score > best_score) {
      best_score <- score
      best_perm <- pp
      best_signs <- signs
    }
  }

  Vout <- Vhat[, best_perm, drop = FALSE]
  for (m in seq_len(M)) Vout[, m] <- best_signs[m] * Vout[, m]

  list(V = Vout, perm = best_perm, signs = best_signs)
}

extract_v_hat <- function(fit, block_sizes, num_pcs) {
  D <- sum(block_sizes)
  Vhat <- matrix(NA_real_, nrow = D, ncol = num_pcs)

  for (m in seq_len(num_pcs)) {
    vv <- unlist(lapply(fit$PCFunctions[[m]], function(z) as.numeric(z)))
    Vhat[, m] <- normalize_vec(vv)
  }

  Vhat
}

reconstruct_from_fit <- function(fit, num_pcs) {
  Reduce("+", lapply(seq_len(num_pcs), function(k) {
    as.matrix(fit$ReconstructedData[[k]])
  }))
}

mrae <- function(Xhat, Xclean_centered) {
  num <- sqrt(rowSums((Xhat - Xclean_centered)^2))
  den <- sqrt(rowSums(Xclean_centered^2))
  mean(num / pmax(den, 1e-12))
}

projection_reconstruction <- function(Xobs_centered, Vhat) {
  Vhat <- as.matrix(Vhat)
  G <- crossprod(Vhat)
  ridge <- 1e-8 * diag(ncol(G))
  Xobs_centered %*% Vhat %*% solve(G + ridge) %*% t(Vhat)
}

score_f1 <- function(true_active, estimated_active) {
  TP <- sum(true_active & estimated_active)
  FP <- sum(!true_active & estimated_active)
  FN <- sum(true_active & !estimated_active)

  precision <- TP / max(TP + FP, 1)
  recall <- TP / max(TP + FN, 1)

  if ((precision + recall) == 0) return(0)
  2 * precision * recall / (precision + recall)
}

make_rough_fd_noise <- function(N, t, white_sd, rough_sd) {
  J <- length(t)

  white <- matrix(rnorm(N * J, sd = white_sd), N, J)

  A1 <- rnorm(N, sd = rough_sd)
  A2 <- rnorm(N, sd = rough_sd)
  A3 <- rnorm(N, sd = rough_sd)
  ph1 <- runif(N, 0, 2 * pi)
  ph2 <- runif(N, 0, 2 * pi)
  ph3 <- runif(N, 0, 2 * pi)

  rough <- matrix(0, N, J)
  for (i in seq_len(N)) {
    rough[i, ] <-
      A1[i] * sin(12 * pi * t + ph1[i]) +
      A2[i] * cos(18 * pi * t + ph2[i]) +
      A3[i] * sin(24 * pi * t + ph3[i])
  }

  white + rough
}

make_image_noise <- function(N, img_nrow, img_ncol, white_sd, rough_sd) {
  d <- img_nrow * img_ncol
  white <- matrix(rnorm(N * d, sd = white_sd), N, d)

  xx <- seq(-1, 1, length.out = img_ncol)
  yy <- seq(-1, 1, length.out = img_nrow)
  grd <- expand.grid(y = yy, x = xx)

  checker1 <- sin(5 * pi * grd$x) * cos(5 * pi * grd$y)
  checker2 <- cos(7 * pi * grd$x + 0.4) * sin(6 * pi * grd$y)
  checker1 <- normalize_vec(checker1)
  checker2 <- normalize_vec(checker2)

  A1 <- rnorm(N, sd = rough_sd)
  A2 <- rnorm(N, sd = rough_sd)

  rough <- A1 %*% t(checker1) + A2 %*% t(checker2)

  white + rough
}

############################################################
## 3. True hybrid loading system
############################################################

make_true_loadings <- function(t1, t2, t3,
                               d_reg,
                               img_nrow,
                               img_ncol,
                               M = 3) {

  ## Smooth functional blocks.
  FD1 <- cbind(
    sin(pi * t1),
    sin(2 * pi * t1),
    sin(3 * pi * t1)
  )

  FD2 <- cbind(
    2 * t2 - 1,
    cos(2 * pi * t2),
    sin(4 * pi * t2)
  )

  FD3 <- cbind(
    exp(-8 * (t3 - 0.25)^2) - 0.6 * exp(-8 * (t3 - 0.75)^2),
    sin(2 * pi * t3) * exp(-0.6 * t3),
    plogis(12 * (t3 - 0.55)) - 0.5
  )

  ## Sparse regular-data block.
  RD <- matrix(0, d_reg, M)
  RD[1:3, 1] <- c(1.0, -0.7, 0.5)
  RD[5:7, 2] <- c(0.8, 1.0, -0.5)
  RD[9:12, 3] <- c(-0.6, 0.7, 0.9, -0.4)

  ## Sparse local image patterns.
  xx <- seq(-1, 1, length.out = img_ncol)
  yy <- seq(-1, 1, length.out = img_nrow)
  grd <- expand.grid(y = yy, x = xx)

  blob <- function(cx, cy, sx = 0.22, sy = 0.22) {
    exp(-((grd$x - cx)^2 / sx^2 + (grd$y - cy)^2 / sy^2))
  }

  pat1 <- blob(-0.45, 0.40) - 0.65 * blob(0.05, 0.40)
  pat2 <- blob(0.45, -0.15) - 0.75 * blob(0.45, -0.65)
  pat3 <- blob(-0.40, -0.40) - blob(0.40, 0.40)

  ## Force image sparsity by zeroing small background values.
  zero_small <- function(z, q = 0.65) {
    z[abs(z) < quantile(abs(z), q)] <- 0
    z
  }

  IMG <- cbind(
    zero_small(pat1),
    zero_small(pat2),
    zero_small(pat3)
  )

  ## Normalize each block contribution before stacking.
  ## These weights control how much each data type contributes to the PCs.
  block_weight <- c(FD1 = 1.00, FD2 = 1.00, FD3 = 0.85, RD = 0.75, IMG = 1.00)

  blocks <- list(FD1 = FD1, FD2 = FD2, FD3 = FD3, RD = RD, IMG = IMG)

  for (bn in names(blocks)) {
    for (m in seq_len(M)) {
      blocks[[bn]][, m] <- block_weight[bn] * normalize_vec(blocks[[bn]][, m])
    }
  }

  Vtrue <- do.call(rbind, blocks)
  Vtrue <- apply(Vtrue, 2, normalize_vec)
  Vtrue <- as.matrix(Vtrue)

  list(
    Vtrue = Vtrue,
    blocks = blocks
  )
}

############################################################
## 4. Score models
############################################################

make_scores <- function(N, M, score_mode) {
  if (score_mode == "dense") {
    group <- rep("Dense", N)
    Scores <- matrix(rnorm(N * M), N, M)
    Scores[, 1] <- 1.2 * Scores[, 1]
    Scores[, 2] <- 0.9 * Scores[, 2]
    Scores[, 3] <- 0.7 * Scores[, 3]
    return(list(Scores = Scores, group = factor(group)))
  }

  ## Figure-6.3-like sparse score structure, with a noise-only group.
  probs <- c(S0 = 0.22, S1 = 0.24, S2 = 0.27, S3 = 0.27)
  Nk <- round(N * probs)
  Nk[length(Nk)] <- N - sum(Nk[-length(Nk)])

  group <- rep(names(Nk), times = Nk)

  lambda_by_group <- rbind(
    S0 = c(0.0, 0.0, 0.0),
    S1 = c(1.4, 0.0, 0.7),
    S2 = c(0.0, 1.2, 0.8),
    S3 = c(1.0, 0.9, 0.0)
  )

  Scores <- matrix(0, nrow = N, ncol = M)
  for (i in seq_len(N)) {
    g <- group[i]
    for (m in seq_len(M)) {
      lam <- lambda_by_group[g, m]
      if (lam > 0) Scores[i, m] <- rnorm(1, sd = sqrt(lam))
    }
  }

  list(Scores = Scores, group = factor(group, levels = names(Nk)))
}

############################################################
## 5. Scenario definitions
############################################################

scenario_pars <- list(
  "A. dense clean baseline" = list(
    score_mode = "dense",
    fd_white_sd = 0.10,
    fd_rough_sd = 0.00,
    rd_sd = 0.08,
    img_white_sd = 0.05,
    img_rough_sd = 0.00
  ),
  "B. smoothness needed" = list(
    score_mode = "dense",
    fd_white_sd = 0.20,
    fd_rough_sd = 0.35,
    rd_sd = 0.10,
    img_white_sd = 0.07,
    img_rough_sd = 0.25
  ),
  "C. sparsity needed" = list(
    score_mode = "sparse",
    fd_white_sd = 0.12,
    fd_rough_sd = 0.00,
    rd_sd = 0.18,
    img_white_sd = 0.11,
    img_rough_sd = 0.00
  ),
  "D. both: smooth plus sparse" = list(
    score_mode = "sparse",
    fd_white_sd = 0.22,
    fd_rough_sd = 0.45,
    rd_sd = 0.18,
    img_white_sd = 0.12,
    img_rough_sd = 0.35
  )
)

scenario_pars <- scenario_pars[SCENARIOS_TO_RUN]

############################################################
## 6. Simulate one hybrid dataset
############################################################

simulate_hybrid_dataset <- function(N, M, pars) {

  J1 <- 90
  J2 <- 90
  J3 <- 70

  t1 <- seq(0, 1, length.out = J1)
  t2 <- seq(0, 1, length.out = J2)
  t3 <- seq(0, 1, length.out = J3)

  d_reg <- 22
  img_nrow <- 14
  img_ncol <- 14
  d_img <- img_nrow * img_ncol

  block_sizes <- c(
    FD1 = J1,
    FD2 = J2,
    FD3 = J3,
    RD = d_reg,
    IMG = d_img
  )

  true_obj <- make_true_loadings(
    t1 = t1,
    t2 = t2,
    t3 = t3,
    d_reg = d_reg,
    img_nrow = img_nrow,
    img_ncol = img_ncol,
    M = M
  )

  Vtrue <- true_obj$Vtrue

  sc <- make_scores(N, M, pars$score_mode)
  Scores <- sc$Scores
  group <- sc$group

  Xclean <- Scores %*% t(Vtrue)
  Xclean_blocks <- split_matrix_blocks(Xclean, block_sizes)

  ## Add noise.
  X1 <- Xclean_blocks$FD1 +
    make_rough_fd_noise(N, t1, pars$fd_white_sd, pars$fd_rough_sd)

  X2 <- Xclean_blocks$FD2 +
    make_rough_fd_noise(N, t2, pars$fd_white_sd, pars$fd_rough_sd)

  X3 <- Xclean_blocks$FD3 +
    make_rough_fd_noise(N, t3, pars$fd_white_sd, pars$fd_rough_sd)

  Xreg <- Xclean_blocks$RD +
    matrix(rnorm(N * d_reg, sd = pars$rd_sd), N, d_reg)

  Ximg <- Xclean_blocks$IMG +
    make_image_noise(
      N = N,
      img_nrow = img_nrow,
      img_ncol = img_ncol,
      white_sd = pars$img_white_sd,
      rough_sd = pars$img_rough_sd
    )

  img_list <- lapply(seq_len(N), function(i) {
    matrix(Ximg[i, ], nrow = img_nrow, ncol = img_ncol, byrow = TRUE)
  })

  Xnoisy <- cbind(X1, X2, X3, Xreg, Ximg)

  list(
    X1 = X1,
    X2 = X2,
    X3 = X3,
    Xreg = Xreg,
    Ximg = Ximg,
    img_list = img_list,
    Xclean = Xclean,
    Xnoisy = Xnoisy,
    Xclean_centered = safe_center(Xclean),
    Xnoisy_centered = safe_center(Xnoisy),
    Vtrue = Vtrue,
    Scores = Scores,
    group = group,
    t1 = t1,
    t2 = t2,
    t3 = t3,
    block_sizes = block_sizes,
    img_nrow = img_nrow,
    img_ncol = img_ncol
  )
}

############################################################
## 7. Build ReMPCA objects
############################################################

make_hd_object <- function(dat) {

  N <- nrow(dat$X1)
  d_reg <- dat$block_sizes["RD"]
  d_img <- dat$block_sizes["IMG"]

  ## Rich enough to let CV choose low or moderate smoothing.
  alpha_fd <- 2^seq(-14, -4, length.out = 6)
  alpha_img <- 2^seq(-10, -3, length.out = 5)

  ## Include zero so the method can choose no sparsity when sparsity is unnecessary.
  sparse_u_grid <- sort(unique(round(c(
    0,
    0.04 * N,
    0.10 * N,
    0.18 * N,
    0.28 * N,
    0.40 * N,
    0.55 * N
  ))))

  sparse_v_reg <- 0:min(12, d_reg - 1)

  sparse_v_img <- sort(unique(round(c(
    0,
    0.15 * d_img,
    0.35 * d_img,
    0.55 * d_img,
    0.70 * d_img,
    0.82 * d_img
  ))))

  fd1_obj <- fdClass(
    data = as.matrix(dat$X1),
    argval = dat$t1,
    Smoothing_parameter = alpha_fd,
    Sparsity_parameter = 0
  )

  fd2_obj <- fdClass(
    data = as.matrix(dat$X2),
    argval = dat$t2,
    Smoothing_parameter = alpha_fd,
    Sparsity_parameter = 0
  )

  fd3_obj <- fdClass(
    data = as.matrix(dat$X3),
    argval = dat$t3,
    Smoothing_parameter = alpha_fd,
    Sparsity_parameter = 0
  )

  rd_obj <- rdClass(
    data = as.matrix(dat$Xreg),
    Sparsity_parameter = sparse_v_reg
  )

  img_obj <- imgClass(
    image = dat$img_list,
    Smoothing_parameter = alpha_img,
    Sparsity_parameter = sparse_v_img
  )

  hd_obj <- hdClass(
    hdlist = list(fd1_obj, fd2_obj, fd3_obj, rd_obj, img_obj),
    argval = seq(0, 1, length.out = N),
    Smoothing_parameter = 0,
    Sparsity_parameter = sparse_u_grid
  )

  n_var <- attr(hd_obj, "n_var")

  list(
    hd_obj = hd_obj,
    n_var = n_var,
    alpha_fd = alpha_fd,
    alpha_img = alpha_img,
    sparse_u_grid = sparse_u_grid,
    sparse_v_reg = sparse_v_reg,
    sparse_v_img = sparse_v_img,
    zero_smooth_v = rep(list(0), n_var),
    zero_sparse_v = rep(list(0), n_var),
    smooth_v = list(alpha_fd, alpha_fd, alpha_fd, 0, alpha_img),
    sparse_v = list(0, 0, 0, sparse_v_reg, sparse_v_img)
  )
}

############################################################
## 8. Fit four methods
############################################################

fit_rempca <- function(hd_stuff,
                       model_name,
                       sparse_u,
                       smooth_v,
                       sparse_v,
                       M,
                       nfolds = 3) {

  message("\n==========================================")
  message("Fitting ", model_name)
  message("==========================================")

  ReMPCA(
    hd = hd_stuff$hd_obj,
    centerhds = TRUE,
    num_pcs = M,
    nfolds_u = nfolds,
    nfolds_v = rep(nfolds, hd_stuff$n_var),
    thresh = 1e-8,
    maxit = 250,
    tuning_iter = 1,
    parallel = FALSE,
    weights = 0,
    smoothness_type = "Second_order",
    sparse_tuning_type = SPARSE_TYPE,
    tuning_order = "Sparsity",
    cv.pick = CV_PICK,
    sparse_tuning_u = sparse_u,
    sparse_tuning_v = sparse_v,
    smooth_tuning_u = 0,
    smooth_tuning_v = smooth_v
  )
}

fit_all_methods <- function(dat, M) {

  hd_stuff <- make_hd_object(dat)

  fit_mfpca <- fit_rempca(
    hd_stuff = hd_stuff,
    model_name = "MFPCA",
    sparse_u = 0,
    smooth_v = hd_stuff$zero_smooth_v,
    sparse_v = hd_stuff$zero_sparse_v,
    M = M
  )

  fit_smooth <- fit_rempca(
    hd_stuff = hd_stuff,
    model_name = "Smoothed HPCA",
    sparse_u = 0,
    smooth_v = hd_stuff$smooth_v,
    sparse_v = hd_stuff$zero_sparse_v,
    M = M
  )

  fit_sparse <- fit_rempca(
    hd_stuff = hd_stuff,
    model_name = "Sparse HPCA",
    sparse_u = hd_stuff$sparse_u_grid,
    smooth_v = hd_stuff$zero_smooth_v,
    sparse_v = hd_stuff$sparse_v,
    M = M
  )

  fit_tw <- fit_rempca(
    hd_stuff = hd_stuff,
    model_name = "Two-way smooth+sparse HPCA",
    sparse_u = hd_stuff$sparse_u_grid,
    smooth_v = hd_stuff$smooth_v,
    sparse_v = hd_stuff$sparse_v,
    M = M
  )

  list(
    "MFPCA" = fit_mfpca,
    "Smoothed HPCA" = fit_smooth,
    "Sparse HPCA" = fit_sparse,
    "Two-way smooth+sparse HPCA" = fit_tw
  )
}

############################################################
## 9. Evaluate one scenario/repetition
############################################################

evaluate_fits <- function(dat, fits, scenario_name, rep_id, M) {

  Vtrue <- dat$Vtrue
  true_active <- abs(dat$Scores) > 1e-12

  out_error <- list()
  out_recon <- list()
  out_zero <- list()
  out_tuning <- list()
  aligned_v <- list()

  for (nm in names(fits)) {
    fit <- fits[[nm]]

    Vhat_raw <- extract_v_hat(fit, dat$block_sizes, M)
    al <- align_columns(Vhat_raw, Vtrue)
    Vhat <- al$V
    aligned_v[[nm]] <- Vhat

    l2 <- sapply(seq_len(M), function(m) {
      sqrt(sum((Vhat[, m] - Vtrue[, m])^2))
    })

    corr <- sapply(seq_len(M), function(m) {
      abs(sum(Vhat[, m] * Vtrue[, m]) /
            sqrt(sum(Vhat[, m]^2) * sum(Vtrue[, m]^2)))
    })

    out_error[[nm]] <- data.frame(
      Scenario = scenario_name,
      Rep = rep_id,
      Method = nm,
      PC = paste0("PC", seq_len(M)),
      L2_error = l2,
      Abs_correlation = corr
    )

    Xhat_raw <- reconstruct_from_fit(fit, M)

    Xhat_projection <- projection_reconstruction(
      Xobs_centered = dat$Xnoisy_centered,
      Vhat = Vhat
    )

    out_recon[[nm]] <- data.frame(
      Scenario = scenario_name,
      Rep = rep_id,
      Method = nm,
      Raw_MRAE_clean = mrae(Xhat_raw, dat$Xclean_centered),
      Projection_MRAE_clean = mrae(Xhat_projection, dat$Xclean_centered),
      Mean_PC_L2_error = mean(l2),
      Mean_PC_abs_correlation = mean(corr)
    )

    U <- as.matrix(fit$PCScores)
    U_aligned <- U[, al$perm, drop = FALSE]
    for (m in seq_len(M)) U_aligned[, m] <- al$signs[m] * U_aligned[, m]

    score_zero_rate <- colMeans(abs(U_aligned) < 1e-8)
    true_zero_rate <- colMeans(!true_active)

    support_f1 <- sapply(seq_len(M), function(m) {
      score_f1(
        true_active = true_active[, m],
        estimated_active = abs(U_aligned[, m]) >= 1e-8
      )
    })

    out_zero[[nm]] <- data.frame(
      Scenario = scenario_name,
      Rep = rep_id,
      Method = nm,
      PC = paste0("PC", seq_len(M)),
      True_zero_rate = true_zero_rate,
      Estimated_zero_rate = score_zero_rate,
      Score_support_F1 = support_f1
    )

    out_tuning[[nm]] <- data.frame(
      Scenario = scenario_name,
      Rep = rep_id,
      Method = nm,
      PC = paste0("PC", seq_len(M)),
      OptimalGammaU = sapply(fit$OptimalGammaU, function(z) paste(z, collapse = ",")),
      OptimalAlphaU = sapply(fit$OptimalAlphaU, function(z) paste(z, collapse = ",")),
      OptimalGammaV = sapply(fit$OptimalGammaV, function(z) paste(unlist(z), collapse = ",")),
      OptimalAlphaV = sapply(fit$OptimalAlphaV, function(z) paste(unlist(z), collapse = ","))
    )
  }

  list(
    error = bind_rows(out_error),
    recon = bind_rows(out_recon),
    zero = bind_rows(out_zero),
    tuning = bind_rows(out_tuning),
    aligned_v = aligned_v
  )
}

############################################################
## 10. Run simulations
############################################################

all_errors <- list()
all_recon <- list()
all_zero <- list()
all_tuning <- list()

saved_plot_object <- NULL
saved_fits <- NULL
saved_eval <- NULL

counter <- 1

for (sc_name in names(scenario_pars)) {
  for (rr in seq_len(NREP)) {

    message("\n\n##################################################")
    message("Scenario: ", sc_name, " | Rep: ", rr)
    message("##################################################")

    dat <- simulate_hybrid_dataset(
      N = N,
      M = M,
      pars = scenario_pars[[sc_name]]
    )

    fits <- fit_all_methods(dat, M = M)

    ev <- evaluate_fits(
      dat = dat,
      fits = fits,
      scenario_name = sc_name,
      rep_id = rr,
      M = M
    )

    all_errors[[counter]] <- ev$error
    all_recon[[counter]] <- ev$recon
    all_zero[[counter]] <- ev$zero
    all_tuning[[counter]] <- ev$tuning

    if (sc_name == PLOT_SCENARIO && rr == 1) {
      saved_plot_object <- dat
      saved_fits <- fits
      saved_eval <- ev
    }

    counter <- counter + 1
  }
}

error_table <- bind_rows(all_errors)
recon_table <- bind_rows(all_recon)
zero_table <- bind_rows(all_zero)
tuning_table <- bind_rows(all_tuning)

summary_table <- recon_table %>%
  group_by(Scenario, Method) %>%
  summarise(
    Mean_PC_L2_error = mean(Mean_PC_L2_error),
    SE_PC_L2_error = sd(Mean_PC_L2_error) / sqrt(n()),
    Mean_projection_MRAE = mean(Projection_MRAE_clean),
    SE_projection_MRAE = sd(Projection_MRAE_clean) / sqrt(n()),
    Mean_raw_MRAE = mean(Raw_MRAE_clean),
    .groups = "drop"
  )

winner_l2 <- summary_table %>%
  group_by(Scenario) %>%
  slice_min(Mean_PC_L2_error, n = 1, with_ties = FALSE) %>%
  ungroup()

winner_projection <- summary_table %>%
  group_by(Scenario) %>%
  slice_min(Mean_projection_MRAE, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\n\n===== Mean results by scenario and method =====\n")
print(summary_table)

cat("\n\n===== Winner by mean PC loading L2 error =====\n")
print(winner_l2)

cat("\n\n===== Winner by clean-signal projection MRAE =====\n")
print(winner_projection)

cat("\n\n===== Score sparsity recovery =====\n")
print(
  zero_table %>%
    group_by(Scenario, Method, PC) %>%
    summarise(
      True_zero_rate = mean(True_zero_rate),
      Estimated_zero_rate = mean(Estimated_zero_rate),
      Score_support_F1 = mean(Score_support_F1),
      .groups = "drop"
    )
)

cat("\n\n===== Selected tuning parameters =====\n")
print(tuning_table)

############################################################
## 11. Main plot: scenario comparison
############################################################

method_order <- c(
  "MFPCA",
  "Smoothed HPCA",
  "Sparse HPCA",
  "Two-way smooth+sparse HPCA"
)

summary_table$Method <- factor(summary_table$Method, levels = method_order)

p_l2 <- ggplot(summary_table,
               aes(x = Method, y = Mean_PC_L2_error, fill = Method)) +
  geom_col(width = 0.75) +
  geom_errorbar(
    aes(
      ymin = Mean_PC_L2_error - SE_PC_L2_error,
      ymax = Mean_PC_L2_error + SE_PC_L2_error
    ),
    width = 0.20
  ) +
  facet_wrap(~ Scenario, scales = "free_y", ncol = 2) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Different scenarios reveal why two-way regularization helps",
    subtitle = "Lower is better. Scenario D is designed to require both smoothness and sparsity.",
    x = NULL,
    y = "Mean functional/hybrid PC loading L2 error"
  )

p_projection <- ggplot(summary_table,
                       aes(x = Method, y = Mean_projection_MRAE, fill = Method)) +
  geom_col(width = 0.75) +
  geom_errorbar(
    aes(
      ymin = Mean_projection_MRAE - SE_projection_MRAE,
      ymax = Mean_projection_MRAE + SE_projection_MRAE
    ),
    width = 0.20
  ) +
  facet_wrap(~ Scenario, scales = "free_y", ncol = 2) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Clean-signal denoising error using the estimated PC subspace",
    subtitle = "Projection MRAE removes the score-shrinkage bias of sparse methods.",
    x = NULL,
    y = "Projection MRAE against clean centered signal"
  )

print(p_l2 / p_projection)

############################################################
## 12. Figure-6.3-style plot: true black vs estimated red
############################################################

make_curve_df <- function(V, dat, Method, Source) {
  block_sizes <- dat$block_sizes
  bl <- split_matrix_blocks(t(V), block_sizes)

  out <- list()

  for (m in seq_len(ncol(V))) {
    out[[length(out) + 1]] <- data.frame(
      Method = Method,
      Source = Source,
      PC = paste0("PC", m),
      Variable = "FD1: paper variable 1",
      Domain = dat$t1,
      Loading = bl$FD1[m, ]
    )

    out[[length(out) + 1]] <- data.frame(
      Method = Method,
      Source = Source,
      PC = paste0("PC", m),
      Variable = "FD2: paper variable 2",
      Domain = dat$t2,
      Loading = bl$FD2[m, ]
    )

    out[[length(out) + 1]] <- data.frame(
      Method = Method,
      Source = Source,
      PC = paste0("PC", m),
      Variable = "FD3: added hybrid functional",
      Domain = dat$t3,
      Loading = bl$FD3[m, ]
    )
  }

  bind_rows(out)
}

dat <- saved_plot_object
ev <- saved_eval

curve_true_all <- bind_rows(lapply(names(saved_fits), function(nm) {
  make_curve_df(
    V = dat$Vtrue,
    dat = dat,
    Method = nm,
    Source = "True"
  )
}))

curve_est_all <- bind_rows(lapply(names(saved_fits), function(nm) {
  make_curve_df(
    V = ev$aligned_v[[nm]],
    dat = dat,
    Method = nm,
    Source = "Estimated"
  )
}))

curve_plot_df <- bind_rows(curve_true_all, curve_est_all) %>%
  group_by(Method, PC, Variable) %>%
  mutate(Loading_scaled = Loading / max(abs(Loading), na.rm = TRUE)) %>%
  ungroup()

curve_plot_df$Method <- factor(curve_plot_df$Method, levels = method_order)

p_fig63 <- ggplot(
  curve_plot_df,
  aes(x = Domain, y = Loading_scaled, color = Source, linetype = Source)
) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "red")) +
  scale_linetype_manual(values = c("True" = "solid", "Estimated" = "solid")) +
  facet_grid(Variable + Method ~ PC, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    strip.text.y = element_text(size = 7)
  ) +
  labs(
    title = paste0("Figure-6.3-style hybrid simulation: ", PLOT_SCENARIO),
    subtitle = "True PCs black; estimated PCs red",
    x = "Domain",
    y = "Panel-scaled PC loading"
  )

print(p_fig63)

############################################################
## 13. Regular-data loading plot: true black vs estimated red
############################################################

make_regular_df <- function(V, dat, Method, Source) {
  idx <- block_indices(dat$block_sizes)$RD
  RD <- V[idx, , drop = FALSE]

  bind_rows(lapply(seq_len(ncol(RD)), function(m) {
    data.frame(
      Method = Method,
      Source = Source,
      PC = paste0("PC", m),
      Feature = paste0("R", seq_len(nrow(RD))),
      Feature_id = seq_len(nrow(RD)),
      Loading = RD[, m]
    )
  }))
}

rd_true_all <- bind_rows(lapply(names(saved_fits), function(nm) {
  make_regular_df(dat$Vtrue, dat, nm, "True")
}))

rd_est_all <- bind_rows(lapply(names(saved_fits), function(nm) {
  make_regular_df(ev$aligned_v[[nm]], dat, nm, "Estimated")
}))

rd_plot_df <- bind_rows(rd_true_all, rd_est_all)
rd_plot_df$Method <- factor(rd_plot_df$Method, levels = method_order)

p_rd <- ggplot(
  rd_plot_df,
  aes(x = Feature_id, y = Loading, color = Source, group = Source)
) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.4) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "red")) +
  facet_grid(Method ~ PC, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = seq_len(dat$block_sizes["RD"]),
    labels = paste0("R", seq_len(dat$block_sizes["RD"]))
  ) +
  labs(
    title = paste0("Hybrid regular-data loadings: ", PLOT_SCENARIO),
    subtitle = "True black; estimated red",
    x = "Regular feature",
    y = "Loading"
  )

print(p_rd)

############################################################
## 14. Image loading plot for the two-way model
############################################################

make_image_df <- function(v, dat, Method, Source) {
  idx <- block_indices(dat$block_sizes)$IMG
  out <- list()

  for (m in seq_len(ncol(v))) {
    z <- v[idx, m]
    mat <- matrix(z, nrow = dat$img_nrow, ncol = dat$img_ncol, byrow = TRUE)

    out[[m]] <- data.frame(
      Method = Method,
      Source = Source,
      PC = paste0("PC", m),
      row = rep(seq_len(dat$img_nrow), times = dat$img_ncol),
      col = rep(seq_len(dat$img_ncol), each = dat$img_nrow),
      Loading = as.vector(mat)
    )
  }

  bind_rows(out)
}

img_plot_df <- bind_rows(
  make_image_df(
    v = dat$Vtrue,
    dat = dat,
    Method = "Two-way smooth+sparse HPCA",
    Source = "True"
  ),
  make_image_df(
    v = ev$aligned_v[["Two-way smooth+sparse HPCA"]],
    dat = dat,
    Method = "Two-way smooth+sparse HPCA",
    Source = "Estimated"
  )
)

p_img <- ggplot(img_plot_df, aes(x = col, y = row, fill = Loading)) +
  geom_raster() +
  scale_y_reverse() +
  scale_fill_gradient2() +
  coord_equal() +
  facet_grid(Source ~ PC) +
  theme_bw(base_size = 10) +
  labs(
    title = paste0("Hybrid image PC loadings for the two-way smooth+sparse HPCA fit"),
    x = "Image column",
    y = "Image row",
    fill = "Loading"
  )

print(p_img)

############################################################
## 15. Score sparsity recovery plot
############################################################

zero_summary <- zero_table %>%
  group_by(Scenario, Method, PC) %>%
  summarise(
    True_zero_rate = mean(True_zero_rate),
    Estimated_zero_rate = mean(Estimated_zero_rate),
    Score_support_F1 = mean(Score_support_F1),
    .groups = "drop"
  )

zero_long <- zero_summary %>%
  pivot_longer(
    cols = c(True_zero_rate, Estimated_zero_rate),
    names_to = "Rate_type",
    values_to = "Rate"
  )

zero_long$Method <- factor(zero_long$Method, levels = method_order)

p_zero <- ggplot(zero_long, aes(x = PC, y = Rate, fill = Rate_type)) +
  geom_col(position = "dodge") +
  facet_grid(Scenario ~ Method) +
  ylim(0, 1) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0)
  ) +
  labs(
    title = "True vs estimated zero-score rates",
    x = NULL,
    y = "Zero-score rate"
  )

p_f1 <- ggplot(zero_summary, aes(x = Method, y = Score_support_F1, fill = Method)) +
  geom_col(width = 0.75) +
  facet_wrap(~ Scenario, ncol = 2) +
  coord_flip() +
  ylim(0, 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "Score support recovery",
    x = NULL,
    y = "F1 score"
  )

print(p_zero)
print(p_f1)

############################################################
## 16. Optional package plots for the key two-way model
############################################################

plot_pc_functions(saved_fits[["Two-way smooth+sparse HPCA"]])
plot_pc_scores(saved_fits[["Two-way smooth+sparse HPCA"]])

cat("\n\nDone.\n")
