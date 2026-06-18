###############################################################################
# UCI HAR ReMPCA APPLICATION -- rewritten corrected version
#
# Goal:
#   Compare no-penalty ReMPCA/MFPCA with a stable two-way ReMPCA fit for the
#   UCI Human Activity Recognition Using Smartphones dataset.
#
# Recommended penalty region used here:
#   - Row/window sparsity:        gamma_u      = 0
#   - Functional smoothness:      alpha_v      around 2^-6 = 0.015625
#   - Functional loading sparsity: gamma_signal around 2--4
#   - Scalar loading sparsity:     gamma_scalar around 1
#
# Why these settings?
#   The activity windows are not a smooth row-time process, so row smoothing is
#   disabled. Row sparsity is also disabled because aggressive row sparsity can
#   collapse the PC scores to one or two selected windows. The useful penalty is
#   mostly on the loading side: moderate second-order smoothness for the 128-point
#   sensor curves and only light sparsity for functional/scalar loadings.
###############################################################################

###############################################################################
# 0. Packages
###############################################################################

library(ReMPCA)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

###############################################################################
# 1. Output folders
###############################################################################

OUT <- "uci_har_rempca_best_penalty_outputs"

dir.create(OUT, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)
dir.create(file.path(OUT, "tables"), showWarnings = FALSE)
dir.create(file.path(OUT, "diagnostics"), showWarnings = FALSE)

###############################################################################
# 2. Data location
###############################################################################

RESULTS_DIR <- "/users/personnel/pmobina/ReMPCA/Regularized-Two-way-MFPCA/ReMPCA/Results"
base_dir <- file.path(RESULTS_DIR, "UCI HAR Dataset")

if (!dir.exists(base_dir)) {
  stop(
    "Could not find the folder 'UCI HAR Dataset' here:\n",
    base_dir,
    "\n\nPlease check the folder name and path."
  )
}

train_signal_dir <- file.path(base_dir, "train", "Inertial Signals")
test_signal_dir  <- file.path(base_dir, "test", "Inertial Signals")

if (!dir.exists(train_signal_dir)) {
  stop("Could not find train/Inertial Signals folder inside UCI HAR Dataset.")
}

if (!dir.exists(test_signal_dir)) {
  stop("Could not find test/Inertial Signals folder inside UCI HAR Dataset.")
}

cat("\nUCI HAR data found successfully:\n")
cat(base_dir, "\n\n")

###############################################################################
# 3. Read raw inertial signals
###############################################################################

activity_labels <- tibble(
  y = 1:6,
  activity = c(
    "WALKING",
    "WALKING_UPSTAIRS",
    "WALKING_DOWNSTAIRS",
    "SITTING",
    "STANDING",
    "LAYING"
  )
)

# Functional blocks used in this hybrid example.
signal_names <- c(
  "body_acc_x",
  "body_acc_y",
  "body_gyro_x"
)

read_signal <- function(split, signal_name) {

  f <- file.path(
    base_dir,
    split,
    "Inertial Signals",
    paste0(signal_name, "_", split, ".txt")
  )

  if (!file.exists(f)) {
    stop("Could not find signal file:\n", f)
  }

  as.matrix(fread(f, header = FALSE))
}

read_split <- function(split) {

  y_file <- file.path(base_dir, split, paste0("y_", split, ".txt"))
  subject_file <- file.path(base_dir, split, paste0("subject_", split, ".txt"))

  if (!file.exists(y_file)) {
    stop("Could not find label file:\n", y_file)
  }

  if (!file.exists(subject_file)) {
    stop("Could not find subject file:\n", subject_file)
  }

  y <- fread(y_file, header = FALSE)[[1]]
  subject <- fread(subject_file, header = FALSE)[[1]]

  meta <- tibble(
    split = split,
    subject = subject,
    y = y
  ) %>%
    left_join(activity_labels, by = "y")

  signals <- lapply(signal_names, function(s) read_signal(split, s))
  names(signals) <- signal_names

  list(
    meta = meta,
    signals = signals
  )
}

train <- read_split("train")
test  <- read_split("test")

meta_all <- bind_rows(train$meta, test$meta) %>%
  mutate(row_id = row_number())

signals_all <- lapply(signal_names, function(s) {
  rbind(train$signals[[s]], test$signals[[s]])
})
names(signals_all) <- signal_names

cat("Full data size:\n")
cat("Number of windows:", nrow(meta_all), "\n")
cat("Number of time points per signal:", ncol(signals_all[[1]]), "\n\n")

###############################################################################
# 4. Select balanced subset
###############################################################################

ACTIVITIES_KEEP <- c(
  "WALKING",
  "WALKING_UPSTAIRS",
  "SITTING",
  "STANDING"
)

# Increase this after the script runs smoothly.
N_PER_ACTIVITY <- 40

set.seed(2026)

selected_meta <- meta_all %>%
  filter(activity %in% ACTIVITIES_KEEP) %>%
  group_by(activity) %>%
  group_modify(~ slice_sample(.x, n = min(nrow(.x), N_PER_ACTIVITY))) %>%
  ungroup() %>%
  arrange(activity, subject, row_id)

idx <- selected_meta$row_id

meta <- selected_meta %>%
  mutate(
    obs_id = row_number(),
    activity = factor(activity, levels = ACTIVITIES_KEEP)
  )

signals_raw <- lapply(signals_all, function(M) {
  M[idx, , drop = FALSE]
})

n_obs <- nrow(meta)
n_time <- ncol(signals_raw[[1]])
time_grid <- seq(0, 1, length.out = n_time)

cat("Subset used for ReMPCA:\n")
cat("Number of rows/windows:", n_obs, "\n")
print(table(meta$activity))
cat("\n")

###############################################################################
# 5. Create non-functional scalar summaries
###############################################################################

make_scalar_summaries <- function(M, prefix) {

  M <- as.matrix(M)

  tibble(
    !!paste0(prefix, "_mean") := rowMeans(M, na.rm = TRUE),
    !!paste0(prefix, "_sd")   := apply(M, 1, sd, na.rm = TRUE),
    !!paste0(prefix, "_rms")  := sqrt(rowMeans(M^2, na.rm = TRUE))
  )
}

scalar_df <- bind_cols(
  lapply(signal_names, function(s) {
    make_scalar_summaries(signals_raw[[s]], s)
  })
)

###############################################################################
# 6. Standardize functional and scalar blocks
###############################################################################

scale_matrix <- function(M) {

  M <- as.matrix(M)

  center <- colMeans(M, na.rm = TRUE)
  scalev <- apply(M, 2, sd, na.rm = TRUE)

  scalev[!is.finite(scalev)] <- 1
  scalev[scalev == 0] <- 1

  M_scaled <- sweep(M, 2, center, "-")
  M_scaled <- sweep(M_scaled, 2, scalev, "/")

  attr(M_scaled, "center") <- center
  attr(M_scaled, "scale") <- scalev

  M_scaled
}

signals_scaled <- lapply(signals_raw, scale_matrix)
scalar_scaled <- scale_matrix(as.matrix(scalar_df))
colnames(scalar_scaled) <- colnames(scalar_df)

X_scaled_full <- do.call(cbind, c(signals_scaled, list(scalar_scaled)))

###############################################################################
# 7. Stable tuning grids
###############################################################################

NUM_PCS_FIT  <- 2
NUM_PCS_SHOW <- 2

# IMPORTANT:
# Do not use the aggressive original row sparsity grid 0, 32, 64, ...
# It can collapse the row scores to only one or two nonzero windows.
# Also avoid length-one alpha_v grids with tuning_iter > 1; some ReMPCA
# versions throw: incorrect number of subscripts on matrix.

# Target values that usually work well for this subset:
#   alpha_v about 2^-6, gamma_signal about 2--4, gamma_scalar about 1.
# The grid below is narrow and safe, so ReMPCA can still choose by GCV/CV.
alpha_v_grid <- c(2^-8, 2^-7, 2^-6, 2^-5)

gamma_u_grid <- 0L

gamma_signal_grid <- c(2L, 4L, 6L)

gamma_scalar_grid <- c(1L, 2L)

# Values to report as the manual recommendation if you want one fixed setting.
alpha_v_recommended      <- 2^-6
gamma_u_recommended      <- 0L
gamma_signal_recommended <- 4L
gamma_scalar_recommended <- 1L

cat("Stable tuning grids used:\n")
cat("alpha_v_grid:\n")
print(alpha_v_grid)
cat("gamma_u_grid:\n")
print(gamma_u_grid)
cat("gamma_signal_grid:\n")
print(gamma_signal_grid)
cat("gamma_scalar_grid:\n")
print(gamma_scalar_grid)
cat("\n")

###############################################################################
# 8. Build hybrid ReMPCA objects
###############################################################################

make_hd_object <- function(
    alpha_v,
    gamma_u,
    gamma_signal,
    gamma_scalar
) {

  fd_list <- lapply(signal_names, function(s) {
    fdClass(
      data = as.matrix(signals_scaled[[s]]),
      argval = time_grid,
      Smoothing_parameter = alpha_v,
      Sparsity_parameter = gamma_signal
    )
  })
  names(fd_list) <- signal_names

  scalar_rd <- rdClass(
    data = as.matrix(scalar_scaled),
    Sparsity_parameter = gamma_scalar
  )

  hdClass(
    hdlist = c(fd_list, list(scalar = scalar_rd)),
    argval = seq(0, 1, length.out = n_obs),

    # Rows are windows grouped by activity, not observations on a smooth time axis.
    Smoothing_parameter = 0,
    Sparsity_parameter = gamma_u
  )
}

###############################################################################
# 9. Output extraction helpers
###############################################################################

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

as_pc_matrix <- function(x, num_pcs) {

  if (is.null(dim(x))) {
    matrix(as.vector(x), ncol = 1)
  } else {
    x <- as.matrix(x)
    x[, seq_len(min(num_pcs, ncol(x))), drop = FALSE]
  }
}

extract_scores <- function(fit, num_pcs) {

  scores <- fit$PCScores %||%
    fit$PC_Scores %||%
    fit$PC_scores %||%
    fit$scores

  if (is.null(scores)) {
    stop("Could not find PC scores in the ReMPCA output.")
  }

  as_pc_matrix(scores, num_pcs)
}

extract_functions <- function(fit, num_pcs, nblocks) {

  funcs <- fit$PCFunctions %||%
    fit$PC_functions %||%
    fit$PC_functions_list %||%
    fit$loadings

  if (is.null(funcs)) {
    stop("Could not find PC functions/loadings in the ReMPCA output.")
  }

  # Format 1: list of PCs, each containing a list of blocks.
  if (length(funcs) >= num_pcs &&
      is.list(funcs[[1]]) &&
      length(funcs[[1]]) == nblocks) {

    block_list <- vector("list", nblocks)

    for (b in seq_len(nblocks)) {
      block_list[[b]] <- do.call(
        cbind,
        lapply(seq_len(num_pcs), function(pc) as.vector(funcs[[pc]][[b]]))
      )
    }

    return(block_list)
  }

  # Format 2: list of blocks.
  if (length(funcs) == nblocks && !is.list(funcs[[1]][[1]])) {
    return(lapply(funcs, as_pc_matrix, num_pcs = num_pcs))
  }

  stop("Unrecognized PC function/loading format.")
}

make_fit_object <- function(fit, method_name) {

  nblocks <- length(signal_names) + 1

  U <- extract_scores(fit, NUM_PCS_FIT)
  V_blocks <- extract_functions(fit, NUM_PCS_FIT, nblocks = nblocks)
  V <- do.call(rbind, V_blocks)

  colnames(U) <- paste0("PC", seq_len(ncol(U)))
  colnames(V) <- paste0("PC", seq_len(ncol(V)))

  list(
    method = method_name,
    fit = fit,
    U = U,
    V_blocks = V_blocks,
    V = V
  )
}

# Align signs to the no-penalty fit. PC signs are arbitrary, so this only makes
# visual comparison easier.
align_to_reference <- function(fit_obj, ref_obj) {

  for (pc in seq_len(min(ncol(fit_obj$V), ncol(ref_obj$V)))) {
    s <- suppressWarnings(cor(fit_obj$V[, pc], ref_obj$V[, pc], use = "complete.obs"))

    if (is.finite(s) && s < 0) {
      fit_obj$U[, pc] <- -fit_obj$U[, pc]
      fit_obj$V[, pc] <- -fit_obj$V[, pc]

      for (b in seq_along(fit_obj$V_blocks)) {
        fit_obj$V_blocks[[b]][, pc] <- -fit_obj$V_blocks[[b]][, pc]
      }
    }
  }

  fit_obj
}

###############################################################################
# 10. Fit no-penalty ReMPCA / MFPCA baseline
###############################################################################

cat("\nFitting no-penalty ReMPCA/MFPCA baseline...\n")

hd_baseline <- make_hd_object(
  alpha_v = 0,
  gamma_u = 0,
  gamma_signal = 0,
  gamma_scalar = 0
)

fit_baseline_raw <- ReMPCA(
  hd = hd_baseline,
  centerhds = FALSE,
  num_pcs = NUM_PCS_FIT,
  nfolds_u = 3,
  nfolds_v = NULL,
  thresh = 1e-8,
  maxit = 300,
  tuning_iter = 1,
  parallel = FALSE,
  weights = 0,
  smoothness_type = "Second_order",
  sparse_tuning_type = "hard",
  tuning_order = "Sparsity",
  cv.pick = "min",
  sparse_tuning_u = 0,
  sparse_tuning_v = list(0, 0, 0, 0),
  smooth_tuning_u = 0,
  smooth_tuning_v = list(0, 0, 0, 0)
)

fit_baseline <- make_fit_object(
  fit = fit_baseline_raw,
  method_name = "No-penalty ReMPCA"
)

###############################################################################
# 11. Fit recommended penalized ReMPCA
###############################################################################

cat("\nFitting smooth + light-sparse ReMPCA with stable grids...\n")

hd_tuned <- make_hd_object(
  alpha_v = alpha_v_grid,
  gamma_u = gamma_u_grid,
  gamma_signal = gamma_signal_grid,
  gamma_scalar = gamma_scalar_grid
)

fit_tuned_raw <- ReMPCA(
  hd = hd_tuned,
  centerhds = FALSE,
  num_pcs = NUM_PCS_FIT,
  nfolds_u = 3,
  nfolds_v = NULL,
  thresh = 1e-8,
  maxit = 300,
  tuning_iter = 1,
  parallel = FALSE,
  weights = 0,
  smoothness_type = "Second_order",
  sparse_tuning_type = "hard",
  tuning_order = "Sparsity",
  cv.pick = "min",

  # Row sparsity: fixed at zero to avoid score collapse.
  sparse_tuning_u = gamma_u_grid,

  # Loading sparsity: light values only, chosen by CV.
  sparse_tuning_v = list(
    gamma_signal_grid,
    gamma_signal_grid,
    gamma_signal_grid,
    gamma_scalar_grid
  ),

  # Row smoothing: fixed at zero.
  smooth_tuning_u = 0,

  # Functional smoothness only, chosen by GCV; scalar block is not smoothed.
  smooth_tuning_v = list(
    alpha_v_grid,
    alpha_v_grid,
    alpha_v_grid,
    0
  )
)

fit_tuned <- make_fit_object(
  fit = fit_tuned_raw,
  method_name = "Smooth light-sparse ReMPCA"
)

fit_tuned <- align_to_reference(fit_tuned, fit_baseline)

models <- list(
  Baseline = fit_baseline,
  Penalized = fit_tuned
)

###############################################################################
# 12. Save selected/fixed tuning parameters
###############################################################################

collapse_tuning <- function(x) {
  if (is.null(x)) return(NA_character_)
  paste(unlist(x), collapse = ";")
}

selected_tuning <- tibble(
  Quantity = c(
    "alpha_v_recommended_manual",
    "gamma_u_recommended_manual",
    "gamma_signal_recommended_manual",
    "gamma_scalar_recommended_manual",
    "OptimalAlphaU_from_fit",
    "OptimalAlphaV_from_fit",
    "OptimalGammaU_from_fit",
    "OptimalGammaV_from_fit"
  ),
  Value = c(
    as.character(alpha_v_recommended),
    as.character(gamma_u_recommended),
    as.character(gamma_signal_recommended),
    as.character(gamma_scalar_recommended),
    collapse_tuning(fit_tuned_raw$OptimalAlphaU),
    collapse_tuning(fit_tuned_raw$OptimalAlphaV),
    collapse_tuning(fit_tuned_raw$OptimalGammaU),
    collapse_tuning(fit_tuned_raw$OptimalGammaV)
  )
)

print(selected_tuning)

write.csv(
  selected_tuning,
  file.path(OUT, "tables", "selected_penalty_parameters.csv"),
  row.names = FALSE
)

###############################################################################
# 13. Save CV/GCV diagnostic plots when available
###############################################################################

save_rempca_plot <- function(fun_name, fit, filename, width = 7, height = 5) {

  pdf(filename, width = width, height = height)

  if (exists(fun_name, mode = "function")) {

    p <- try(get(fun_name)(fit), silent = TRUE)

    if (inherits(p, "try-error")) {
      plot.new()
      title(main = paste("Could not create", fun_name))
      text(0.5, 0.5, as.character(p), cex = 0.8)
    } else if (inherits(p, "ggplot")) {
      print(p)
    } else {
      # Some ReMPCA plotting functions draw directly and return NULL.
      invisible(p)
    }

  } else {
    plot.new()
    title(main = paste(fun_name, "not found"))
  }

  dev.off()
}

save_rempca_plot(
  "plot_cv_u",
  fit_tuned_raw,
  file.path(OUT, "diagnostics", "cv_u_row_sparsity.pdf")
)

save_rempca_plot(
  "plot_cv_v",
  fit_tuned_raw,
  file.path(OUT, "diagnostics", "cv_v_feature_sparsity.pdf")
)

save_rempca_plot(
  "plot_gcv_u",
  fit_tuned_raw,
  file.path(OUT, "diagnostics", "gcv_u_row_smoothness.pdf")
)

save_rempca_plot(
  "plot_gcv_v",
  fit_tuned_raw,
  file.path(OUT, "diagnostics", "gcv_v_functional_smoothness.pdf")
)

###############################################################################
# 14. Reconstruction, roughness, sparsity, and score diagnostics
###############################################################################

project_reconstruct <- function(X, V, num_pcs = NUM_PCS_FIT) {

  V_use <- as.matrix(V[, seq_len(num_pcs), drop = FALSE])
  G <- crossprod(V_use)
  G_inv <- solve(G + diag(1e-8, ncol(G)))

  Scores <- X %*% V_use %*% G_inv
  X_hat <- Scores %*% t(V_use)

  list(
    Scores = Scores,
    X_hat = X_hat
  )
}

rough_vec <- function(x) {

  x <- as.vector(x)

  if (length(x) < 4) {
    return(0)
  }

  sum(diff(x, differences = 2)^2, na.rm = TRUE) /
    (sum(x^2, na.rm = TRUE) + 1e-8)
}

rough_model <- function(fit_obj, num_pcs = NUM_PCS_SHOW) {

  r <- 0

  for (pc in seq_len(num_pcs)) {

    # Only report roughness of the functional loading curves, not rows.
    for (b in seq_along(signal_names)) {
      r <- r + rough_vec(fit_obj$V_blocks[[b]][, pc])
    }
  }

  r
}

nonzero_count <- function(x, tol = 1e-8) {
  sum(abs(as.vector(x)) > tol, na.rm = TRUE)
}

score_activity_r2 <- function(U, activity, pc) {
  df <- tibble(score = as.vector(U[, pc]), activity = activity)
  fit <- lm(score ~ activity, data = df)
  summary(fit)$r.squared
}

quality_tbl <- bind_rows(
  lapply(names(models), function(m) {

    fit_obj <- models[[m]]
    rec <- project_reconstruct(X_scaled_full, fit_obj$V, NUM_PCS_FIT)

    tibble(
      Method = fit_obj$method,
      Relative_reconstruction_error =
        sum((X_scaled_full - rec$X_hat)^2) / sum(X_scaled_full^2),
      Functional_loading_roughness = rough_model(fit_obj, NUM_PCS_SHOW),
      Nonzero_U_PC1 = nonzero_count(fit_obj$U[, 1]),
      Nonzero_U_PC2 = nonzero_count(fit_obj$U[, 2]),
      Nonzero_V_PC1 = nonzero_count(fit_obj$V[, 1]),
      Nonzero_V_PC2 = nonzero_count(fit_obj$V[, 2]),
      Activity_R2_PC1 = score_activity_r2(fit_obj$U, meta$activity, 1),
      Activity_R2_PC2 = score_activity_r2(fit_obj$U, meta$activity, 2)
    )
  })
)

print(quality_tbl)

write.csv(
  quality_tbl,
  file.path(OUT, "tables", "model_quality_summary.csv"),
  row.names = FALSE
)

###############################################################################
# 15. PC score plots
###############################################################################

score_df <- bind_rows(
  lapply(names(models), function(m) {

    fit_obj <- models[[m]]

    as_tibble(fit_obj$U[, seq_len(NUM_PCS_SHOW), drop = FALSE]) %>%
      setNames(paste0("PC", seq_len(NUM_PCS_SHOW))) %>%
      mutate(
        obs_id = meta$obs_id,
        subject = meta$subject,
        activity = meta$activity,
        Method = fit_obj$method
      )
  })
)

p_scores <- score_df %>%
  ggplot(aes(x = PC1, y = PC2, color = activity, shape = activity)) +
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~ Method, nrow = 1, scales = "free") +
  labs(
    title = "PC scores from smartphone activity data",
    subtitle = "No-penalty ReMPCA versus smooth light-sparse ReMPCA",
    x = "PC1 score",
    y = "PC2 score",
    color = "Activity",
    shape = "Activity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p_scores)

ggsave(
  file.path(OUT, "figures", "01_pc_scores_activity.pdf"),
  p_scores,
  width = 10,
  height = 5
)

###############################################################################
# 16. Row score plot
###############################################################################

row_order_tbl <- meta %>%
  arrange(activity, subject, obs_id) %>%
  mutate(order_index = row_number()) %>%
  select(obs_id, order_index, activity)

row_score_long <- score_df %>%
  left_join(row_order_tbl, by = c("obs_id", "activity")) %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "score"
  )

p_row_scores <- row_score_long %>%
  ggplot(aes(x = order_index, y = score, color = activity)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(alpha = 0.7, linewidth = 0.7) +
  facet_grid(PC ~ Method, scales = "free_y") +
  labs(
    title = "Row scores across activity windows",
    subtitle = "Recommended fit uses no row sparsity, avoiding one-window score collapse",
    x = "Windows grouped by activity",
    y = "Row score",
    color = "Activity"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p_row_scores)

ggsave(
  file.path(OUT, "figures", "02_row_scores_no_collapse.pdf"),
  p_row_scores,
  width = 10,
  height = 7
)

###############################################################################
# 17. Functional loading comparison
###############################################################################

pc_labs <- paste0("PC", seq_len(NUM_PCS_SHOW))

loading_df <- bind_rows(
  lapply(names(models), function(m) {

    fit_obj <- models[[m]]

    bind_rows(
      lapply(seq_along(signal_names), function(b) {

        as_tibble(fit_obj$V_blocks[[b]][, seq_len(NUM_PCS_SHOW), drop = FALSE]) %>%
          setNames(pc_labs) %>%
          mutate(
            time = time_grid,
            signal = signal_names[b],
            Method = fit_obj$method
          ) %>%
          pivot_longer(
            cols = starts_with("PC"),
            names_to = "PC",
            values_to = "loading"
          )
      })
    )
  })
)

p_loadings <- loading_df %>%
  ggplot(aes(x = time, y = loading, color = Method, linetype = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_line(linewidth = 0.95) +
  facet_grid(signal ~ PC, scales = "free_y") +
  labs(
    title = "Functional loadings for inertial signals",
    subtitle = "Recommended ReMPCA: second-order smoothness plus light loading sparsity",
    x = "Normalized time within window",
    y = "Loading",
    color = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p_loadings)

ggsave(
  file.path(OUT, "figures", "03_functional_loading_comparison.pdf"),
  p_loadings,
  width = 10,
  height = 7
)

###############################################################################
# 18. Non-functional scalar loading plot
###############################################################################

scalar_block_index <- length(signal_names) + 1

scalar_loading_df <- bind_rows(
  lapply(names(models), function(m) {

    fit_obj <- models[[m]]

    as_tibble(fit_obj$V_blocks[[scalar_block_index]][, seq_len(NUM_PCS_SHOW), drop = FALSE]) %>%
      setNames(pc_labs) %>%
      mutate(
        covariate = colnames(scalar_scaled),
        Method = fit_obj$method
      ) %>%
      pivot_longer(
        cols = starts_with("PC"),
        names_to = "PC",
        values_to = "loading"
      )
  })
)

p_scalar <- scalar_loading_df %>%
  ggplot(aes(x = covariate, y = loading, fill = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_x") +
  labs(
    title = "Non-functional scalar loadings",
    subtitle = "Scalar block contains mean, standard deviation, and RMS summaries",
    x = NULL,
    y = "Loading",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p_scalar)

ggsave(
  file.path(OUT, "figures", "04_scalar_loading_comparison.pdf"),
  p_scalar,
  width = 9,
  height = 6
)

###############################################################################
# 19. Save LaTeX table for the paper
###############################################################################

latex_table <- "
\\begin{table}[!htbp]
\\centering
\\caption{Hybrid structure of the smartphone activity data.}
\\label{tab:har_hybrid_structure}
\\begin{tabular}{llll}
\\toprule
Data component & Meaning & Type & Model role \\\\
\\midrule
Activity window & One sensor window & Row direction & $\\bm u$ \\\\
Body acceleration, x-axis & Curve over 128 time points & Functional & $\\bm v_1$ \\\\
Body acceleration, y-axis & Curve over 128 time points & Functional & $\\bm v_2$ \\\\
Body gyroscope, x-axis & Curve over 128 time points & Functional & $\\bm v_3$ \\\\
Signal summaries & Mean, standard deviation, and RMS & Non-functional & $\\bm w_1$ \\\\
Activity label & Walking, upstairs walking, sitting, standing & External grouping & Used for plots \\\\
\\bottomrule
\\end{tabular}
\\end{table}
"

writeLines(
  latex_table,
  con = file.path(OUT, "tables", "har_hybrid_structure_table.tex")
)

###############################################################################
# 20. Done
###############################################################################

cat("\n====================================================\n")
cat("Done. Outputs saved in:", OUT, "\n\n")

cat("Manual recommended penalties:\n")
cat("alpha_v      =", alpha_v_recommended, "\n")
cat("gamma_u      =", gamma_u_recommended, "\n")
cat("gamma_signal =", gamma_signal_recommended, "\n")
cat("gamma_scalar =", gamma_scalar_recommended, "\n\n")
cat("Selected by ReMPCA from stable grids:\n")
print(selected_tuning)
cat("\n")

cat("Main figures:\n")
cat("figures/01_pc_scores_activity.pdf\n")
cat("figures/02_row_scores_no_collapse.pdf\n")
cat("figures/03_functional_loading_comparison.pdf\n")
cat("figures/04_scalar_loading_comparison.pdf\n\n")

cat("Tables:\n")
cat("tables/selected_penalty_parameters.csv\n")
cat("tables/model_quality_summary.csv\n")
cat("tables/har_hybrid_structure_table.tex\n")
cat("====================================================\n")
