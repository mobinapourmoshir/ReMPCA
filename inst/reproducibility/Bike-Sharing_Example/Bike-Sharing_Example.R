###############################################################################
# UCI Bike Sharing Dataset
# No-penalty baseline versus PC-specific two-way smooth/sparse ReMPCA
#
#   1) loads bike_day.rda and bike_hour.rda directly from data/
#   2) prepares the hybrid data
#   3) fits the two-PC no-penalty baseline
#   4) fits PC1 on the original data with penalties
#   5) deflates the data and fits PC2 with penalties
#   6) reproduces the comparison plots
###############################################################################

## 0. Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ReMPCA)
  library(tidyverse)
  library(patchwork)
  library(ggh4x)
})

out_dir <- file.path("inst", "reproducibility", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t_start <- Sys.time()

## 1. Load data ---------------------------------------------------------------

load(file.path("data", "bike_day.rda"))
load(file.path("data", "bike_hour.rda"))

day_dat  <- bike_day  |> mutate(dteday = as.Date(dteday))
hour_dat <- bike_hour |> mutate(dteday = as.Date(dteday))

## 2. Daily functional curves (24 hourly values per day) ----------------------

# Keep only days with all 24 hourly records.
complete_days <- hour_dat |>
  group_by(dteday) |>
  summarize(n_hours = n_distinct(hr), .groups = "drop") |>
  filter(n_hours == 24) |>
  pull(dteday)

hour_complete <- hour_dat |> filter(dteday %in% complete_days)
day_complete  <- day_dat  |> filter(dteday %in% complete_days) |> arrange(dteday)

# One row per day, one column per hour (h00, ..., h23).
hourly_curves <- function(variable) {
  hour_complete |>
    select(dteday, hr, value = all_of(variable)) |>
    arrange(dteday, hr) |>
    mutate(hr = paste0("h", sprintf("%02d", hr))) |>
    pivot_wider(names_from = hr, values_from = value) |>
    arrange(dteday) |>
    select(-dteday) |>
    as.matrix()
}

X_raw <- list(
  Count       = hourly_curves("cnt"),
  Temperature = hourly_curves("temp"),
  Humidity    = hourly_curves("hum"),
  Wind        = hourly_curves("windspeed")
)

## 3. Non-functional daily covariates -----------------------------------------

# Dummy-coded calendar variables plus the four daily weather averages
# (30 columns). Figure 4 (right panel) displays the nine covariates listed
# in Table 1 of the paper; all 30 loadings are written to the CSV file.
regular_df <- day_complete |>
  transmute(
    season        = factor(season),
    year          = factor(yr),
    month         = factor(mnth),
    holiday       = factor(holiday),
    weekday       = factor(weekday),
    workingday    = factor(workingday),
    weathersit    = factor(weathersit),
    temp_day      = temp,
    atemp_day     = atemp,
    hum_day       = hum,
    windspeed_day = windspeed
  )

X_regular_raw <- model.matrix(
  ~ season + year + month + holiday + weekday + workingday + weathersit +
    temp_day + atemp_day + hum_day + windspeed_day - 1,
  data = regular_df
)
regular_feature_names <- colnames(X_regular_raw)

displayed_features <- c(
  "temp_day", "atemp_day", "hum_day", "windspeed_day",
  "season1", "season2", "season3", "season4", "workingday1"
)

## 4. Scaling -----------------------------------------------------------------

# Functional blocks: one global mean and SD per block.
global_scale <- function(X) (X - mean(X)) / sd(as.vector(X))

# Non-functional block: columnwise standardisation.
X_regular <- scale(X_regular_raw, center = TRUE, scale = TRUE)
X_regular[is.na(X_regular)] <- 0

blocks_original <- list(
  Count_curve       = global_scale(X_raw$Count),
  Temperature_curve = global_scale(X_raw$Temperature),
  Humidity_curve    = global_scale(X_raw$Humidity),
  Wind_curve        = global_scale(X_raw$Wind),
  Daily_covariates  = as.matrix(X_regular)
)

n_days    <- nrow(blocks_original$Count_curve)
hour_grid <- seq(0, 23, length.out = 24)
day_grid  <- seq(0, 1, length.out = n_days)

## 5. Tuning parameters -------------------------------------------------------

# alpha = smoothness parameter; gamma = number of entries set to zero
# (hard thresholding). v entries are ordered Count, Temperature, Humidity,
# Wind, Daily covariates; the non-functional block is never smoothed.
tuning <- list(
  HPCA = list(
    alpha_u = 0, gamma_u = 0L,
    alpha_v = c(0, 0, 0, 0, 0), gamma_v = c(0L, 0L, 0L, 0L, 0L)
  ),
  ReHPCA_PC1 = list(
    alpha_u = 2^-6, gamma_u = 1L,
    alpha_v = c(2^-2, 2^-6, 2^-5, 2^-4, 0), gamma_v = c(0L, 2L, 0L, 0L, 8L)
  ),
  ReHPCA_PC2 = list(
    alpha_u = 2^-7, gamma_u = 0L,
    alpha_v = c(2^-3, 2^-4, 2^-4, 2^-4, 0), gamma_v = c(1L, 0L, 10L, 0L, 11L)
  )
)

# hdClass object for the five blocks under one tuning set.
make_hd <- function(blocks, tp) {
  fd_blocks <- lapply(1:4, function(j) {
    fdClass(
      data = blocks[[j]],
      argval = hour_grid,
      Smoothing_parameter = tp$alpha_v[j],
      Sparsity_parameter = tp$gamma_v[j]
    )
  })
  names(fd_blocks) <- names(blocks)[1:4]

  rd_block <- rdClass(data = blocks[[5]], Sparsity_parameter = tp$gamma_v[5])

  hdClass(
    hdlist = c(fd_blocks, list(Daily_covariates = rd_block)),
    argval = day_grid,
    Smoothing_parameter = tp$alpha_u,
    Sparsity_parameter = tp$gamma_u
  )
}

## 6. HPCA: no-penalty baseline with two PCs ----------------------------------

set.seed(2026)
fit_hpca <- ReMPCA(
  hd = make_hd(blocks_original, tuning$HPCA),
  centerhds = FALSE,
  num_pcs = 2,
  nfolds_u = 5,
  nfolds_v = NULL,
  thresh = 1e-6,
  maxit = 500,
  tuning_iter = 1,
  parallel = FALSE,
  weights = NULL,
  smoothness_type = "Second_order",
  sparse_tuning_type = "soft",
  tuning_order = "Sparsity",
  cv.pick = "1se",
  sparse_tuning_u = 0,
  smooth_tuning_u = 0,
  sparse_tuning_v = list(0, 0, 0, 0, 0),
  smooth_tuning_v = list(0, 0, 0, 0, 0)
)

# The sign of a principal component is not identified, and the sign returned
# can depend on the linear-algebra library. Orient each HPCA component so that
# its scores sum to a negative value (the orientation shown in the paper).
for (pc in 1:2) {
  if (sum(fit_hpca$PCScores[, pc]) > 0) {
    fit_hpca$PCScores[, pc] <- -fit_hpca$PCScores[, pc]
    fit_hpca$PCFunctions[[pc]] <- lapply(fit_hpca$PCFunctions[[pc]], `-`)
  }
}

## 7. ReHPCA: PC1 on the original data, PC2 on the deflated data --------------

fit_rehpca_one_pc <- function(blocks, tp) {
  set.seed(2026)
  ReMPCA(
    hd = make_hd(blocks, tp),
    centerhds = FALSE,
    num_pcs = 1,
    nfolds_u = 2,
    nfolds_v = NULL,
    thresh = 1e-4,
    maxit = 140,
    tuning_iter = 1,
    parallel = FALSE,
    weights = NULL,
    smoothness_type = "Second_order",
    sparse_tuning_type = "hard",
    tuning_order = "Sparsity",
    cv.pick = "min",
    sparse_tuning_u = tp$gamma_u,
    smooth_tuning_u = tp$alpha_u,
    sparse_tuning_v = as.list(tp$gamma_v),
    smooth_tuning_v = as.list(tp$alpha_v)
  )
}

fit_pc1 <- fit_rehpca_one_pc(blocks_original, tuning$ReHPCA_PC1)

# Deflation: X_b <- X_b - u1 v1b' for every block b.
u_pc1 <- as.numeric(fit_pc1$PCScores[, 1])
blocks_residual <- Map(
  function(X, v) X - outer(u_pc1, as.numeric(v)),
  blocks_original,
  fit_pc1$PCFunctions[[1]]
)

fit_pc2 <- fit_rehpca_one_pc(blocks_residual, tuning$ReHPCA_PC2)

# Combine both components and align each sign with HPCA.
fit_rehpca <- list(
  PCScores = cbind(PC1 = u_pc1, PC2 = as.numeric(fit_pc2$PCScores[, 1])),
  PCFunctions = list(fit_pc1$PCFunctions[[1]], fit_pc2$PCFunctions[[1]])
)

for (pc in 1:2) {
  if (cor(fit_hpca$PCScores[, pc], fit_rehpca$PCScores[, pc]) < 0) {
    fit_rehpca$PCScores[, pc] <- -fit_rehpca$PCScores[, pc]
    fit_rehpca$PCFunctions[[pc]] <- lapply(fit_rehpca$PCFunctions[[pc]], `-`)
  }
}

fits <- list(HPCA = fit_hpca, ReHPCA = fit_rehpca)

## 8. Tables: tuning parameters and non-functional loadings -------------------

tuning_table <- imap_dfr(tuning, function(tp, fit_name) {
  tibble(
    fit = fit_name,
    direction_or_block = c("u (days)", "Count", "Temperature", "Humidity",
                           "Wind", "Daily covariates"),
    alpha_smoothness = c(tp$alpha_u, tp$alpha_v),
    gamma_sparsity = c(tp$gamma_u, tp$gamma_v),
    n_entries = c(n_days, sapply(blocks_original, ncol))
  )
})
write_csv(tuning_table, file.path(out_dir, "bike_tuning_parameters.csv"))

regular_loadings <- imap_dfr(fits, function(fit, method) {
  map_dfr(1:2, function(pc) {
    tibble(
      method = method,
      PC = paste0("PC", pc),
      feature = regular_feature_names,
      loading = as.numeric(fit$PCFunctions[[pc]][[5]])
    )
  })
}) |>
  mutate(
    is_zero = loading == 0,
    shown_in_figure4 = feature %in% displayed_features
  )
write_csv(regular_loadings, file.path(out_dir, "bike_nonfunctional_loadings.csv"))

zeroed <- regular_loadings |>
  filter(method == "ReHPCA", is_zero) |>
  group_by(PC) |>
  summarize(features = paste(feature, collapse = ", "), .groups = "drop")
cat("Non-functional covariates set to zero by ReHPCA:\n")
cat(sprintf("  %s: %s\n", zeroed$PC, zeroed$features), sep = "")

## 9. Plot data ---------------------------------------------------------------

method_cols <- c("HPCA" = "#D95F02", "ReHPCA" = "#0072B2")
functional_blocks <- c("Count", "Temperature", "Humidity", "Wind")

score_long <- imap_dfr(fits, function(fit, method) {
  tibble(
    date = day_complete$dteday,
    Method = method,
    PC1 = as.numeric(fit$PCScores[, 1]),
    PC2 = as.numeric(fit$PCScores[, 2])
  )
}) |>
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "Score")

functional_long <- imap_dfr(fits, function(fit, method) {
  map_dfr(1:2, function(pc) {
    map_dfr(1:4, function(j) {
      tibble(
        Method = method,
        PC = paste0("PC", pc),
        Functional_Block = functional_blocks[j],
        Hour = hour_grid,
        Loading = as.numeric(fit$PCFunctions[[pc]][[j]])
      )
    })
  })
}) |>
  mutate(Functional_Block = factor(Functional_Block, levels = functional_blocks))

# Displayed covariates, ordered by their largest absolute HPCA loading.
# The day-type row shows two bars: workingday1 and its complement
# (non-working day), whose loading is the negative of workingday1.
feature_order <- regular_loadings |>
  filter(method == "HPCA", shown_in_figure4) |>
  group_by(feature) |>
  summarize(max_abs = max(abs(loading)), .groups = "drop") |>
  arrange(desc(max_abs)) |>
  pull(feature)

regular_plot_df <- regular_loadings |>
  filter(shown_in_figure4) |>
  bind_rows(
    regular_loadings |>
      filter(feature == "workingday1") |>
      mutate(loading = -loading)
  ) |>
  mutate(feature = factor(feature, levels = rev(feature_order)))

## 10. Figure 4: u scores, functional loadings, non-functional loadings -------

compact_theme <- theme_bw(base_size = 15) +
  theme(
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey92", color = "grey45", linewidth = 0.40),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", linewidth = 0.35),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.5, color = "black")
  )

p_u_time <- ggplot(score_long, aes(date, Score, color = Method, group = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.30, color = "grey40") +
  geom_line(linewidth = 0.85, alpha = 0.92) +
  facet_wrap(~ PC, nrow = 1, scales = "fixed") +
  scale_color_manual(values = method_cols, guide = "none") +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(title = "u direction: PC scores over ordered days", x = NULL, y = "u score") +
  compact_theme +
  theme(
    plot.title = element_text(face = "bold", size = 19, margin = margin(b = 4)),
    strip.text.x = element_text(face = "bold", size = 16),
    axis.title.y.left = element_text(face = "bold", size = 16, margin = margin(r = 5)),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    panel.spacing.x = unit(0.05, "cm"),
    plot.margin = margin(2, 2, 1, 0)
  )

p_v_function <- ggplot(functional_long, aes(Hour, Loading, color = Method, group = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.30, color = "grey40") +
  geom_line(linewidth = 0.95, alpha = 0.95) +
  facet_grid(
    Functional_Block ~ PC,
    scales = "free_y",
    switch = "y",
    labeller = labeller(Functional_Block = c(
      Count = "Count", Temperature = "Temp.", Humidity = "Humidity", Wind = "Wind"
    ))
  ) +
  scale_color_manual(values = method_cols, guide = "none") +
  scale_x_continuous(breaks = c(0, 6, 12, 18),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                     labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0.08, 0.08))) +
  labs(title = "v direction: functional PC loading curves",
       x = "Hour of day", y = "PC loading") +
  compact_theme +
  theme(
    plot.title = element_text(face = "bold", size = 19, hjust = 0.5, margin = margin(b = 5)),
    plot.title.position = "panel",
    strip.text.x = element_text(face = "bold", size = 16, margin = margin(3, 3, 3, 3)),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 16,
                                     margin = margin(2, 2, 2, 2)),
    strip.placement = "outside",
    strip.switch.pad.grid = unit(0, "pt"),
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.04, "cm"),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 4)),
    axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 3)),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.margin = margin(2, 2, 0, 2)
  )

p_v_regular <- ggplot(regular_plot_df, aes(loading, feature, fill = method)) +
  geom_vline(xintercept = 0, linewidth = 0.30, color = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.95) +
  facet_grid(. ~ PC, scales = "free_x") +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2),
                     labels = scales::label_number(accuracy = 0.01),
                     guide = guide_axis(check.overlap = TRUE)) +
  labs(title = "v direction: regular covariate loadings", x = "PC loading", y = NULL) +
  compact_theme +
  theme(
    plot.title = element_text(face = "bold", size = 19, hjust = 0.5, margin = margin(b = 5)),
    plot.title.position = "panel",
    strip.text.x = element_text(face = "bold", size = 16, margin = margin(3, 3, 3, 3)),
    panel.spacing.x = unit(0.20, "cm"),
    axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 8)),
    axis.text.x = element_text(size = 11, color = "black", margin = margin(t = 3)),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.ticks.length.x = unit(1.5, "pt"),
    plot.margin = margin(2, 2, 0, 2)
  )

# Method legend drawn as a small stand-alone panel (upper left).
p_legend <- ggplot(tibble(Method = names(method_cols), y = c(0.65, 0.35))) +
  geom_segment(aes(x = 0.04, xend = 0.34, y = y, yend = y, color = Method),
               linewidth = 1.7) +
  geom_text(aes(x = 0.42, y = y, label = Method),
            hjust = 0, size = 5.5, fontface = "bold") +
  scale_color_manual(values = method_cols, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

fig4 <- ((p_legend | p_u_time) + plot_layout(widths = c(0.17, 1))) /
  ((wrap_elements(full = p_v_function) | wrap_elements(full = p_v_regular)) +
     plot_layout(widths = c(1.65, 1.00))) +
  plot_layout(heights = c(0.55, 1.80))

ggsave(file.path(out_dir, "bike_compare_u_and_v_HPCA_ReHPCA1.pdf"), fig4,
       width = 16, height = 10, units = "in")

## 11. Figure 5: fitted daily-hour heatmaps of bike counts ---------------------

# Two-PC reconstruction of the scaled Count block, back-transformed to counts.
count_mean <- mean(X_raw$Count)
count_sd   <- sd(as.vector(X_raw$Count))

fitted_counts <- function(fit) {
  rec <- outer(fit$PCScores[, 1], as.numeric(fit$PCFunctions[[1]][[1]])) +
    outer(fit$PCScores[, 2], as.numeric(fit$PCFunctions[[2]][[1]]))
  pmax(rec * count_sd + count_mean, 0)
}

heatmap_df <- imap_dfr(fits, function(fit, method) {
  M <- fitted_counts(fit)
  tibble(
    date = rep(day_complete$dteday, each = ncol(M)),
    hour = rep(hour_grid, times = nrow(M)),
    log_count = log1p(as.vector(t(M))),
    Method = paste(method, "fitted")
  )
})

fill_upper <- quantile(heatmap_df$log_count, probs = 0.99)

heat_palette <- colorRampPalette(c(
  "#2D004B", "#542788", "#8073AC", "#B358A0",
  "#D65F5F", "#F07C3E", "#F6B44B", "#F4D35E", "#FFF3A3"
))(256)

fig5 <- ggplot(heatmap_df, aes(hour, date, fill = log_count)) +
  geom_tile(width = 1, height = 1) +
  facet_wrap(~ Method, nrow = 1) +
  scale_fill_gradientn(colors = heat_palette, limits = c(0, fill_upper),
                       oob = scales::squish, name = "log(1 + count)") +
  scale_x_continuous(breaks = seq(0, 23, by = 4), expand = c(0, 0)) +
  scale_y_date(date_breaks = "3 months", date_labels = "%Y-%m",
               expand = expansion(mult = c(0.002, 0.002))) +
  labs(title = "Daily-hour fitted heatmaps of bike-rental demand",
       x = "Hour of day", y = "Date") +
  theme_minimal(base_size = 17) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_rect(fill = "grey92", color = "grey55", linewidth = 0.25),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(margin = margin(t = 4)),
    axis.text.y = element_text(margin = margin(r = 4)),
    legend.position = "bottom",
    legend.title = element_text(size = 17),
    legend.text = element_text(size = 16),
    legend.key.width = unit(1.7, "cm"),
    legend.key.height = unit(0.30, "cm"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.25, "cm"),
    plot.margin = margin(5, 8, 5, 8),
    aspect.ratio = 0.55
  )

ggsave(file.path(out_dir, "heatmap_HPCA_vs_ReHPCA_fitted1.pdf"), fig5,
       width = 16, height = 7, units = "in")

## 12. Figure 6: mean curves plus/minus the ReHPCA functional PCs --------------

# Raw-scale mean curve +/- sd(u_k) * v_kb * sd(block b), for k = 1, 2.
score_sd <- apply(fit_rehpca$PCScores, 2, sd)

mean_pm_df <- map_dfr(1:2, function(k) {
  map_dfr(1:4, function(j) {
    mu  <- colMeans(X_raw[[j]])
    dev <- score_sd[k] * as.numeric(fit_rehpca$PCFunctions[[k]][[j]]) *
      sd(as.vector(X_raw[[j]]))
    tibble(
      PC = paste0("PC", k),
      Functional_Block = functional_blocks[j],
      Hour = hour_grid,
      `Mean - PC` = mu - dev,
      Mean = mu,
      `Mean + PC` = mu + dev
    )
  })
}) |>
  pivot_longer(c(`Mean - PC`, Mean, `Mean + PC`),
               names_to = "Curve", values_to = "Value") |>
  mutate(
    Curve = factor(Curve, levels = c("Mean - PC", "Mean", "Mean + PC")),
    Functional_Block = factor(Functional_Block, levels = functional_blocks)
  )

# "-" and "+" markers along the perturbed curves (65 points per curve).
dense_hours <- seq(0, 23, length.out = 65)
pm_symbol_df <- mean_pm_df |>
  filter(Curve != "Mean") |>
  group_by(PC, Functional_Block, Curve) |>
  group_modify(~ tibble(
    Hour = dense_hours,
    Value = approx(.x$Hour, .x$Value, xout = dense_hours, rule = 2)$y,
    Symbol = if (.y$Curve == "Mean - PC") "-" else "+"
  )) |>
  ungroup()

pm_cols <- c("Mean - PC" = "#D62728", "Mean" = "black", "Mean + PC" = "#009E73")

fig6 <- ggplot(mean_pm_df, aes(Hour, Value, color = Curve, group = Curve)) +
  geom_line(linewidth = 1.15, lineend = "round") +
  geom_text(data = pm_symbol_df, aes(label = Symbol), size = 3.1, show.legend = FALSE) +
  facet_grid2(rows = vars(PC), cols = vars(Functional_Block),
              scales = "free_y", independent = "y", switch = "y", axes = "y") +
  scale_color_manual(values = pm_cols) +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 23),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(0.10, 0.10))) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.8, alpha = 1))) +
  labs(title = "Mean curves plus and minus ReHPCA functional PCs",
       x = "Hour of day", y = NULL) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 16)),
    strip.text.x = element_text(face = "bold", size = 18, color = "black",
                                margin = margin(5, 5, 5, 5)),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 20, color = "black",
                                     margin = margin(4, 4, 4, 4)),
    strip.background = element_rect(fill = "grey92", color = "grey45", linewidth = 0.5),
    strip.placement = "outside",
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", linewidth = 0.5),
    panel.spacing.x = unit(0.08, "cm"),
    panel.spacing.y = unit(0.25, "cm"),
    axis.title.x = element_text(size = 19, face = "bold", margin = margin(t = 12)),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    legend.key.width = unit(1.6, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    plot.margin = margin(12, 12, 12, 12)
  )

ggsave(file.path(out_dir, "bike_mean_ReHPCA1.pdf"), fig6,
       width = 17, height = 8.5, units = "in")

## 13. Session information ----------------------------------------------------

cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "bike_sessionInfo.txt"))
