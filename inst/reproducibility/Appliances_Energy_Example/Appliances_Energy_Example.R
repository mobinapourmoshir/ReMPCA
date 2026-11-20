###############################################################################
# UCI Appliances Energy Prediction Dataset
###############################################################################

## 0. Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ReMPCA)
  library(tidyverse)
  library(lubridate)
  library(patchwork)
})

out_dir <- file.path("inst", "reproducibility", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

t_start <- Sys.time()

## 1. Load data ---------------------------------------------------------------

load(file.path("data", "energydata.rda"))

indoor_temp_cols <- paste0("T", c(1, 2, 3, 4, 5, 7, 8, 9))
indoor_hum_cols  <- paste0("RH_", c(1, 2, 3, 4, 5, 7, 8, 9))

appliance_dat <- appliance_dat |>
  mutate(
    date = ymd_hms(date),
    day = as.Date(date),
    slot = hour(date) * 6 + minute(date) / 10,
    hour_of_day = slot / 6,
    Indoor_temperature = rowMeans(as.matrix(across(all_of(indoor_temp_cols))), na.rm = TRUE),
    Indoor_humidity = rowMeans(as.matrix(across(all_of(indoor_hum_cols))), na.rm = TRUE),
    Outdoor_temperature = T_out,
    Outdoor_humidity = RH_out
  )

## 2. Daily functional curves (144 ten-minute values per day) -----------------

# Keep only calendar days with all 144 ten-minute measurements.
complete_days <- appliance_dat |>
  group_by(day) |>
  summarize(n_slots = n_distinct(slot), .groups = "drop") |>
  filter(n_slots == 144) |>
  pull(day)

appliance_complete <- appliance_dat |>
  filter(day %in% complete_days) |>
  arrange(day, slot)

day_complete <- appliance_complete |> distinct(day) |> arrange(day)

n_days    <- nrow(day_complete)
hour_grid <- sort(unique(appliance_complete$hour_of_day))
day_grid  <- seq(0, 1, length.out = n_days)

# One row per day, one column per ten-minute slot.
daily_curves <- function(variable) {
  M <- appliance_complete |>
    select(day, slot, value = all_of(variable)) |>
    arrange(day, slot) |>
    mutate(slot = paste0("s", sprintf("%03d", slot))) |>
    pivot_wider(names_from = slot, values_from = value) |>
    arrange(day) |>
    select(-day) |>
    as.matrix()
  rownames(M) <- as.character(day_complete$day)
  M
}

functional_names <- c("Appliances", "Lights", "Indoor temp.", "Indoor humidity",
                      "Outdoor temp.", "Outdoor humidity")

X_raw <- list(
  daily_curves("Appliances"),
  daily_curves("lights"),
  daily_curves("Indoor_temperature"),
  daily_curves("Indoor_humidity"),
  daily_curves("Outdoor_temperature"),
  daily_curves("Outdoor_humidity")
)
names(X_raw) <- functional_names

## 3. Non-functional daily covariates -----------------------------------------

regular_df <- appliance_complete |>
  group_by(day) |>
  summarize(
    month_num = month(first(day)),
    weekday_num = wday(first(day), week_start = 1),
    Day_type = if_else(weekday_num %in% c(6, 7), "Weekend", "Weekday"),
    Pressure_day = mean(Press_mm_hg, na.rm = TRUE),
    Windspeed_day = mean(Windspeed, na.rm = TRUE),
    Visibility_day = mean(Visibility, na.rm = TRUE),
    Tdewpoint_day = mean(Tdewpoint, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(day) |>
  mutate(month_factor = factor(month_num, levels = sort(unique(month_num))))

# Month indicators for the observed months and both day-type indicators.
X_month <- model.matrix(~ month_factor - 1, data = regular_df)
colnames(X_month) <- paste0("month_", levels(regular_df$month_factor))

X_regular_raw <- cbind(
  as.matrix(regular_df[, c("Pressure_day", "Windspeed_day",
                           "Visibility_day", "Tdewpoint_day")]),
  X_month,
  daytype_weekday = as.integer(regular_df$Day_type == "Weekday"),
  daytype_weekend = as.integer(regular_df$Day_type == "Weekend")
)
regular_feature_names <- colnames(X_regular_raw)

## 4. Scaling -----------------------------------------------------------------

# Functional blocks: one global mean and SD per block.
global_scale <- function(X) (X - mean(X)) / sd(as.vector(X))

# Non-functional block: columnwise standardisation.
X_regular <- scale(X_regular_raw, center = TRUE, scale = TRUE)
X_regular[is.na(X_regular)] <- 0

blocks_original <- c(
  lapply(X_raw, global_scale),
  list(Daily_covariates = as.matrix(X_regular))
)
names(blocks_original) <- c(
  "Appliance_curve", "Lights_curve",
  "Indoor_temperature_curve", "Indoor_humidity_curve",
  "Outdoor_temperature_curve", "Outdoor_humidity_curve",
  "Daily_covariates"
)

## 5. Tuning parameters -------------------------------------------------------

# alpha = smoothness parameter; gamma = number of entries set to zero
# (hard thresholding). v entries follow the block order above; the
# non-functional block is never smoothed.
tuning <- list(
  HPCA = list(
    alpha_u = 0, gamma_u = 0L,
    alpha_v = rep(0, 7), gamma_v = rep(0L, 7)
  ),
  ReHPCA_PC1 = list(
    alpha_u = 2^-12, gamma_u = 45L,
    alpha_v = c(2^0, 2^1, 2^-3, 2^-3, 2^-3, 2^-3, 0),
    gamma_v = c(8L, 50L, 0L, 0L, 0L, 0L, 3L)
  ),
  ReHPCA_PC2 = list(
    alpha_u = 2^-10, gamma_u = 55L,
    alpha_v = c(2^2, 2^2, 2^-2, 2^-2, 2^-2, 2^-2, 0),
    gamma_v = c(20L, 60L, 0L, 0L, 0L, 0L, 5L)
  )
)

# hdClass object for the seven blocks under one tuning set.
make_hd <- function(blocks, tp) {
  fd_blocks <- lapply(1:6, function(j) {
    fdClass(
      data = blocks[[j]],
      argval = hour_grid,
      Smoothing_parameter = tp$alpha_v[j],
      Sparsity_parameter = tp$gamma_v[j]
    )
  })
  names(fd_blocks) <- names(blocks)[1:6]

  rd_block <- rdClass(data = blocks[[7]], Sparsity_parameter = tp$gamma_v[7])

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
  nfolds_u = 2,
  nfolds_v = NULL,
  thresh = 1e-6,
  maxit = 500,
  tuning_iter = 1,
  parallel = FALSE,
  weights = NULL,
  smoothness_type = "Second_order",
  sparse_tuning_type = "hard",
  tuning_order = "Smoothness",
  cv.pick = "min",
  sparse_tuning_u = 0L,
  smooth_tuning_u = 0,
  sparse_tuning_v = as.list(tuning$HPCA$gamma_v),
  smooth_tuning_v = as.list(tuning$HPCA$alpha_v)
)

# The sign of a principal component is not identified, and the sign returned
# can depend on the linear-algebra library. Orient each HPCA component so that
# its scores sum to a positive value (the orientation shown in the paper).
for (pc in 1:2) {
  if (sum(fit_hpca$PCScores[, pc]) < 0) {
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
    thresh = 1e-5,
    maxit = 300,
    tuning_iter = 1,
    parallel = FALSE,
    weights = NULL,
    smoothness_type = "Second_order",
    sparse_tuning_type = "hard",
    tuning_order = "Smoothness",
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
    direction_or_block = c("u (days)", functional_names, "Daily covariates"),
    alpha_smoothness = c(tp$alpha_u, tp$alpha_v),
    gamma_sparsity = c(tp$gamma_u, tp$gamma_v),
    n_entries = c(n_days, sapply(blocks_original, ncol))
  )
})
write_csv(tuning_table, file.path(out_dir, "appliances_tuning_parameters.csv"))

regular_loadings <- imap_dfr(fits, function(fit, method) {
  map_dfr(1:2, function(pc) {
    tibble(
      method = method,
      PC = paste0("PC", pc),
      feature = regular_feature_names,
      loading = as.numeric(fit$PCFunctions[[pc]][[7]])
    )
  })
}) |>
  mutate(is_zero = loading == 0)
write_csv(regular_loadings, file.path(out_dir, "appliances_nonfunctional_loadings.csv"))

zeroed <- regular_loadings |>
  filter(method == "ReHPCA", is_zero) |>
  group_by(PC) |>
  summarize(features = paste(feature, collapse = ", "), .groups = "drop")
cat("Non-functional covariates set to zero by ReHPCA:\n")
cat(sprintf("  %s: %s\n", zeroed$PC, zeroed$features), sep = "")

## 9. Plot data -----------------------------------------------------------------

method_cols <- c("HPCA" = "#D95F02", "ReHPCA" = "#0072B2")

score_long <- imap_dfr(fits, function(fit, method) {
  tibble(
    date = day_complete$day,
    Method = method,
    PC1 = as.numeric(fit$PCScores[, 1]),
    PC2 = as.numeric(fit$PCScores[, 2])
  )
}) |>
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "Score")

functional_long <- imap_dfr(fits, function(fit, method) {
  map_dfr(1:2, function(pc) {
    map_dfr(1:6, function(j) {
      tibble(
        Method = method,
        PC = paste0("PC", pc),
        Functional_Block = functional_names[j],
        Hour = hour_grid,
        Loading = as.numeric(fit$PCFunctions[[pc]][[j]])
      )
    })
  })
}) |>
  mutate(Functional_Block = factor(Functional_Block, levels = functional_names))

# Display adjustment carried over from the original analysis script: in the
# PC2 indoor-humidity panel the HPCA curve is shifted vertically by the mean
# difference (ReHPCA - HPCA) so that the two curve shapes can be compared.
# All other curves are plotted as estimated.
hin_pc2_shift <- functional_long |>
  filter(PC == "PC2", Functional_Block == "Indoor humidity") |>
  select(Hour, Method, Loading) |>
  pivot_wider(names_from = Method, values_from = Loading) |>
  summarize(shift = mean(ReHPCA - HPCA)) |>
  pull(shift)

functional_long <- functional_long |>
  mutate(Loading_plot = if_else(
    PC == "PC2" & Functional_Block == "Indoor humidity" & Method == "HPCA",
    Loading + hin_pc2_shift,
    Loading
  ))

# Readable covariate labels, grouped (weather, month, day type) and ordered
# within group by the largest absolute loading over both methods and PCs.
regular_plot_df <- regular_loadings |>
  mutate(
    month_id = suppressWarnings(as.integer(str_extract(feature, "\\d+"))),
    label = case_when(
      feature == "Pressure_day" ~ "Daily pressure",
      feature == "Windspeed_day" ~ "Daily wind speed",
      feature == "Visibility_day" ~ "Daily visibility",
      feature == "Tdewpoint_day" ~ "Daily dew point",
      feature == "daytype_weekday" ~ "Weekday",
      feature == "daytype_weekend" ~ "Weekend",
      str_detect(feature, "^month_") ~ paste0("Month: ", month.abb[month_id])
    ),
    group = case_when(
      str_detect(feature, "_day$") ~ 1L,
      str_detect(feature, "^month_") ~ 2L,
      TRUE ~ 3L
    )
  )

label_order <- regular_plot_df |>
  group_by(feature, label, group, month_id) |>
  summarize(max_abs = max(abs(loading)), .groups = "drop") |>
  arrange(group, desc(max_abs), month_id, feature) |>
  pull(label)

regular_plot_df <- regular_plot_df |>
  mutate(label = factor(label, levels = rev(label_order)))

## 10. Figure S1: u scores, functional loadings, non-functional loadings ------

compact_theme <- theme_bw(base_size = 9.5) +
  theme(
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey92", color = "grey45", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", linewidth = 0.3),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.5, color = "black")
  )

u_bound <- ceiling(max(abs(score_long$Score)) / 10) * 10

p_u_time <- ggplot(score_long, aes(date, Score, color = Method, group = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.22, color = "grey45") +
  geom_line(linewidth = 0.65, alpha = 0.9) +
  facet_wrap(~ PC, nrow = 1, scales = "fixed") +
  scale_color_manual(values = method_cols, guide = "none") +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-u_bound, u_bound),
                     breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "u direction: PC scores over ordered days", x = NULL, y = "u score") +
  compact_theme +
  theme(
    plot.title = element_text(face = "bold", size = 21, margin = margin(b = 4)),
    strip.text.x = element_text(face = "bold", size = 18),
    axis.title.y.left = element_text(face = "bold", size = 18, margin = margin(r = 5)),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.y.left = element_text(margin = margin(r = 1)),
    axis.ticks.length.y = grid::unit(1.2, "pt"),
    panel.spacing.x = unit(0.05, "cm"),
    plot.margin = margin(2, 2, 1, 0)
  )

p_v_function <- ggplot(functional_long,
                       aes(Hour, Loading_plot, color = Method, group = Method)) +
  geom_hline(yintercept = 0, linewidth = 0.22, color = "grey40") +
  geom_line(linewidth = 0.70, alpha = 0.92) +
  facet_grid(
    Functional_Block ~ PC,
    scales = "free_y",
    switch = "y",
    labeller = labeller(Functional_Block = c(
      "Appliances"       = "Appliance",
      "Lights"           = "Lights",
      "Indoor temp."     = "Indoor\ntemp.",
      "Indoor humidity"  = "Indoor\nhumidity",
      "Outdoor temp."    = "Outdoor\ntemp.",
      "Outdoor humidity" = "Outdoor\nhumidity"
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
    plot.title = element_text(face = "bold", size = 21, hjust = 0.5, margin = margin(b = 5)),
    plot.title.position = "panel",
    strip.text.x = element_text(face = "bold", size = 18, margin = margin(3, 3, 3, 3)),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 17,
                                     margin = margin(4, 5, 4, 5)),
    strip.placement = "outside",
    strip.switch.pad.grid = unit(3, "pt"),
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.04, "cm"),
    axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = 4)),
    axis.title.y.left = element_text(angle = 90, face = "bold", size = 18,
                                     margin = margin(r = 12)),
    axis.text.x = element_text(size = 13, color = "black", margin = margin(t = 3)),
    axis.text.y = element_text(size = 13, color = "black"),
    plot.margin = margin(2, 4, 0, 4)
  )

p_v_regular <- ggplot(regular_plot_df, aes(loading, label, fill = method)) +
  geom_vline(xintercept = 0, linewidth = 0.22, color = "grey40") +
  geom_col(position = position_dodge(width = 0.68), width = 0.58, alpha = 0.92) +
  facet_grid(. ~ PC, scales = "free_x") +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 2),
                     labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0.15, 0.15)),
                     guide = guide_axis(check.overlap = TRUE)) +
  labs(title = "v direction: regular covariate loadings", x = "PC loading", y = NULL) +
  compact_theme +
  theme(
    plot.title = element_text(face = "bold", size = 21, hjust = 0.5, margin = margin(b = 5)),
    plot.title.position = "panel",
    strip.text.x = element_text(face = "bold", size = 18, margin = margin(3, 3, 3, 3)),
    panel.spacing.x = unit(0.35, "cm"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.20),
    axis.title.x = element_text(face = "bold", size = 18, margin = margin(t = 8)),
    axis.text.x = element_text(size = 12, color = "black", margin = margin(t = 4)),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.ticks.length.x = unit(1.5, "pt"),
    plot.margin = margin(2, 4, 0, 3)
  )

# Method legend drawn as a small stand-alone panel (upper left).
p_legend <- ggplot(tibble(Method = names(method_cols), y = c(0.65, 0.35))) +
  geom_segment(aes(x = 0.04, xend = 0.34, y = y, yend = y, color = Method),
               linewidth = 1.7) +
  geom_text(aes(x = 0.42, y = y, label = Method),
            hjust = 0, size = 6.0, fontface = "bold") +
  scale_color_manual(values = method_cols, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

fig_s1 <- ((p_legend | p_u_time) + plot_layout(widths = c(0.17, 1))) /
  ((wrap_elements(full = p_v_function) | wrap_elements(full = p_v_regular)) +
     plot_layout(widths = c(1.65, 1.00))) +
  plot_layout(heights = c(0.55, 1.80))

ggsave(file.path(out_dir, "appliances_compare_u_and_v.pdf"),
       fig_s1, width = 19, height = 13, units = "in")

## 11. Figure S2: mean curves plus/minus the ReHPCA functional PCs -------------

# Raw-scale mean curve +/- sd(u_k) * v_kb * sd(block b), for k = 1, 2.
# Energy curves (Appliances, Lights) are truncated at zero.
mean_pm_df <- map_dfr(1:2, function(k) {
  u_sd <- sd(fit_rehpca$PCScores[, k])
  map_dfr(1:6, function(j) {
    mu  <- colMeans(X_raw[[j]])
    dev <- u_sd * sd(as.vector(X_raw[[j]])) * as.numeric(fit_rehpca$PCFunctions[[k]][[j]])
    tibble(
      Variable = functional_names[j],
      PC = k,
      Hour = hour_grid,
      Mean = mu,
      `Mean + PC` = mu + dev,
      `Mean - PC` = mu - dev
    )
  })
}) |>
  pivot_longer(c(Mean, `Mean + PC`, `Mean - PC`),
               names_to = "Curve", values_to = "Value") |>
  mutate(
    Curve = factor(Curve, levels = c("Mean - PC", "Mean", "Mean + PC")),
    Value = if_else(Variable %in% c("Appliances", "Lights"), pmax(Value, 0), Value)
  )

# Row labels: share of the total sum of squares of the six scaled functional
# blocks captured by the rank-one term u_k v_k' of each ReHPCA component.
total_ss <- sum(sapply(blocks_original[1:6], function(M) sum(M^2)))
pc_pct <- sapply(1:2, function(k) {
  u <- fit_rehpca$PCScores[, k]
  100 * sum(sapply(1:6, function(j) {
    sum(outer(u, as.numeric(fit_rehpca$PCFunctions[[k]][[j]]))^2)
  })) / total_ss
})
pc_labels <- sprintf("PC %d (%.2f%%)", 1:2, pc_pct)

curve_cols <- c("Mean - PC" = "#D62728", "Mean" = "black", "Mean + PC" = "#009E73")
curve_linewidths <- c("Mean - PC" = 0.80, "Mean" = 1.05, "Mean + PC" = 0.80)

# One panel per (PC, variable): variable names on the top row only, PC
# labels on the first column only, x tick marks on the bottom row only.
make_panel <- function(k, j) {
  ggplot(filter(mean_pm_df, PC == k, Variable == functional_names[j]),
         aes(Hour, Value, color = Curve, linewidth = Curve, group = Curve)) +
    geom_line(alpha = 0.95) +
    scale_color_manual(values = curve_cols) +
    scale_linewidth_manual(values = curve_linewidths, guide = "none") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 23)) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.8))) +
    labs(
      title = if (k == 1) functional_names[j] else NULL,
      x = NULL,
      y = if (j == 1) pc_labels[k] else NULL,
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 19, hjust = 0.5, margin = margin(b = 8)),
      axis.title.y = element_text(face = "bold", size = 18, margin = margin(r = 8)),
      axis.text = element_text(size = 15, color = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 21),
      legend.key.width = unit(2, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    if (k == 1) {
      # top row: no tick marks; labels placed as in the published figure
      theme(axis.ticks.x = element_blank(),
            axis.ticks.length.x = unit(0, "pt"),
            axis.text.x = element_text(size = 15, color = "black", margin = margin(t = 1.7, b = -1.7)))
    } else theme()
}

fig_s2 <- wrap_plots(
  lapply(seq_len(12), function(i) make_panel(k = (i - 1) %/% 6 + 1, j = (i - 1) %% 6 + 1)),
  ncol = 6,
  guides = "collect"
) +
  plot_annotation(
    title = "Mean curves plus and minus ReHPCA functional PCs",
    theme = theme(plot.title = element_text(face = "bold", size = 27, hjust = 0.5))
  ) &
  theme(legend.position = "bottom", legend.justification = "center")

ggsave(file.path(out_dir, "appliances_ReHPCA_mean_plus_minus_PC1_PC2.pdf"),
       fig_s2, width = 23, height = 11, units = "in")

## 12. Figure S3: u scores coloured by month, sized by daily energy ------------

month_cols <- c(
  "Jan" = "#E76F51", "Feb" = "#66A61E", "Mar" = "#00A6B4", "Apr" = "#B77CFF",
  "May" = "#E69F00", "Jun" = "#56B4E9", "Jul" = "#009E73", "Aug" = "#CC79A7",
  "Sep" = "#999999", "Oct" = "#D55E00", "Nov" = "#0072B2", "Dec" = "#F0E442"
)

daily_energy_df <- appliance_complete |>
  group_by(date = day) |>
  summarize(
    Daily_energy = sum(Appliances),
    Month = month(first(day), label = TRUE, abbr = TRUE),
    .groups = "drop"
  )

u_score_df <- score_long |>
  left_join(daily_energy_df, by = "date") |>
  mutate(Month = factor(Month, levels = unique(Month)))

fig_s3 <- ggplot(u_score_df, aes(date, Score)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_line(aes(group = interaction(Method, PC)),
            color = "grey45", linewidth = 0.35, alpha = 0.65) +
  geom_point(aes(color = Month, size = Daily_energy), alpha = 0.85) +
  facet_grid(PC ~ Method, scales = "free_y") +
  scale_color_manual(values = month_cols, drop = FALSE,
                     guide = guide_legend(order = 1)) +
  scale_size_continuous(range = c(1.0, 3.5),
                        breaks = pretty(u_score_df$Daily_energy, n = 4),
                        name = "Daily energy", guide = guide_legend(order = 2)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m",
               expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Ordered-day scores in the u direction",
       x = "Date", y = "u score", color = "Month") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 30, hjust = 0.5, margin = margin(b = 15)),
    strip.text = element_text(face = "bold", size = 26),
    axis.title.x = element_text(face = "bold", size = 26, margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", size = 26, margin = margin(r = 12)),
    axis.text.x = element_text(size = 21, color = "black"),
    axis.text.y = element_text(size = 21, color = "black"),
    legend.title = element_text(face = "bold", size = 25),
    legend.text = element_text(size = 22),
    legend.key.height = unit(0.9, "cm"),
    legend.key.width = unit(1.0, "cm"),
    legend.position = "right",
    panel.spacing = unit(0.25, "cm"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(out_dir, "appliances_u_scores_HPCA_vs_ReHPCA.pdf"),
       fig_s3, width = 19, height = 12, units = "in")

## 13. Session information ----------------------------------------------------

cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "s\n")
writeLines(capture.output(sessionInfo()), file.path(out_dir, "appliances_sessionInfo.txt"))
