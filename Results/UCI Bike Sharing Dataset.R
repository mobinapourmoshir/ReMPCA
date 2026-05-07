###############################################################################
# Real multivariate two-way hybrid data example using ReMPCA
#
# Dataset:
#   UCI Bike Sharing Dataset
#
# Hybrid structure:
#   Rows:
#       Days, ordered in time.
#
#   Functional variables:
#       1. 24-hour bike count curve
#       2. 24-hour temperature curve
#       3. 24-hour humidity curve
#       4. 24-hour wind-speed curve
#
#   Non-functional variables:
#       Daily season, month, weekday, holiday, working day, weather situation,
#       and daily aggregated weather variables.
#
# Two-way structure:
#   u direction:
#       ordered days
#
#   v direction:
#       functional 24-hour grids + regular daily covariates
###############################################################################

###############################################################################
# 0. Install and load packages
###############################################################################

packages_needed <- c(
  "remotes",
  "tidyverse",
  "lubridate",
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

install_if_missing(packages_needed)

if (!requireNamespace("ReMPCA", quietly = TRUE)) {
  remotes::install_github("mobinapourmoshir/ReMPCA", force = TRUE)
}

library(ReMPCA)
library(tidyverse)
library(lubridate)
library(gt)
library(patchwork)

###############################################################################
# 1. Create output folders
###############################################################################

dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

###############################################################################
# 2. Download the UCI Bike Sharing Dataset
###############################################################################

zip_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00275/Bike-Sharing-Dataset.zip"
zip_file <- "data/Bike-Sharing-Dataset.zip"

if (!file.exists(zip_file)) {
  download.file(zip_url, destfile = zip_file, mode = "wb")
}

unzip(zip_file, exdir = "data/Bike-Sharing-Dataset")

hour_file <- "data/Bike-Sharing-Dataset/hour.csv"
day_file  <- "data/Bike-Sharing-Dataset/day.csv"

hour_dat <- read.csv(hour_file)
day_dat  <- read.csv(day_file)

hour_dat$dteday <- as.Date(hour_dat$dteday)
day_dat$dteday  <- as.Date(day_dat$dteday)

cat("Hourly data dimensions:\n")
print(dim(hour_dat))

cat("Daily data dimensions:\n")
print(dim(day_dat))

###############################################################################
# 3. Build daily functional curves from hourly data
###############################################################################

# We only keep complete days with all 24 hours available.
complete_days <- hour_dat |>
  group_by(dteday) |>
  summarize(n_hours = n_distinct(hr), .groups = "drop") |>
  filter(n_hours == 24) |>
  pull(dteday)

hour_complete <- hour_dat |>
  filter(dteday %in% complete_days)

day_complete <- day_dat |>
  filter(dteday %in% complete_days) |>
  arrange(dteday)

# Helper function: convert a long hourly variable into a day x 24 matrix.
make_daily_curve_matrix <- function(data, variable_name) {
  data |>
    select(dteday, hr, value = all_of(variable_name)) |>
    arrange(dteday, hr) |>
    mutate(hr = paste0("h", sprintf("%02d", hr))) |>
    pivot_wider(
      names_from = hr,
      values_from = value
    ) |>
    arrange(dteday) |>
    select(-dteday) |>
    as.matrix()
}

X_count_raw <- make_daily_curve_matrix(hour_complete, "cnt")
X_temp_raw  <- make_daily_curve_matrix(hour_complete, "temp")
X_hum_raw   <- make_daily_curve_matrix(hour_complete, "hum")
X_wind_raw  <- make_daily_curve_matrix(hour_complete, "windspeed")

# Check dimensions
cat("\nFunctional block dimensions:\n")
print(dim(X_count_raw))
print(dim(X_temp_raw))
print(dim(X_hum_raw))
print(dim(X_wind_raw))

###############################################################################
# 4. Build non-functional daily covariate matrix
###############################################################################

# Use daily variables as regular/non-functional covariates.
# Categorical variables are one-hot encoded.
regular_df <- day_complete |>
  transmute(
    date = dteday,
    season = factor(season),
    year = factor(yr),
    month = factor(mnth),
    holiday = factor(holiday),
    weekday = factor(weekday),
    workingday = factor(workingday),
    weathersit = factor(weathersit),
    temp_day = temp,
    atemp_day = atemp,
    hum_day = hum,
    windspeed_day = windspeed
  )

X_regular_raw <- model.matrix(
  ~ season + year + month + holiday + weekday + workingday + weathersit +
    temp_day + atemp_day + hum_day + windspeed_day - 1,
  data = regular_df
)

regular_feature_names <- colnames(X_regular_raw)

cat("\nRegular block dimensions:\n")
print(dim(X_regular_raw))

###############################################################################
# 5. Scale blocks
###############################################################################

# For hybrid PCA, scaling is important so that bike count does not dominate
# weather variables just because it has a larger numerical range.

scale_matrix_global <- function(M) {
  M <- as.matrix(M)
  (M - mean(M, na.rm = TRUE)) / sd(as.vector(M), na.rm = TRUE)
}

scale_matrix_columns <- function(M) {
  M <- as.matrix(M)
  M_scaled <- scale(M, center = TRUE, scale = TRUE)
  M_scaled[is.na(M_scaled)] <- 0
  as.matrix(M_scaled)
}

X_count <- scale_matrix_global(X_count_raw)
X_temp  <- scale_matrix_global(X_temp_raw)
X_hum   <- scale_matrix_global(X_hum_raw)
X_wind  <- scale_matrix_global(X_wind_raw)

X_regular <- scale_matrix_columns(X_regular_raw)

n_days <- nrow(X_count)
hour_grid <- seq(0, 23, length.out = 24)
day_grid <- seq(0, 1, length.out = n_days)

###############################################################################
# 6. Create ReMPCA objects
###############################################################################

# Tuning grids.
# These are moderate grids, useful for a real-data example.
# You can make them denser later for the final paper.
smooth_grid_v <- 2^seq(-12, 4, length.out = 10)
smooth_grid_u <- 2^seq(-12, 4, length.out = 10)

sparse_grid_u <- unique(round(seq(0, floor(0.50 * n_days), length.out = 10)))

sparse_grid_24 <- unique(round(seq(0, 23, length.out = 10)))
sparse_grid_regular <- unique(round(seq(0, ncol(X_regular) - 1, length.out = 12)))

fd_count <- fdClass(
  data = X_count,
  argval = hour_grid,
  Smoothing_parameter = smooth_grid_v,
  Sparsity_parameter = sparse_grid_24
)

fd_temp <- fdClass(
  data = X_temp,
  argval = hour_grid,
  Smoothing_parameter = smooth_grid_v,
  Sparsity_parameter = sparse_grid_24
)

fd_hum <- fdClass(
  data = X_hum,
  argval = hour_grid,
  Smoothing_parameter = smooth_grid_v,
  Sparsity_parameter = sparse_grid_24
)

fd_wind <- fdClass(
  data = X_wind,
  argval = hour_grid,
  Smoothing_parameter = smooth_grid_v,
  Sparsity_parameter = sparse_grid_24
)

rd_daily <- rdClass(
  data = X_regular,
  Sparsity_parameter = sparse_grid_regular
)

bike_hd <- hdClass(
  hdlist = list(
    Count_curve = fd_count,
    Temperature_curve = fd_temp,
    Humidity_curve = fd_hum,
    Wind_curve = fd_wind,
    Daily_covariates = rd_daily
  ),
  argval = day_grid,
  Smoothing_parameter = smooth_grid_u,
  Sparsity_parameter = sparse_grid_u
)

###############################################################################
# 7. Fit two-way multivariate hybrid ReMPCA
###############################################################################

set.seed(2026)

bike_fit <- ReMPCA(
  hd = bike_hd,
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
  sparse_tuning_v = NULL,
  smooth_tuning_u = NULL,
  smooth_tuning_v = NULL
)

print(bike_fit)

###############################################################################
# 8. Extract results
###############################################################################

scores <- as.data.frame(bike_fit$PCScores)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))

score_df <- bind_cols(
  tibble(
    date = day_complete$dteday,
    season = factor(day_complete$season),
    month = factor(month(day_complete$dteday)),
    weekday = factor(wday(day_complete$dteday, label = TRUE)),
    workingday = factor(day_complete$workingday),
    holiday = factor(day_complete$holiday)
  ),
  scores
)

# Extract PC functions/loadings.
# bike_fit$PCFunctions[[pc]][[block]]
get_pc_block <- function(pc, block) {
  as.numeric(bike_fit$PCFunctions[[pc]][[block]])
}

pc1_count <- get_pc_block(1, 1)
pc1_temp  <- get_pc_block(1, 2)
pc1_hum   <- get_pc_block(1, 3)
pc1_wind  <- get_pc_block(1, 4)
pc1_reg   <- get_pc_block(1, 5)

pc2_count <- get_pc_block(2, 1)
pc2_temp  <- get_pc_block(2, 2)
pc2_hum   <- get_pc_block(2, 3)
pc2_wind  <- get_pc_block(2, 4)
pc2_reg   <- get_pc_block(2, 5)

###############################################################################
# 9. Tables
###############################################################################

safe_collapse <- function(x, digits = 4) {
  if (is.null(x)) return(NA_character_)
  x <- unlist(x)
  if (length(x) == 0) return(NA_character_)
  paste(signif(as.numeric(x), digits), collapse = ", ")
}

block_summary <- tibble(
  Block = c(
    "Hourly count curve",
    "Hourly temperature curve",
    "Hourly humidity curve",
    "Hourly wind-speed curve",
    "Daily regular covariates"
  ),
  Type = c("Functional", "Functional", "Functional", "Functional", "Regular"),
  Rows = n_days,
  Columns = c(
    ncol(X_count),
    ncol(X_temp),
    ncol(X_hum),
    ncol(X_wind),
    ncol(X_regular)
  ),
  Row_Domain = "Ordered days",
  Column_Domain = c(
    "Hour of day",
    "Hour of day",
    "Hour of day",
    "Hour of day",
    "Daily scalar/categorical variables"
  )
)

tuning_summary <- tibble(
  Quantity = c(
    "OptimalAlphaU",
    "OptimalGammaU",
    "OptimalAlphaV",
    "OptimalGammaV",
    "VarianceExplained"
  ),
  Value = c(
    safe_collapse(bike_fit$OptimalAlphaU),
    safe_collapse(bike_fit$OptimalGammaU),
    safe_collapse(bike_fit$OptimalAlphaV),
    safe_collapse(bike_fit$OptimalGammaV),
    safe_collapse(bike_fit$VarianceExplained)
  )
)

top_regular_pc1 <- tibble(
  Feature = regular_feature_names,
  Loading = pc1_reg,
  AbsLoading = abs(pc1_reg)
) |>
  arrange(desc(AbsLoading)) |>
  slice_head(n = 15)

top_regular_pc2 <- tibble(
  Feature = regular_feature_names,
  Loading = pc2_reg,
  AbsLoading = abs(pc2_reg)
) |>
  arrange(desc(AbsLoading)) |>
  slice_head(n = 15)

write.csv(block_summary, "tables/bike_block_summary.csv", row.names = FALSE)
write.csv(tuning_summary, "tables/bike_tuning_summary.csv", row.names = FALSE)
write.csv(top_regular_pc1, "tables/bike_top_regular_loadings_PC1.csv", row.names = FALSE)
write.csv(top_regular_pc2, "tables/bike_top_regular_loadings_PC2.csv", row.names = FALSE)
write.csv(score_df, "tables/bike_scores_by_day.csv", row.names = FALSE)

block_summary |>
  gt() |>
  tab_header(
    title = "Hybrid Data Structure",
    subtitle = "UCI Bike Sharing Dataset transformed into multivariate two-way hybrid data"
  ) |>
  gtsave("tables/bike_block_summary.html")

tuning_summary |>
  gt() |>
  tab_header(
    title = "Selected ReMPCA Tuning Parameters"
  ) |>
  gtsave("tables/bike_tuning_summary.html")

top_regular_pc1 |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  gt() |>
  tab_header(
    title = "Top Daily Covariates for PC1"
  ) |>
  gtsave("tables/bike_top_regular_loadings_PC1.html")

top_regular_pc2 |>
  mutate(across(where(is.numeric), ~ round(.x, 4))) |>
  gt() |>
  tab_header(
    title = "Top Daily Covariates for PC2"
  ) |>
  gtsave("tables/bike_top_regular_loadings_PC2.html")

###############################################################################
# 10. Plot PC scores over calendar time
###############################################################################

p_score_time_pc1 <- ggplot(score_df, aes(x = date, y = PC1)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(color = season), size = 1.4, alpha = 0.8) +
  labs(
    title = "PC1 Scores Over Time",
    subtitle = "Row direction is ordered by day",
    x = "Date",
    y = "PC1 score",
    color = "Season"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_score_time_pc2 <- ggplot(score_df, aes(x = date, y = PC2)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(color = season), size = 1.4, alpha = 0.8) +
  labs(
    title = "PC2 Scores Over Time",
    subtitle = "Second hybrid component",
    x = "Date",
    y = "PC2 score",
    color = "Season"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  "figures/bike_PC1_scores_over_time.png",
  p_score_time_pc1,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/bike_PC2_scores_over_time.png",
  p_score_time_pc2,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  "figures/bike_PC_scores_over_time.pdf",
  p_score_time_pc1 / p_score_time_pc2,
  width = 10,
  height = 8
)

###############################################################################
# 11. Plot score distributions by season and working day
###############################################################################

p_pc1_season <- ggplot(score_df, aes(x = season, y = PC1, fill = season)) +
  geom_boxplot(alpha = 0.75) +
  labs(
    title = "PC1 Score Distribution by Season",
    x = "Season",
    y = "PC1 score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

p_pc1_workingday <- ggplot(score_df, aes(x = workingday, y = PC1, fill = workingday)) +
  geom_boxplot(alpha = 0.75) +
  labs(
    title = "PC1 Score Distribution by Working Day",
    x = "Working day",
    y = "PC1 score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(
  "figures/bike_PC1_score_boxplots.pdf",
  p_pc1_season + p_pc1_workingday,
  width = 11,
  height = 5
)

ggsave(
  "figures/bike_PC1_score_boxplots.png",
  p_pc1_season + p_pc1_workingday,
  width = 11,
  height = 5,
  dpi = 300
)

###############################################################################
# 12. Plot functional loadings for PC1 and PC2
###############################################################################

loading_df <- bind_rows(
  tibble(
    PC = "PC1",
    Hour = hour_grid,
    Count = pc1_count,
    Temperature = pc1_temp,
    Humidity = pc1_hum,
    Wind = pc1_wind
  ),
  tibble(
    PC = "PC2",
    Hour = hour_grid,
    Count = pc2_count,
    Temperature = pc2_temp,
    Humidity = pc2_hum,
    Wind = pc2_wind
  )
) |>
  pivot_longer(
    cols = c(Count, Temperature, Humidity, Wind),
    names_to = "Functional_Block",
    values_to = "Loading"
  )

p_loading <- ggplot(
  loading_df,
  aes(x = Hour, y = Loading, linetype = PC)
) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Functional_Block, scales = "free_y", ncol = 2) +
  labs(
    title = "Estimated Functional Loadings",
    subtitle = "Hourly loading curves for each functional block",
    x = "Hour of day",
    y = "Loading"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  "figures/bike_functional_loadings_PC1_PC2.pdf",
  p_loading,
  width = 10,
  height = 7
)

ggsave(
  "figures/bike_functional_loadings_PC1_PC2.png",
  p_loading,
  width = 10,
  height = 7,
  dpi = 300
)

###############################################################################
# 13. Plot regular covariate loadings
###############################################################################

regular_loading_df <- bind_rows(
  tibble(
    PC = "PC1",
    Feature = regular_feature_names,
    Loading = pc1_reg
  ),
  tibble(
    PC = "PC2",
    Feature = regular_feature_names,
    Loading = pc2_reg
  )
) |>
  group_by(PC) |>
  mutate(AbsLoading = abs(Loading)) |>
  arrange(PC, desc(AbsLoading)) |>
  slice_head(n = 15) |>
  ungroup() |>
  mutate(Feature = reorder(Feature, AbsLoading))

p_regular <- ggplot(
  regular_loading_df,
  aes(x = Feature, y = Loading, fill = PC)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(
    title = "Top Non-Functional Covariate Loadings",
    x = "Daily covariate",
    y = "Loading"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  "figures/bike_top_regular_covariate_loadings.pdf",
  p_regular,
  width = 10,
  height = 7
)

ggsave(
  "figures/bike_top_regular_covariate_loadings.png",
  p_regular,
  width = 10,
  height = 7,
  dpi = 300
)

###############################################################################
# 14. Plot original functional heatmaps
###############################################################################

matrix_to_heatmap_df <- function(M, block_name) {
  as.data.frame(as.table(M)) |>
    as_tibble() |>
    transmute(
      Day_Index = as.integer(Var1),
      Hour = as.integer(Var2) - 1,
      Value = as.numeric(Freq),
      Block = block_name
    )
}

heat_df <- bind_rows(
  matrix_to_heatmap_df(X_count, "Count"),
  matrix_to_heatmap_df(X_temp, "Temperature"),
  matrix_to_heatmap_df(X_hum, "Humidity"),
  matrix_to_heatmap_df(X_wind, "Wind speed")
)

p_heat <- ggplot(heat_df, aes(x = Hour, y = Day_Index, fill = Value)) +
  geom_raster() +
  scale_y_reverse() +
  facet_wrap(~ Block, scales = "free", ncol = 2) +
  labs(
    title = "Functional Blocks as Day-by-Hour Heatmaps",
    subtitle = "Scaled data used in the hybrid PCA model",
    x = "Hour of day",
    y = "Day index",
    fill = "Scaled value"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "figures/bike_functional_block_heatmaps.pdf",
  p_heat,
  width = 11,
  height = 8
)

ggsave(
  "figures/bike_functional_block_heatmaps.png",
  p_heat,
  width = 11,
  height = 8,
  dpi = 300
)

###############################################################################
# 15. Reconstruct first two PCs and plot reconstruction heatmaps
###############################################################################

reconstruct_block <- function(block_index, npc = 2) {
  reconstructed <- matrix(
    0,
    nrow = n_days,
    ncol = length(bike_fit$PCFunctions[[1]][[block_index]])
  )

  for (pc in seq_len(npc)) {
    u_pc <- as.numeric(bike_fit$PCScores[, pc])
    v_pc <- as.numeric(bike_fit$PCFunctions[[pc]][[block_index]])
    reconstructed <- reconstructed + outer(u_pc, v_pc)
  }

  reconstructed
}

X_count_recon <- reconstruct_block(1, npc = 2)
X_temp_recon  <- reconstruct_block(2, npc = 2)
X_hum_recon   <- reconstruct_block(3, npc = 2)
X_wind_recon  <- reconstruct_block(4, npc = 2)

recon_heat_df <- bind_rows(
  matrix_to_heatmap_df(X_count, "Count") |> mutate(Version = "Observed"),
  matrix_to_heatmap_df(X_count_recon, "Count") |> mutate(Version = "Reconstructed"),
  matrix_to_heatmap_df(X_temp, "Temperature") |> mutate(Version = "Observed"),
  matrix_to_heatmap_df(X_temp_recon, "Temperature") |> mutate(Version = "Reconstructed"),
  matrix_to_heatmap_df(X_hum, "Humidity") |> mutate(Version = "Observed"),
  matrix_to_heatmap_df(X_hum_recon, "Humidity") |> mutate(Version = "Reconstructed"),
  matrix_to_heatmap_df(X_wind, "Wind speed") |> mutate(Version = "Observed"),
  matrix_to_heatmap_df(X_wind_recon, "Wind speed") |> mutate(Version = "Reconstructed")
)

p_recon_heat <- ggplot(recon_heat_df, aes(x = Hour, y = Day_Index, fill = Value)) +
  geom_raster() +
  scale_y_reverse() +
  facet_grid(Block ~ Version, scales = "free") +
  labs(
    title = "Observed and Reconstructed Functional Blocks",
    subtitle = "Reconstruction using the first two hybrid principal components",
    x = "Hour of day",
    y = "Day index",
    fill = "Scaled value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "figures/bike_observed_vs_reconstructed_heatmaps.pdf",
  p_recon_heat,
  width = 12,
  height = 10
)

ggsave(
  "figures/bike_observed_vs_reconstructed_heatmaps.png",
  p_recon_heat,
  width = 12,
  height = 10,
  dpi = 300
)

###############################################################################
# 16. Reconstruction error table
###############################################################################

rmse <- function(A, B) {
  sqrt(mean((as.numeric(A) - as.numeric(B))^2))
}

reconstruction_table <- tibble(
  Block = c("Count", "Temperature", "Humidity", "Wind speed"),
  RMSE_Observed_vs_Reconstructed = c(
    rmse(X_count, X_count_recon),
    rmse(X_temp, X_temp_recon),
    rmse(X_hum, X_hum_recon),
    rmse(X_wind, X_wind_recon)
  )
)

write.csv(
  reconstruction_table,
  "tables/bike_reconstruction_error_table.csv",
  row.names = FALSE
)

reconstruction_table |>
  mutate(RMSE_Observed_vs_Reconstructed = round(RMSE_Observed_vs_Reconstructed, 4)) |>
  gt() |>
  tab_header(
    title = "Functional Block Reconstruction Error",
    subtitle = "Observed versus two-component ReMPCA reconstruction"
  ) |>
  gtsave("tables/bike_reconstruction_error_table.html")

###############################################################################
# 17. Built-in ReMPCA plots
###############################################################################

pdf("figures/bike_builtin_ReMPCA_plots.pdf", width = 11, height = 8.5)

try(plot_pc_functions(bike_fit), silent = TRUE)
try(plot_pc_scores(bike_fit), silent = TRUE)
try(plot_cv_u(bike_fit), silent = TRUE)
try(plot_cv_v(bike_fit), silent = TRUE)
try(plot_gcv_u(bike_fit), silent = TRUE)
try(plot_gcv_v(bike_fit), silent = TRUE)

dev.off()

###############################################################################
# 18. Final output
###############################################################################

cat("\nDone.\n")
cat("Figures saved in: figures/\n")
cat("Tables saved in: tables/\n")
cat("\nMain fitted object: bike_fit\n")
cat("Main hybrid object: bike_hd\n")
cat("Main score table: score_df\n")
cat("Main block table: block_summary\n")
cat("Main tuning table: tuning_summary\n")
