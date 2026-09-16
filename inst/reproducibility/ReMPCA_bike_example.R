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

###############################################################################
# 0. Load packages
###############################################################################

library(ReMPCA)
library(tidyverse)
library(lubridate)
library(patchwork)

###############################################################################
# 1. Load package data
###############################################################################

load(file.path("data", "bike_day.rda"))
load(file.path("data", "bike_hour.rda"))

day_dat  <- bike_day
hour_dat <- bike_hour

day_dat$dteday  <- as.Date(day_dat$dteday)
hour_dat$dteday <- as.Date(hour_dat$dteday)

cat("Hourly data dimensions:\n")
print(dim(hour_dat))

cat("Daily data dimensions:\n")
print(dim(day_dat))

###############################################################################
# 2. Build daily functional curves from hourly data
###############################################################################

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

# Count
X_count_raw <- hour_complete |>
  select(dteday, hr, value = cnt) |>
  arrange(dteday, hr) |>
  mutate(hr = paste0("h", sprintf("%02d", hr))) |>
  pivot_wider(names_from = hr, values_from = value) |>
  arrange(dteday) |>
  select(-dteday) |>
  as.matrix()

# Temperature
X_temp_raw <- hour_complete |>
  select(dteday, hr, value = temp) |>
  arrange(dteday, hr) |>
  mutate(hr = paste0("h", sprintf("%02d", hr))) |>
  pivot_wider(names_from = hr, values_from = value) |>
  arrange(dteday) |>
  select(-dteday) |>
  as.matrix()

# Humidity
X_hum_raw <- hour_complete |>
  select(dteday, hr, value = hum) |>
  arrange(dteday, hr) |>
  mutate(hr = paste0("h", sprintf("%02d", hr))) |>
  pivot_wider(names_from = hr, values_from = value) |>
  arrange(dteday) |>
  select(-dteday) |>
  as.matrix()

# Wind speed
X_wind_raw <- hour_complete |>
  select(dteday, hr, value = windspeed) |>
  arrange(dteday, hr) |>
  mutate(hr = paste0("h", sprintf("%02d", hr))) |>
  pivot_wider(names_from = hr, values_from = value) |>
  arrange(dteday) |>
  select(-dteday) |>
  as.matrix()

cat("\nFunctional block dimensions:\n")
print(dim(X_count_raw))
print(dim(X_temp_raw))
print(dim(X_hum_raw))
print(dim(X_wind_raw))

###############################################################################
# 3. Build non-functional daily covariate matrix
###############################################################################

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
# 4. Scale blocks
###############################################################################

# Global scaling for functional blocks

X_count <- as.matrix(X_count_raw)
count_sd <- sd(as.vector(X_count), na.rm = TRUE)
if (!is.finite(count_sd) || count_sd == 0) count_sd <- 1
X_count <- (X_count - mean(X_count, na.rm = TRUE)) / count_sd

X_temp <- as.matrix(X_temp_raw)
temp_sd <- sd(as.vector(X_temp), na.rm = TRUE)
if (!is.finite(temp_sd) || temp_sd == 0) temp_sd <- 1
X_temp <- (X_temp - mean(X_temp, na.rm = TRUE)) / temp_sd

X_hum <- as.matrix(X_hum_raw)
hum_sd <- sd(as.vector(X_hum), na.rm = TRUE)
if (!is.finite(hum_sd) || hum_sd == 0) hum_sd <- 1
X_hum <- (X_hum - mean(X_hum, na.rm = TRUE)) / hum_sd

X_wind <- as.matrix(X_wind_raw)
wind_sd <- sd(as.vector(X_wind), na.rm = TRUE)
if (!is.finite(wind_sd) || wind_sd == 0) wind_sd <- 1
X_wind <- (X_wind - mean(X_wind, na.rm = TRUE)) / wind_sd

# Columnwise scaling for the non-functional block.
X_regular <- scale(
  as.matrix(X_regular_raw),
  center = TRUE,
  scale = TRUE
)
X_regular[is.na(X_regular)] <- 0
X_regular <- as.matrix(X_regular)

n_days <- nrow(X_count)

hour_grid <- seq(0, 23, length.out = 24)
day_grid  <- seq(0, 1, length.out = n_days)

blocks_original <- list(
  Count_curve = X_count,
  Temperature_curve = X_temp,
  Humidity_curve = X_hum,
  Wind_curve = X_wind,
  Daily_covariates = X_regular
)


###############################################################################
# 5. Fit strict no-penalty baseline with two PCs
###############################################################################

fd_count_base <- fdClass(
  data = blocks_original$Count_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_temp_base <- fdClass(
  data = blocks_original$Temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_hum_base <- fdClass(
  data = blocks_original$Humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_wind_base <- fdClass(
  data = blocks_original$Wind_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

rd_daily_base <- rdClass(
  data = blocks_original$Daily_covariates,
  Sparsity_parameter = 0L
)

bike_hd_base <- hdClass(
  hdlist = list(
    Count_curve = fd_count_base,
    Temperature_curve = fd_temp_base,
    Humidity_curve = fd_hum_base,
    Wind_curve = fd_wind_base,
    Daily_covariates = rd_daily_base
  ),
  argval = day_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

set.seed(2026)

cat("\nFitting strict no-penalty baseline with two PCs...\n")

bike_fit_base <- ReMPCA(
  hd = bike_hd_base,
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

###############################################################################
# 6. Fit PC1 on the original data
###############################################################################

fd_count_pc1 <- fdClass(
  data = blocks_original$Count_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-2,
  Sparsity_parameter = 0L
)

fd_temp_pc1 <- fdClass(
  data = blocks_original$Temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-6,
  Sparsity_parameter = 2L
)

fd_hum_pc1 <- fdClass(
  data = blocks_original$Humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-5,
  Sparsity_parameter = 0L
)

fd_wind_pc1 <- fdClass(
  data = blocks_original$Wind_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-4,
  Sparsity_parameter = 0L
)

rd_daily_pc1 <- rdClass(
  data = blocks_original$Daily_covariates,
  Sparsity_parameter = 8L
)

bike_hd_pc1 <- hdClass(
  hdlist = list(
    Count_curve = fd_count_pc1,
    Temperature_curve = fd_temp_pc1,
    Humidity_curve = fd_hum_pc1,
    Wind_curve = fd_wind_pc1,
    Daily_covariates = rd_daily_pc1
  ),
  argval = day_grid,
  Smoothing_parameter = 2^-6,
  Sparsity_parameter = 1L
)

set.seed(2026)

cat("\nFitting PC1: two-way smooth/sparse ReMPCA...\n")

bike_fit_pc1 <- ReMPCA(
  hd = bike_hd_pc1,
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
  sparse_tuning_u = 1L,
  smooth_tuning_u = 2^-6,
  sparse_tuning_v = list(
    0L,
    2L,
    0L,
    0L,
    8L
  ),
  smooth_tuning_v = list(
    2^-2,
    2^-6,
    2^-5,
    2^-4,
    0
  )
)

###############################################################################
# 7. Deflate the original data using PC1
###############################################################################

u_pc1 <- as.numeric(bike_fit_pc1$PCScores[, 1])

pc1_reconstruction <- list(
  Count_curve = outer(
    u_pc1,
    as.numeric(bike_fit_pc1$PCFunctions[[1]][[1]])
  ),
  Temperature_curve = outer(
    u_pc1,
    as.numeric(bike_fit_pc1$PCFunctions[[1]][[2]])
  ),
  Humidity_curve = outer(
    u_pc1,
    as.numeric(bike_fit_pc1$PCFunctions[[1]][[3]])
  ),
  Wind_curve = outer(
    u_pc1,
    as.numeric(bike_fit_pc1$PCFunctions[[1]][[4]])
  ),
  Daily_covariates = outer(
    u_pc1,
    as.numeric(bike_fit_pc1$PCFunctions[[1]][[5]])
  )
)

blocks_residual_after_pc1 <- list(
  Count_curve =
    blocks_original$Count_curve - pc1_reconstruction$Count_curve,
  Temperature_curve =
    blocks_original$Temperature_curve - pc1_reconstruction$Temperature_curve,
  Humidity_curve =
    blocks_original$Humidity_curve - pc1_reconstruction$Humidity_curve,
  Wind_curve =
    blocks_original$Wind_curve - pc1_reconstruction$Wind_curve,
  Daily_covariates =
    blocks_original$Daily_covariates - pc1_reconstruction$Daily_covariates
)

###############################################################################
# 8. Fit PC2 on the residual data
###############################################################################

fd_count_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Count_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-3,
  Sparsity_parameter = 1L
)

fd_temp_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-4,
  Sparsity_parameter = 0L
)

fd_hum_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-4,
  Sparsity_parameter = 10L
)

fd_wind_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Wind_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-4,
  Sparsity_parameter = 0L
)

rd_daily_pc2 <- rdClass(
  data = blocks_residual_after_pc1$Daily_covariates,
  Sparsity_parameter = 11L
)

bike_hd_pc2 <- hdClass(
  hdlist = list(
    Count_curve = fd_count_pc2,
    Temperature_curve = fd_temp_pc2,
    Humidity_curve = fd_hum_pc2,
    Wind_curve = fd_wind_pc2,
    Daily_covariates = rd_daily_pc2
  ),
  argval = day_grid,
  Smoothing_parameter = 2^-7,
  Sparsity_parameter = 0L
)

set.seed(2026)

cat("\nFitting PC2: residual two-way smooth/sparse ReMPCA...\n")

bike_fit_pc2_residual <- ReMPCA(
  hd = bike_hd_pc2,
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
  sparse_tuning_u = 0L,
  smooth_tuning_u = 2^-7,
  sparse_tuning_v = list(
    1L,
    0L,
    10L,
    0L,
    11L
  ),
  smooth_tuning_v = list(
    2^-3,
    2^-4,
    2^-4,
    2^-4,
    0
  )
)

###############################################################################
# 9. Combine PC1 and PC2 and align signs with the baseline
###############################################################################

bike_fit_pc_specific <- list(
  PCScores = cbind(
    PC1 = as.numeric(bike_fit_pc1$PCScores[, 1]),
    PC2 = as.numeric(bike_fit_pc2_residual$PCScores[, 1])
  ),
  PCFunctions = list(
    bike_fit_pc1$PCFunctions[[1]],
    bike_fit_pc2_residual$PCFunctions[[1]]
  )
)

for (pc in seq_len(2)) {

  cc <- suppressWarnings(
    cor(
      as.numeric(bike_fit_base$PCScores[, pc]),
      as.numeric(bike_fit_pc_specific$PCScores[, pc]),
      use = "complete.obs"
    )
  )

  if (is.finite(cc) && cc < 0) {

    bike_fit_pc_specific$PCScores[, pc] <-
      -bike_fit_pc_specific$PCScores[, pc]

    for (b in seq_along(bike_fit_pc_specific$PCFunctions[[pc]])) {
      bike_fit_pc_specific$PCFunctions[[pc]][[b]] <-
        -bike_fit_pc_specific$PCFunctions[[pc]][[b]]
    }
  }
}

models <- list(
  "No penalty" = bike_fit_base,
  "PC-specific two-way" = bike_fit_pc_specific
)

###############################################################################
# 10. Extract scores and loadings
###############################################################################

score_base <- as.data.frame(
  bike_fit_base$PCScores[, seq_len(2), drop = FALSE]
)
colnames(score_base) <- paste0("PC", seq_len(2))

score_specific <- as.data.frame(
  bike_fit_pc_specific$PCScores[, seq_len(2), drop = FALSE]
)
colnames(score_specific) <- paste0("PC", seq_len(2))

score_df_base <- bind_cols(
  tibble(
    date = day_complete$dteday,
    season = factor(day_complete$season),
    month = factor(month(day_complete$dteday)),
    weekday = factor(wday(day_complete$dteday, label = TRUE)),
    workingday = factor(day_complete$workingday),
    holiday = factor(day_complete$holiday),
    Method = "No penalty"
  ),
  score_base
)

score_df_specific <- bind_cols(
  tibble(
    date = day_complete$dteday,
    season = factor(day_complete$season),
    month = factor(month(day_complete$dteday)),
    weekday = factor(wday(day_complete$dteday, label = TRUE)),
    workingday = factor(day_complete$workingday),
    holiday = factor(day_complete$holiday),
    Method = "PC-specific two-way"
  ),
  score_specific
)

score_df <- bind_rows(
  score_df_base,
  score_df_specific
)

# Functional loading curves: no-penalty baseline.
functional_loading_base_pc1 <- tibble(
  Method = "No penalty",
  PC = "PC1",
  Hour = hour_grid,
  Count = as.numeric(bike_fit_base$PCFunctions[[1]][[1]]),
  Temperature = as.numeric(bike_fit_base$PCFunctions[[1]][[2]]),
  Humidity = as.numeric(bike_fit_base$PCFunctions[[1]][[3]]),
  Wind = as.numeric(bike_fit_base$PCFunctions[[1]][[4]])
)

functional_loading_base_pc2 <- tibble(
  Method = "No penalty",
  PC = "PC2",
  Hour = hour_grid,
  Count = as.numeric(bike_fit_base$PCFunctions[[2]][[1]]),
  Temperature = as.numeric(bike_fit_base$PCFunctions[[2]][[2]]),
  Humidity = as.numeric(bike_fit_base$PCFunctions[[2]][[3]]),
  Wind = as.numeric(bike_fit_base$PCFunctions[[2]][[4]])
)

# Functional loading curves
functional_loading_specific_pc1 <- tibble(
  Method = "PC-specific two-way",
  PC = "PC1",
  Hour = hour_grid,
  Count = as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[1]]),
  Temperature = as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[2]]),
  Humidity = as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[3]]),
  Wind = as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[4]])
)

functional_loading_specific_pc2 <- tibble(
  Method = "PC-specific two-way",
  PC = "PC2",
  Hour = hour_grid,
  Count = as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[1]]),
  Temperature = as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[2]]),
  Humidity = as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[3]]),
  Wind = as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[4]])
)

functional_loading_df <- bind_rows(
  functional_loading_base_pc1,
  functional_loading_base_pc2,
  functional_loading_specific_pc1,
  functional_loading_specific_pc2
) |>
  pivot_longer(
    cols = c(Count, Temperature, Humidity, Wind),
    names_to = "Functional_Block",
    values_to = "Loading"
  )

# Regular covariate loadings.
regular_loading_all <- bind_rows(
  tibble(
    Method = "No penalty",
    PC = "PC1",
    Feature = regular_feature_names,
    Loading = as.numeric(bike_fit_base$PCFunctions[[1]][[5]])
  ),
  tibble(
    Method = "No penalty",
    PC = "PC2",
    Feature = regular_feature_names,
    Loading = as.numeric(bike_fit_base$PCFunctions[[2]][[5]])
  ),
  tibble(
    Method = "PC-specific two-way",
    PC = "PC1",
    Feature = regular_feature_names,
    Loading = as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[5]])
  ),
  tibble(
    Method = "PC-specific two-way",
    PC = "PC2",
    Feature = regular_feature_names,
    Loading = as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[5]])
  )
)

regular_feature_order <- regular_loading_all |>
  filter(Method == "No penalty") |>
  group_by(Feature) |>
  summarize(
    MaxAbsLoading = max(abs(Loading), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(MaxAbsLoading)) |>
  pull(Feature)

regular_loading_df <- regular_loading_all |>
  mutate(
    Feature = factor(
      Feature,
      levels = rev(regular_feature_order)
    )
  )

###############################################################################
# 11. Penalty and sparsity summaries
###############################################################################

penalty_by_pc_block <- tibble(
  Method = "PC-specific two-way",
  PC = rep(c("PC1", "PC2"), each = 6),
  Direction_or_Block = rep(
    c(
      "u direction: ordered days",
      "v block: Count curve",
      "v block: Temperature curve",
      "v block: Humidity curve",
      "v block: Wind curve",
      "v block: Daily regular covariates"
    ),
    times = 2
  ),
  Alpha_smoothness = c(
    2^-6,
    2^-2,
    2^-6,
    2^-5,
    2^-4,
    0,
    2^-7,
    2^-3,
    2^-4,
    2^-4,
    2^-4,
    0
  ),
  Gamma_sparsity = c(
    1L,
    0L,
    2L,
    0L,
    0L,
    8L,
    0L,
    1L,
    0L,
    10L,
    0L,
    11L
  )
)

#print(penalty_by_pc_block)

zero_regular_summary <- regular_loading_all |>
  group_by(Method, PC) |>
  summarize(
    Total_regular_covariates = n(),
    Zero_abs_lt_1e_6 = sum(abs(Loading) < 1e-6),
    Near_zero_abs_lt_1e_3 = sum(abs(Loading) < 1e-3),
    Nonzero_abs_ge_1e_3 = sum(abs(Loading) >= 1e-3),
    .groups = "drop"
  )

#print(zero_regular_summary)

zero_functional_summary <- functional_loading_df |>
  group_by(Method, PC, Functional_Block) |>
  summarize(
    Total_hourly_points = n(),
    Zero_abs_lt_1e_6 = sum(abs(Loading) < 1e-6),
    Near_zero_abs_lt_1e_3 = sum(abs(Loading) < 1e-3),
    Nonzero_abs_ge_1e_3 = sum(abs(Loading) >= 1e-3),
    .groups = "drop"
  )

#print(zero_functional_summary)

# Fit-output summary, without defining a helper function.
base_alpha_u <- if (is.null(bike_fit_base$OptimalAlphaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_base$OptimalAlphaU)), 4), collapse = ", ")
}

pc1_alpha_u <- if (is.null(bike_fit_pc1$OptimalAlphaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc1$OptimalAlphaU)), 4), collapse = ", ")
}

pc2_alpha_u <- if (is.null(bike_fit_pc2_residual$OptimalAlphaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc2_residual$OptimalAlphaU)), 4), collapse = ", ")
}

base_gamma_u <- if (is.null(bike_fit_base$OptimalGammaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_base$OptimalGammaU)), 4), collapse = ", ")
}

pc1_gamma_u <- if (is.null(bike_fit_pc1$OptimalGammaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc1$OptimalGammaU)), 4), collapse = ", ")
}

pc2_gamma_u <- if (is.null(bike_fit_pc2_residual$OptimalGammaU)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc2_residual$OptimalGammaU)), 4), collapse = ", ")
}

base_alpha_v <- if (is.null(bike_fit_base$OptimalAlphaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_base$OptimalAlphaV)), 4), collapse = ", ")
}

pc1_alpha_v <- if (is.null(bike_fit_pc1$OptimalAlphaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc1$OptimalAlphaV)), 4), collapse = ", ")
}

pc2_alpha_v <- if (is.null(bike_fit_pc2_residual$OptimalAlphaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc2_residual$OptimalAlphaV)), 4), collapse = ", ")
}

base_gamma_v <- if (is.null(bike_fit_base$OptimalGammaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_base$OptimalGammaV)), 4), collapse = ", ")
}

pc1_gamma_v <- if (is.null(bike_fit_pc1$OptimalGammaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc1$OptimalGammaV)), 4), collapse = ", ")
}

pc2_gamma_v <- if (is.null(bike_fit_pc2_residual$OptimalGammaV)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc2_residual$OptimalGammaV)), 4), collapse = ", ")
}

base_var <- if (is.null(bike_fit_base$VarianceExplained)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_base$VarianceExplained)), 4), collapse = ", ")
}

pc1_var <- if (is.null(bike_fit_pc1$VarianceExplained)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc1$VarianceExplained)), 4), collapse = ", ")
}

pc2_var <- if (is.null(bike_fit_pc2_residual$VarianceExplained)) {
  NA_character_
} else {
  paste(signif(as.numeric(unlist(bike_fit_pc2_residual$VarianceExplained)), 4), collapse = ", ")
}

fit_output_summary <- tibble(
  Object = c(
    "No penalty two-PC fit",
    "PC1 two-way fit",
    "PC2 residual two-way fit"
  ),
  OptimalAlphaU = c(
    base_alpha_u,
    pc1_alpha_u,
    pc2_alpha_u
  ),
  OptimalGammaU = c(
    base_gamma_u,
    pc1_gamma_u,
    pc2_gamma_u
  ),
  OptimalAlphaV = c(
    base_alpha_v,
    pc1_alpha_v,
    pc2_alpha_v
  ),
  OptimalGammaV = c(
    base_gamma_v,
    pc1_gamma_v,
    pc2_gamma_v
  ),
  VarianceExplained = c(
    base_var,
    pc1_var,
    pc2_var
  )
)

#print(fit_output_summary)

###############################################################################
# 12. Plot setup
###############################################################################

method_levels <- c("HPCA", "ReHPCA")

method_cols <- c(
  "HPCA" = "#D95F02",
  "ReHPCA" = "#0072B2"
)

pc_levels <- c("PC1", "PC2")

functional_block_names <- c(
  "Count",
  "Temperature",
  "Humidity",
  "Wind"
)

score_plot_df <- score_df |>
  mutate(
    Method = recode(
      Method,
      "No penalty" = "HPCA",
      "PC-specific two-way" = "ReHPCA"
    ),
    Method = factor(Method, levels = method_levels)
  )

functional_plot_df <- functional_loading_df |>
  mutate(
    Method = recode(
      as.character(Method),
      "No penalty" = "HPCA",
      "PC-specific two-way" = "ReHPCA"
    ),
    Method = factor(Method, levels = method_levels),
    PC = factor(PC, levels = pc_levels),
    Functional_Block = factor(
      Functional_Block,
      levels = functional_block_names
    )
  )

regular_plot_all <- regular_loading_df |>
  mutate(
    Method = recode(
      as.character(Method),
      "No penalty" = "HPCA",
      "PC-specific two-way" = "ReHPCA"
    ),
    Method = factor(Method, levels = method_levels),
    PC = factor(PC, levels = pc_levels),
    Feature_raw = as.character(Feature)
  )

compact_theme <- theme_bw(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey45",
      linewidth = 0.25
    ),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.22
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey40",
      linewidth = 0.30
    ),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.5, color = "black"),
    plot.margin = margin(3, 5, 3, 5)
  )

###############################################################################
# 13. HPCA versus ReHPCA comparison plots
###############################################################################

###############################################################################
# 13.1 u-direction scores
###############################################################################

score_long <- score_plot_df |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "Score"
  ) |>
  mutate(
    PC = factor(PC, levels = pc_levels)
  )

u_max <- max(abs(score_long$Score), na.rm = TRUE)

if (!is.finite(u_max) || u_max == 0) {
  u_lim <- c(-1, 1)
} else {
  u_step <- 10^floor(log10(u_max))
  u_bound <- ceiling(u_max / u_step) * u_step
  u_lim <- c(-u_bound, u_bound)
}

p_u_time <- ggplot(
  score_long,
  aes(
    x = date,
    y = Score,
    color = Method,
    group = Method
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.22,
    color = "grey45"
  ) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.90
  ) +
  facet_wrap(
    ~ PC,
    nrow = 1,
    scales = "fixed"
  ) +
  scale_color_manual(
    values = method_cols,
    breaks = method_levels
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y-%m",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = u_lim,
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  guides(
    color = guide_legend(
      nrow = 1,
      override.aes = list(linewidth = 0.9)
    )
  ) +
  labs(
    title = "u direction: PC scores over ordered days",
    x = NULL,
    y = NULL
  ) +
  compact_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.title.y.left = element_blank(),
    axis.text.y.left = element_text(margin = margin(r = 1)),
    axis.ticks.length.y = grid::unit(1.2, "pt"),
    plot.margin = margin(3, 5, 3, 0)
  )

p_u_label <- ggplot() +
  annotate(
    "text",
    x = 1,
    y = 0.5,
    label = "u score",
    angle = 90,
    size = 3.0,
    color = "black"
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(0, -8, 0, 0)
  )

###############################################################################
# 13.2 v-direction functional loading curves
###############################################################################

p_v_function <- ggplot(
  functional_plot_df,
  aes(
    x = Hour,
    y = Loading,
    color = Method,
    group = Method
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.22,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.70,
    alpha = 0.92
  ) +
  facet_grid(
    Functional_Block ~ PC,
    scales = "free_y",
    switch = "y"
  ) +
  scale_color_manual(
    values = method_cols,
    breaks = method_levels
  ) +
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 23),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  guides(color = "none") +
  labs(
    title = "v direction: functional PC loading curves",
    x = "Hour of day",
    y = "PC loading"
  ) +
  compact_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )

###############################################################################
# 13.3 v-direction regular covariate loadings
###############################################################################

regular_plot_all <- regular_plot_all |>
  mutate(
    Feature_label = case_when(
      Feature_raw == "temp_day" ~ "Daily temp.",
      Feature_raw == "atemp_day" ~ "Daily apparent temp.",
      Feature_raw == "hum_day" ~ "Daily humidity",
      Feature_raw == "windspeed_day" ~ "Daily wind speed",

      Feature_raw == "season1" ~ "Season: Spring",
      Feature_raw == "season2" ~ "Season: Summer",
      Feature_raw == "season3" ~ "Season: Fall",
      Feature_raw == "season4" ~ "Season: Winter",

      Feature_raw == "year0" ~ "Year: 2011",
      Feature_raw == "year1" ~ "Year: 2012",

      Feature_raw == "month1" ~ "Month: Jan",
      Feature_raw == "month2" ~ "Month: Feb",
      Feature_raw == "month3" ~ "Month: Mar",
      Feature_raw == "month4" ~ "Month: Apr",
      Feature_raw == "month5" ~ "Month: May",
      Feature_raw == "month6" ~ "Month: Jun",
      Feature_raw == "month7" ~ "Month: Jul",
      Feature_raw == "month8" ~ "Month: Aug",
      Feature_raw == "month9" ~ "Month: Sep",
      Feature_raw == "month10" ~ "Month: Oct",
      Feature_raw == "month11" ~ "Month: Nov",
      Feature_raw == "month12" ~ "Month: Dec",

      Feature_raw == "holiday0" ~ "Not holiday",
      Feature_raw == "holiday1" ~ "Holiday",

      Feature_raw == "weekday0" ~ "Sunday",
      Feature_raw == "weekday1" ~ "Monday",
      Feature_raw == "weekday2" ~ "Tuesday",
      Feature_raw == "weekday3" ~ "Wednesday",
      Feature_raw == "weekday4" ~ "Thursday",
      Feature_raw == "weekday5" ~ "Friday",
      Feature_raw == "weekday6" ~ "Saturday",

      Feature_raw == "workingday0" ~ "Non-working day",
      Feature_raw == "workingday1" ~ "Working day",

      Feature_raw == "weathersit1" ~ "Weather: clear",
      Feature_raw == "weathersit2" ~ "Weather: mist/cloudy",
      Feature_raw == "weathersit3" ~ "Weather: light rain/snow",
      Feature_raw == "weathersit4" ~ "Weather: heavy rain/snow",

      TRUE ~ Feature_raw
    ),
    Feature_group = case_when(
      Feature_raw %in% c(
        "temp_day",
        "atemp_day",
        "hum_day",
        "windspeed_day"
      ) ~ "Daily weather",
      str_detect(Feature_raw, "^season") ~ "Season",
      str_detect(Feature_raw, "^year") ~ "Year",
      str_detect(Feature_raw, "^month") ~ "Month",
      str_detect(Feature_raw, "^holiday") ~ "Holiday",
      str_detect(Feature_raw, "^weekday") ~ "Weekday",
      str_detect(Feature_raw, "^workingday") ~ "Working day",
      str_detect(Feature_raw, "^weathersit") ~ "Weather category",
      TRUE ~ "Other"
    )
  )

regular_feature_order_plot <- regular_plot_all |>
  group_by(
    Feature_raw,
    Feature_label,
    Feature_group
  ) |>
  summarize(
    MaxAbsLoading = max(abs(Loading), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(
    Feature_group,
    desc(MaxAbsLoading)
  )

regular_plot_df <- regular_plot_all |>
  mutate(
    Feature_label = factor(
      Feature_label,
      levels = rev(
        regular_feature_order_plot$Feature_label
      )
    )
  )

p_v_regular <- ggplot(
  regular_plot_df,
  aes(
    x = Loading,
    y = Feature_label,
    fill = Method
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.22,
    color = "grey40"
  ) +
  geom_col(
    position = position_dodge(width = 0.68),
    width = 0.58,
    alpha = 0.92
  ) +
  facet_grid(
    . ~ PC,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = method_cols,
    breaks = method_levels
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 3),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  guides(fill = "none") +
  labs(
    title = "v direction: regular covariate loadings",
    x = "PC loading",
    y = NULL
  ) +
  compact_theme +
  theme(
    axis.text.y = element_text(size = 6.3),
    axis.text.x = element_text(size = 7),
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.20
    ),
    legend.position = "none"
  )

###############################################################################
# 13.4 Combined layout
###############################################################################

p_u_top <- (p_u_label | p_u_time) +
  plot_layout(
    widths = c(0.035, 1)
  )

p_bottom_v <- (p_v_function | p_v_regular) +
  plot_layout(
    widths = c(1.65, 1.00),
    guides = "collect"
  )

p_bike_combined <- (p_u_top / p_bottom_v) +
  plot_layout(
    heights = c(0.55, 1.80),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

###############################################################################
# 13.5 Display all comparison plots
###############################################################################

print(p_u_time)
print(p_v_function)
print(p_v_regular)
print(p_bike_combined)

###############################################################################
# 14. Mean curves plus/minus ReHPCA functional PCs
###############################################################################


functional_order <- c(
  "Count",
  "Temperature",
  "Humidity",
  "Wind"
)

pc_order <- c(
  "PC1",
  "PC2"
)

# Raw-scale means.
mean_count <- colMeans(X_count_raw, na.rm = TRUE)
mean_temp  <- colMeans(X_temp_raw, na.rm = TRUE)
mean_hum   <- colMeans(X_hum_raw, na.rm = TRUE)
mean_wind  <- colMeans(X_wind_raw, na.rm = TRUE)

# Raw-scale standard deviations used to undo the global functional scaling.
sd_count_raw <- sd(as.vector(X_count_raw), na.rm = TRUE)
sd_temp_raw  <- sd(as.vector(X_temp_raw), na.rm = TRUE)
sd_hum_raw   <- sd(as.vector(X_hum_raw), na.rm = TRUE)
sd_wind_raw  <- sd(as.vector(X_wind_raw), na.rm = TRUE)

if (!is.finite(sd_count_raw) || sd_count_raw == 0) sd_count_raw <- 1
if (!is.finite(sd_temp_raw)  || sd_temp_raw  == 0) sd_temp_raw  <- 1
if (!is.finite(sd_hum_raw)   || sd_hum_raw   == 0) sd_hum_raw   <- 1
if (!is.finite(sd_wind_raw)  || sd_wind_raw  == 0) sd_wind_raw  <- 1

# Score scales.
score_sd_pc1 <- sd(
  as.numeric(bike_fit_pc_specific$PCScores[, 1]),
  na.rm = TRUE
)

score_sd_pc2 <- sd(
  as.numeric(bike_fit_pc_specific$PCScores[, 2]),
  na.rm = TRUE
)

if (!is.finite(score_sd_pc1) || score_sd_pc1 == 0) score_sd_pc1 <- 1
if (!is.finite(score_sd_pc2) || score_sd_pc2 == 0) score_sd_pc2 <- 1

# Raw-scale PC deviations.
dev_pc1_count <- score_sd_pc1 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[1]]) *
  sd_count_raw

dev_pc1_temp <- score_sd_pc1 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[2]]) *
  sd_temp_raw

dev_pc1_hum <- score_sd_pc1 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[3]]) *
  sd_hum_raw

dev_pc1_wind <- score_sd_pc1 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[1]][[4]]) *
  sd_wind_raw

dev_pc2_count <- score_sd_pc2 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[1]]) *
  sd_count_raw

dev_pc2_temp <- score_sd_pc2 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[2]]) *
  sd_temp_raw

dev_pc2_hum <- score_sd_pc2 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[3]]) *
  sd_hum_raw

dev_pc2_wind <- score_sd_pc2 *
  as.numeric(bike_fit_pc_specific$PCFunctions[[2]][[4]]) *
  sd_wind_raw

mean_pm_pc_df <- bind_rows(
  tibble(
    PC = "PC1",
    Functional_Block = "Count",
    Hour = hour_grid,
    Mean = mean_count,
    Mean_minus_PC = mean_count - dev_pc1_count,
    Mean_plus_PC = mean_count + dev_pc1_count
  ),
  tibble(
    PC = "PC1",
    Functional_Block = "Temperature",
    Hour = hour_grid,
    Mean = mean_temp,
    Mean_minus_PC = mean_temp - dev_pc1_temp,
    Mean_plus_PC = mean_temp + dev_pc1_temp
  ),
  tibble(
    PC = "PC1",
    Functional_Block = "Humidity",
    Hour = hour_grid,
    Mean = mean_hum,
    Mean_minus_PC = mean_hum - dev_pc1_hum,
    Mean_plus_PC = mean_hum + dev_pc1_hum
  ),
  tibble(
    PC = "PC1",
    Functional_Block = "Wind",
    Hour = hour_grid,
    Mean = mean_wind,
    Mean_minus_PC = mean_wind - dev_pc1_wind,
    Mean_plus_PC = mean_wind + dev_pc1_wind
  ),
  tibble(
    PC = "PC2",
    Functional_Block = "Count",
    Hour = hour_grid,
    Mean = mean_count,
    Mean_minus_PC = mean_count - dev_pc2_count,
    Mean_plus_PC = mean_count + dev_pc2_count
  ),
  tibble(
    PC = "PC2",
    Functional_Block = "Temperature",
    Hour = hour_grid,
    Mean = mean_temp,
    Mean_minus_PC = mean_temp - dev_pc2_temp,
    Mean_plus_PC = mean_temp + dev_pc2_temp
  ),
  tibble(
    PC = "PC2",
    Functional_Block = "Humidity",
    Hour = hour_grid,
    Mean = mean_hum,
    Mean_minus_PC = mean_hum - dev_pc2_hum,
    Mean_plus_PC = mean_hum + dev_pc2_hum
  ),
  tibble(
    PC = "PC2",
    Functional_Block = "Wind",
    Hour = hour_grid,
    Mean = mean_wind,
    Mean_minus_PC = mean_wind - dev_pc2_wind,
    Mean_plus_PC = mean_wind + dev_pc2_wind
  )
) |>
  mutate(
    PC = factor(PC, levels = pc_order),
    Functional_Block = factor(
      Functional_Block,
      levels = functional_order
    )
  )

mean_line_df <- mean_pm_pc_df |>
  transmute(
    PC,
    Functional_Block,
    Hour,
    Curve = "Mean",
    Value = Mean
  ) |>
  mutate(
    Curve = factor(
      Curve,
      levels = c(
        "Mean - PC",
        "Mean",
        "Mean + PC"
      )
    )
  )

pm_symbol_df_raw <- bind_rows(
  mean_pm_pc_df |>
    transmute(
      PC,
      Functional_Block,
      Hour,
      Curve = "Mean - PC",
      Value = Mean_minus_PC,
      Symbol = "\u2212"
    ),
  mean_pm_pc_df |>
    transmute(
      PC,
      Functional_Block,
      Hour,
      Curve = "Mean + PC",
      Value = Mean_plus_PC,
      Symbol = "+"
    )
)

# Dense interpolation for the + and - symbol-lines.
dense_hour_grid <- seq(
  0,
  23,
  length.out = 200
)

pm_symbol_df <- tibble()

for (pc_i in pc_order) {

  for (block_i in functional_order) {

    for (curve_i in c("Mean - PC", "Mean + PC")) {

      tmp <- pm_symbol_df_raw |>
        filter(
          PC == pc_i,
          Functional_Block == block_i,
          Curve == curve_i
        )

      tmp_interp <- approx(
        x = tmp$Hour,
        y = tmp$Value,
        xout = dense_hour_grid,
        rule = 2
      )

      pm_symbol_df <- bind_rows(
        pm_symbol_df,
        tibble(
          PC = pc_i,
          Functional_Block = block_i,
          Hour = dense_hour_grid,
          Curve = curve_i,
          Value = tmp_interp$y,
          Symbol = if_else(
            curve_i == "Mean - PC",
            "\u2212",
            "+"
          )
        )
      )
    }
  }
}

pm_symbol_df <- pm_symbol_df |>
  mutate(
    PC = factor(PC, levels = pc_order),
    Functional_Block = factor(
      Functional_Block,
      levels = functional_order
    ),
    Curve = factor(
      Curve,
      levels = c(
        "Mean - PC",
        "Mean",
        "Mean + PC"
      )
    )
  )

legend_line_df <- tibble(
  Hour = c(0, 1, 0, 1, 0, 1),
  Value = c(0, 0, 0, 0, 0, 0),
  Curve = factor(
    rep(
      c(
        "Mean - PC",
        "Mean",
        "Mean + PC"
      ),
      each = 2
    ),
    levels = c(
      "Mean - PC",
      "Mean",
      "Mean + PC"
    )
  ),
  PC = factor(
    "PC1",
    levels = pc_order
  ),
  Functional_Block = factor(
    "Count",
    levels = functional_order
  ),
  legend_group = rep(
    c(
      "Mean - PC",
      "Mean",
      "Mean + PC"
    ),
    each = 2
  )
)

pm_cols <- c(
  "Mean - PC" = "#D62728",
  "Mean" = "black",
  "Mean + PC" = "#009E73"
)

pm_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey45",
      linewidth = 0.25
    ),
    strip.placement = "outside",
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "grey40",
      linewidth = 0.30
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(1.15, "cm"),
    legend.spacing.x = grid::unit(0.35, "cm"),
    plot.margin = margin(6, 8, 6, 8)
  )

p_mean_pm_rehpca <- ggplot() +

  geom_line(
    data = legend_line_df,
    aes(
      x = Hour,
      y = Value,
      color = Curve,
      group = legend_group
    ),
    linewidth = 1.35,
    alpha = 0,
    show.legend = TRUE
  ) +

  geom_line(
    data = mean_line_df,
    aes(
      x = Hour,
      y = Value,
      color = Curve
    ),
    linewidth = 1.10,
    lineend = "round",
    show.legend = FALSE
  ) +

  geom_text(
    data = pm_symbol_df,
    aes(
      x = Hour,
      y = Value,
      color = Curve,
      label = Symbol
    ),
    size = 2.8,
    fontface = "plain",
    show.legend = FALSE
  ) +

  facet_wrap(
    vars(PC, Functional_Block),
    nrow = 2,
    scales = "free_y",
    labeller = labeller(
      PC = label_value,
      Functional_Block = label_value
    )
  ) +

  scale_color_manual(
    values = pm_cols,
    breaks = c(
      "Mean - PC",
      "Mean",
      "Mean + PC"
    ),
    drop = FALSE
  ) +

  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 23),
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 4),
    expand = expansion(mult = c(0.10, 0.10))
  ) +

  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(
        alpha = 1,
        linewidth = 1.5,
        linetype = "solid"
      )
    )
  ) +

  labs(
    title = "Mean curves plus and minus ReHPCA functional PCs",
    x = "Hour of day",
    y = NULL
  ) +

  pm_theme

print(p_mean_pm_rehpca)

###############################################################################
# 15. Daily-hour fitted heatmaps: HPCA versus ReHPCA
###############################################################################

# Reconstruct the fitted Count block for HPCA.
rec_hpca_count_scaled <- matrix(
  0,
  nrow = n_days,
  ncol = ncol(X_count)
)

for (pc in seq_len(2)) {

  u_pc <- as.numeric(
    bike_fit_base$PCScores[, pc]
  )

  v_count_pc <- as.numeric(
    bike_fit_base$PCFunctions[[pc]][[1]]
  )

  rec_hpca_count_scaled <-
    rec_hpca_count_scaled +
    outer(
      u_pc,
      v_count_pc
    )
}

# Reconstruct the fitted Count block for ReHPCA.
rec_rehpca_count_scaled <- matrix(
  0,
  nrow = n_days,
  ncol = ncol(X_count)
)

for (pc in seq_len(2)) {

  u_pc <- as.numeric(
    bike_fit_pc_specific$PCScores[, pc]
  )

  v_count_pc <- as.numeric(
    bike_fit_pc_specific$PCFunctions[[pc]][[1]]
  )

  rec_rehpca_count_scaled <-
    rec_rehpca_count_scaled +
    outer(
      u_pc,
      v_count_pc
    )
}

count_raw_mean <- mean(
  as.vector(X_count_raw),
  na.rm = TRUE
)

count_raw_sd <- sd(
  as.vector(X_count_raw),
  na.rm = TRUE
)

if (!is.finite(count_raw_sd) || count_raw_sd == 0) {
  count_raw_sd <- 1
}

rec_hpca_count_raw <- pmax(
  rec_hpca_count_scaled *
    count_raw_sd +
    count_raw_mean,
  0
)

rec_rehpca_count_raw <- pmax(
  rec_rehpca_count_scaled *
    count_raw_sd +
    count_raw_mean,
  0
)

# Convert HPCA fitted matrix to long format.
colnames(rec_hpca_count_raw) <-
  as.character(hour_grid)

hpca_heatmap_df <- as_tibble(
  rec_hpca_count_raw,
  .name_repair = "minimal"
) |>
  mutate(
    date = day_complete$dteday
  ) |>
  pivot_longer(
    cols = -date,
    names_to = "hour",
    values_to = "count"
  ) |>
  mutate(
    hour = as.numeric(hour),
    Method = "HPCA fitted",
    count = pmax(count, 0),
    log_count = log1p(count)
  )

# Convert ReHPCA fitted matrix to long format.
colnames(rec_rehpca_count_raw) <-
  as.character(hour_grid)

rehpca_heatmap_df <- as_tibble(
  rec_rehpca_count_raw,
  .name_repair = "minimal"
) |>
  mutate(
    date = day_complete$dteday
  ) |>
  pivot_longer(
    cols = -date,
    names_to = "hour",
    values_to = "count"
  ) |>
  mutate(
    hour = as.numeric(hour),
    Method = "ReHPCA fitted",
    count = pmax(count, 0),
    log_count = log1p(count)
  )

heatmap_fit_df <- bind_rows(
  hpca_heatmap_df,
  rehpca_heatmap_df
) |>
  mutate(
    Method = factor(
      Method,
      levels = c(
        "HPCA fitted",
        "ReHPCA fitted"
      )
    )
  )

my_heat_palette <- colorRampPalette(c(
  "#2D004B",
  "#542788",
  "#8073AC",
  "#B358A0",
  "#D65F5F",
  "#F07C3E",
  "#F6B44B",
  "#F4D35E",
  "#FFF3A3"
))(256)

fill_upper <- quantile(
  heatmap_fit_df$log_count,
  probs = 0.99,
  na.rm = TRUE
)

if (!is.finite(fill_upper) || fill_upper <= 0) {
  fill_upper <- max(
    heatmap_fit_df$log_count,
    na.rm = TRUE
  )
}

fill_limits <- c(
  0,
  fill_upper
)

heatmap_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey55",
      linewidth = 0.25
    ),
    axis.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      size = 8,
      color = "black"
    ),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(
      0.18,
      "cm"
    ),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.width = grid::unit(
      1.7,
      "cm"
    ),
    legend.key.height = grid::unit(
      0.23,
      "cm"
    ),
    plot.margin = margin(5, 8, 5, 8)
  )

p_fitted_heatmaps <- ggplot(
  heatmap_fit_df,
  aes(
    x = hour,
    y = date,
    fill = log_count
  )
) +
  geom_tile(
    width = 1,
    height = 1
  ) +
  facet_wrap(
    ~ Method,
    nrow = 1
  ) +
  scale_fill_gradientn(
    colors = my_heat_palette,
    limits = fill_limits,
    oob = scales::squish,
    name = expression(log(1 + count))
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      23,
      by = 4
    ),
    expand = c(0, 0)
  ) +
  scale_y_date(
    date_breaks = "3 months",
    date_labels = "%Y-%m",
    expand = expansion(
      mult = c(0.002, 0.002)
    )
  ) +
  labs(
    title = "Daily-hour fitted heatmaps of bike-rental demand",
    x = "Hour of day",
    y = "Date"
  ) +
  heatmap_theme

print(p_fitted_heatmaps)

###############################################################################
# 16. Built-in ReMPCA plots
###############################################################################

try(
  plot_pc_functions(
    bike_fit_base
  ),
  silent = TRUE
)

try(
  plot_pc_scores(
    bike_fit_base
  ),
  silent = TRUE
)

try(
  plot_pc_functions(
    bike_fit_pc1
  ),
  silent = TRUE
)

try(
  plot_pc_scores(
    bike_fit_pc1
  ),
  silent = TRUE
)

try(
  plot_pc_functions(
    bike_fit_pc2_residual
  ),
  silent = TRUE
)

try(
  plot_pc_scores(
    bike_fit_pc2_residual
  ),
  silent = TRUE
)
