###############################################################################
# UCI Appliances Energy Prediction Dataset
###############################################################################

###############################################################################
# 0. Load packages
###############################################################################

library(ReMPCA)
library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)

###############################################################################
# 1. Load the UCI Appliances Energy Prediction data
###############################################################################

load(file.path("data", "energydata.rda"))

appliance_dat$date <- ymd_hms(appliance_dat$date)

cat("Raw data dimensions:\n")
print(dim(appliance_dat))

###############################################################################
# 2. Prepare the daily data and derived functional variables
###############################################################################

indoor_temp_cols <- paste0("T", c(1, 2, 3, 4, 5, 7, 8, 9))
indoor_hum_cols  <- paste0("RH_", c(1, 2, 3, 4, 5, 7, 8, 9))

appliance_dat <- appliance_dat |>
  mutate(
    day = as.Date(date),
    slot = hour(date) * 6 + minute(date) / 10,
    hour_of_day = slot / 6,
    Indoor_temperature = rowMeans(
      as.matrix(across(all_of(indoor_temp_cols))),
      na.rm = TRUE
    ),
    Indoor_humidity = rowMeans(
      as.matrix(across(all_of(indoor_hum_cols))),
      na.rm = TRUE
    ),
    Outdoor_temperature = T_out,
    Outdoor_humidity = RH_out
  )

# Keep only complete calendar days with all 144 ten-minute measurements.
complete_days <- appliance_dat |>
  group_by(day) |>
  summarize(n_slots = n_distinct(slot), .groups = "drop") |>
  filter(n_slots == 144) |>
  pull(day)

appliance_complete <- appliance_dat |>
  filter(day %in% complete_days) |>
  arrange(day, slot)

day_complete <- appliance_complete |>
  distinct(day) |>
  arrange(day)

n_days <- nrow(day_complete)
hour_grid <- sort(unique(appliance_complete$hour_of_day))
day_grid <- seq(0, 1, length.out = n_days)

cat("\nComplete days:\n")
print(n_days)

cat("\nTime points per day:\n")
print(length(hour_grid))

###############################################################################
# 3. Build daily functional curves from the ten-minute data
###############################################################################

# Appliance energy use.
X_appliance_raw <- appliance_complete |>
  select(day, slot, value = Appliances) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Lights energy use.
X_lights_raw <- appliance_complete |>
  select(day, slot, value = lights) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Indoor temperature.
X_tin_raw <- appliance_complete |>
  select(day, slot, value = Indoor_temperature) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Indoor humidity.
X_hin_raw <- appliance_complete |>
  select(day, slot, value = Indoor_humidity) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Outdoor temperature.
X_tout_raw <- appliance_complete |>
  select(day, slot, value = Outdoor_temperature) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Outdoor humidity.
X_hout_raw <- appliance_complete |>
  select(day, slot, value = Outdoor_humidity) |>
  arrange(day, slot) |>
  mutate(slot = paste0("s", sprintf("%03d", slot))) |>
  pivot_wider(names_from = slot, values_from = value) |>
  arrange(day) |>
  select(-day) |>
  as.matrix()

# Preserve the calendar-day labels used in the original script.
rownames(X_appliance_raw) <- as.character(day_complete$day)
rownames(X_lights_raw)    <- as.character(day_complete$day)
rownames(X_tin_raw)       <- as.character(day_complete$day)
rownames(X_hin_raw)       <- as.character(day_complete$day)
rownames(X_tout_raw)      <- as.character(day_complete$day)
rownames(X_hout_raw)      <- as.character(day_complete$day)

cat("\nFunctional block dimensions:\n")
print(dim(X_appliance_raw))
print(dim(X_lights_raw))
print(dim(X_tin_raw))
print(dim(X_hin_raw))
print(dim(X_tout_raw))
print(dim(X_hout_raw))

###############################################################################
# 4. Build the daily non-functional covariate matrix
###############################################################################

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
  mutate(
    month_factor = factor(
      month_num,
      levels = sort(unique(month_num))
    ),
    Day_type = factor(
      Day_type,
      levels = c("Weekday", "Weekend")
    )
  )

# Month dummy variables for the months observed in the data.
X_month <- model.matrix(~ month_factor - 1, data = regular_df)
colnames(X_month) <- paste0("month_", levels(regular_df$month_factor))

X_daytype <- regular_df |>
  transmute(
    daytype_weekday = as.integer(Day_type == "Weekday"),
    daytype_weekend = as.integer(Day_type == "Weekend")
  ) |>
  as.matrix()

# Daily weather variables.
X_weather <- regular_df |>
  transmute(
    Pressure_day,
    Windspeed_day,
    Visibility_day,
    Tdewpoint_day
  ) |>
  as.matrix()

X_regular_raw <- cbind(
  X_weather,
  X_month,
  X_daytype
)

regular_feature_names <- colnames(X_regular_raw)

cat("\nRegular block dimensions:\n")
print(dim(X_regular_raw))

cat("\nRegular feature names:\n")
print(regular_feature_names)

###############################################################################
# 5. Scale the functional and non-functional blocks
###############################################################################

# Global scaling for the functional blocks.
X_appliance <- as.matrix(X_appliance_raw)
appliance_sd <- sd(as.vector(X_appliance), na.rm = TRUE)
if (!is.finite(appliance_sd) || appliance_sd == 0) appliance_sd <- 1
X_appliance <- (X_appliance - mean(X_appliance, na.rm = TRUE)) / appliance_sd

X_lights <- as.matrix(X_lights_raw)
lights_sd <- sd(as.vector(X_lights), na.rm = TRUE)
if (!is.finite(lights_sd) || lights_sd == 0) lights_sd <- 1
X_lights <- (X_lights - mean(X_lights, na.rm = TRUE)) / lights_sd

X_tin <- as.matrix(X_tin_raw)
tin_sd <- sd(as.vector(X_tin), na.rm = TRUE)
if (!is.finite(tin_sd) || tin_sd == 0) tin_sd <- 1
X_tin <- (X_tin - mean(X_tin, na.rm = TRUE)) / tin_sd

X_hin <- as.matrix(X_hin_raw)
hin_sd <- sd(as.vector(X_hin), na.rm = TRUE)
if (!is.finite(hin_sd) || hin_sd == 0) hin_sd <- 1
X_hin <- (X_hin - mean(X_hin, na.rm = TRUE)) / hin_sd

X_tout <- as.matrix(X_tout_raw)
tout_sd <- sd(as.vector(X_tout), na.rm = TRUE)
if (!is.finite(tout_sd) || tout_sd == 0) tout_sd <- 1
X_tout <- (X_tout - mean(X_tout, na.rm = TRUE)) / tout_sd

X_hout <- as.matrix(X_hout_raw)
hout_sd <- sd(as.vector(X_hout), na.rm = TRUE)
if (!is.finite(hout_sd) || hout_sd == 0) hout_sd <- 1
X_hout <- (X_hout - mean(X_hout, na.rm = TRUE)) / hout_sd

# Columnwise scaling for the non-functional block.
X_regular <- scale(
  as.matrix(X_regular_raw),
  center = TRUE,
  scale = TRUE
)
X_regular[is.na(X_regular)] <- 0
X_regular <- as.matrix(X_regular)

blocks_original <- list(
  Appliance_curve = X_appliance,
  Lights_curve = X_lights,
  Indoor_temperature_curve = X_tin,
  Indoor_humidity_curve = X_hin,
  Outdoor_temperature_curve = X_tout,
  Outdoor_humidity_curve = X_hout,
  Daily_covariates = X_regular
)

NUM_PCS <- 2

###############################################################################
# 6. Fit the strict no-penalty HPCA baseline with two PCs
###############################################################################

fd_appliance_base <- fdClass(
  data = blocks_original$Appliance_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_lights_base <- fdClass(
  data = blocks_original$Lights_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_tin_base <- fdClass(
  data = blocks_original$Indoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_hin_base <- fdClass(
  data = blocks_original$Indoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_tout_base <- fdClass(
  data = blocks_original$Outdoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

fd_hout_base <- fdClass(
  data = blocks_original$Outdoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

rd_daily_base <- rdClass(
  data = blocks_original$Daily_covariates,
  Sparsity_parameter = 0L
)

appliances_hd_base <- hdClass(
  hdlist = list(
    Appliance_curve = fd_appliance_base,
    Lights_curve = fd_lights_base,
    Indoor_temperature_curve = fd_tin_base,
    Indoor_humidity_curve = fd_hin_base,
    Outdoor_temperature_curve = fd_tout_base,
    Outdoor_humidity_curve = fd_hout_base,
    Daily_covariates = rd_daily_base
  ),
  argval = day_grid,
  Smoothing_parameter = 0,
  Sparsity_parameter = 0L
)

set.seed(2026)

cat("\nFitting HPCA: strict no-penalty baseline with two PCs...\n")

fit_hpca <- ReMPCA(
  hd = appliances_hd_base,
  centerhds = FALSE,
  num_pcs = NUM_PCS,
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
  sparse_tuning_v = list(0L, 0L, 0L, 0L, 0L, 0L, 0L),
  smooth_tuning_v = list(0, 0, 0, 0, 0, 0, 0)
)

###############################################################################
# 7. Fit ReHPCA PC1 on the original data
###############################################################################

fd_appliance_pc1 <- fdClass(
  data = blocks_original$Appliance_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^0,
  Sparsity_parameter = 8L
)

fd_lights_pc1 <- fdClass(
  data = blocks_original$Lights_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^1,
  Sparsity_parameter = 50L
)

fd_tin_pc1 <- fdClass(
  data = blocks_original$Indoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-3,
  Sparsity_parameter = 0L
)

fd_hin_pc1 <- fdClass(
  data = blocks_original$Indoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-3,
  Sparsity_parameter = 0L
)

fd_tout_pc1 <- fdClass(
  data = blocks_original$Outdoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-3,
  Sparsity_parameter = 0L
)

fd_hout_pc1 <- fdClass(
  data = blocks_original$Outdoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-3,
  Sparsity_parameter = 0L
)

rd_daily_pc1 <- rdClass(
  data = blocks_original$Daily_covariates,
  Sparsity_parameter = 3L
)

appliances_hd_pc1 <- hdClass(
  hdlist = list(
    Appliance_curve = fd_appliance_pc1,
    Lights_curve = fd_lights_pc1,
    Indoor_temperature_curve = fd_tin_pc1,
    Indoor_humidity_curve = fd_hin_pc1,
    Outdoor_temperature_curve = fd_tout_pc1,
    Outdoor_humidity_curve = fd_hout_pc1,
    Daily_covariates = rd_daily_pc1
  ),
  argval = day_grid,
  Smoothing_parameter = 2^-12,
  Sparsity_parameter = 45L
)

set.seed(2026)

cat("\nFitting ReHPCA PC1 on the original data...\n")

fit_rehpca_pc1 <- ReMPCA(
  hd = appliances_hd_pc1,
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
  sparse_tuning_u = 45L,
  smooth_tuning_u = 2^-12,
  sparse_tuning_v = list(
    8L,
    50L,
    0L,
    0L,
    0L,
    0L,
    3L
  ),
  smooth_tuning_v = list(
    2^0,
    2^1,
    2^-3,
    2^-3,
    2^-3,
    2^-3,
    0
  )
)

###############################################################################
# 8. Deflate the original data using ReHPCA PC1
###############################################################################

u_pc1 <- as.numeric(fit_rehpca_pc1$PCScores[, 1])

pc1_reconstruction <- list(
  Appliance_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[1]])
  ),
  Lights_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[2]])
  ),
  Indoor_temperature_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[3]])
  ),
  Indoor_humidity_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[4]])
  ),
  Outdoor_temperature_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[5]])
  ),
  Outdoor_humidity_curve = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[6]])
  ),
  Daily_covariates = outer(
    u_pc1,
    as.numeric(fit_rehpca_pc1$PCFunctions[[1]][[7]])
  )
)

blocks_residual_after_pc1 <- list(
  Appliance_curve =
    blocks_original$Appliance_curve - pc1_reconstruction$Appliance_curve,
  Lights_curve =
    blocks_original$Lights_curve - pc1_reconstruction$Lights_curve,
  Indoor_temperature_curve =
    blocks_original$Indoor_temperature_curve - pc1_reconstruction$Indoor_temperature_curve,
  Indoor_humidity_curve =
    blocks_original$Indoor_humidity_curve - pc1_reconstruction$Indoor_humidity_curve,
  Outdoor_temperature_curve =
    blocks_original$Outdoor_temperature_curve - pc1_reconstruction$Outdoor_temperature_curve,
  Outdoor_humidity_curve =
    blocks_original$Outdoor_humidity_curve - pc1_reconstruction$Outdoor_humidity_curve,
  Daily_covariates =
    blocks_original$Daily_covariates - pc1_reconstruction$Daily_covariates
)

###############################################################################
# 9. Fit ReHPCA PC2
###############################################################################

fd_appliance_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Appliance_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^2,
  Sparsity_parameter = 20L
)

fd_lights_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Lights_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^2,
  Sparsity_parameter = 60L
)

fd_tin_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Indoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-2,
  Sparsity_parameter = 0L
)

fd_hin_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Indoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-2,
  Sparsity_parameter = 0L
)

fd_tout_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Outdoor_temperature_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-2,
  Sparsity_parameter = 0L
)

fd_hout_pc2 <- fdClass(
  data = blocks_residual_after_pc1$Outdoor_humidity_curve,
  argval = hour_grid,
  Smoothing_parameter = 2^-2,
  Sparsity_parameter = 0L
)

rd_daily_pc2 <- rdClass(
  data = blocks_residual_after_pc1$Daily_covariates,
  Sparsity_parameter = 5L
)

appliances_hd_pc2 <- hdClass(
  hdlist = list(
    Appliance_curve = fd_appliance_pc2,
    Lights_curve = fd_lights_pc2,
    Indoor_temperature_curve = fd_tin_pc2,
    Indoor_humidity_curve = fd_hin_pc2,
    Outdoor_temperature_curve = fd_tout_pc2,
    Outdoor_humidity_curve = fd_hout_pc2,
    Daily_covariates = rd_daily_pc2
  ),
  argval = day_grid,
  Smoothing_parameter = 2^-10,
  Sparsity_parameter = 55L
)

set.seed(2026)

cat("\nFitting ReHPCA PC2 on the residual data...\n")

fit_rehpca_pc2 <- ReMPCA(
  hd = appliances_hd_pc2,
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
  sparse_tuning_u = 55L,
  smooth_tuning_u = 2^-10,
  sparse_tuning_v = list(
    20L,
    60L,
    0L,
    0L,
    0L,
    0L,
    5L
  ),
  smooth_tuning_v = list(
    2^2,
    2^2,
    2^-2,
    2^-2,
    2^-2,
    2^-2,
    0
  )
)

###############################################################################
# 10.  Align signs with HPCA
###############################################################################

fit_rehpca <- list(
  PCScores = cbind(
    PC1 = as.numeric(fit_rehpca_pc1$PCScores[, 1]),
    PC2 = as.numeric(fit_rehpca_pc2$PCScores[, 1])
  ),
  PCFunctions = list(
    fit_rehpca_pc1$PCFunctions[[1]],
    fit_rehpca_pc2$PCFunctions[[1]]
  )
)

for (pc in seq_len(NUM_PCS)) {

  cc <- suppressWarnings(
    cor(
      as.numeric(fit_hpca$PCScores[, pc]),
      as.numeric(fit_rehpca$PCScores[, pc]),
      use = "complete.obs"
    )
  )

  if (is.finite(cc) && cc < 0) {

    fit_rehpca$PCScores[, pc] <- -fit_rehpca$PCScores[, pc]

    for (b in seq_along(fit_rehpca$PCFunctions[[pc]])) {
      fit_rehpca$PCFunctions[[pc]][[b]] <-
        -fit_rehpca$PCFunctions[[pc]][[b]]
    }
  }
}

models <- list(
  "HPCA" = fit_hpca,
  "ReHPCA" = fit_rehpca
)

###############################################################################
# 11. Extract scores and loadings
###############################################################################

# u-direction scores.
score_hpca <- as.data.frame(
  fit_hpca$PCScores[, seq_len(NUM_PCS), drop = FALSE]
)
colnames(score_hpca) <- paste0("PC", seq_len(NUM_PCS))

score_rehpca <- as.data.frame(
  fit_rehpca$PCScores[, seq_len(NUM_PCS), drop = FALSE]
)
colnames(score_rehpca) <- paste0("PC", seq_len(NUM_PCS))

score_df_hpca <- bind_cols(
  tibble(
    date = day_complete$day,
    Method = "HPCA"
  ),
  score_hpca
)

score_df_rehpca <- bind_cols(
  tibble(
    date = day_complete$day,
    Method = "ReHPCA"
  ),
  score_rehpca
)

score_df <- bind_rows(
  score_df_hpca,
  score_df_rehpca
)

# Functional v-direction loadings.
functional_loading_hpca_pc1 <- tibble(
  Method = "HPCA",
  PC = "PC1",
  Hour = hour_grid,
  Appliances = as.numeric(fit_hpca$PCFunctions[[1]][[1]]),
  Lights = as.numeric(fit_hpca$PCFunctions[[1]][[2]]),
  `Indoor temp.` = as.numeric(fit_hpca$PCFunctions[[1]][[3]]),
  `Indoor humidity` = as.numeric(fit_hpca$PCFunctions[[1]][[4]]),
  `Outdoor temp.` = as.numeric(fit_hpca$PCFunctions[[1]][[5]]),
  `Outdoor humidity` = as.numeric(fit_hpca$PCFunctions[[1]][[6]])
)

functional_loading_hpca_pc2 <- tibble(
  Method = "HPCA",
  PC = "PC2",
  Hour = hour_grid,
  Appliances = as.numeric(fit_hpca$PCFunctions[[2]][[1]]),
  Lights = as.numeric(fit_hpca$PCFunctions[[2]][[2]]),
  `Indoor temp.` = as.numeric(fit_hpca$PCFunctions[[2]][[3]]),
  `Indoor humidity` = as.numeric(fit_hpca$PCFunctions[[2]][[4]]),
  `Outdoor temp.` = as.numeric(fit_hpca$PCFunctions[[2]][[5]]),
  `Outdoor humidity` = as.numeric(fit_hpca$PCFunctions[[2]][[6]])
)

functional_loading_rehpca_pc1 <- tibble(
  Method = "ReHPCA",
  PC = "PC1",
  Hour = hour_grid,
  Appliances = as.numeric(fit_rehpca$PCFunctions[[1]][[1]]),
  Lights = as.numeric(fit_rehpca$PCFunctions[[1]][[2]]),
  `Indoor temp.` = as.numeric(fit_rehpca$PCFunctions[[1]][[3]]),
  `Indoor humidity` = as.numeric(fit_rehpca$PCFunctions[[1]][[4]]),
  `Outdoor temp.` = as.numeric(fit_rehpca$PCFunctions[[1]][[5]]),
  `Outdoor humidity` = as.numeric(fit_rehpca$PCFunctions[[1]][[6]])
)

functional_loading_rehpca_pc2 <- tibble(
  Method = "ReHPCA",
  PC = "PC2",
  Hour = hour_grid,
  Appliances = as.numeric(fit_rehpca$PCFunctions[[2]][[1]]),
  Lights = as.numeric(fit_rehpca$PCFunctions[[2]][[2]]),
  `Indoor temp.` = as.numeric(fit_rehpca$PCFunctions[[2]][[3]]),
  `Indoor humidity` = as.numeric(fit_rehpca$PCFunctions[[2]][[4]]),
  `Outdoor temp.` = as.numeric(fit_rehpca$PCFunctions[[2]][[5]]),
  `Outdoor humidity` = as.numeric(fit_rehpca$PCFunctions[[2]][[6]])
)

functional_block_names <- c(
  "Appliances",
  "Lights",
  "Indoor temp.",
  "Indoor humidity",
  "Outdoor temp.",
  "Outdoor humidity"
)

functional_loading_df <- bind_rows(
  functional_loading_hpca_pc1,
  functional_loading_hpca_pc2,
  functional_loading_rehpca_pc1,
  functional_loading_rehpca_pc2
) |>
  pivot_longer(
    cols = all_of(functional_block_names),
    names_to = "Functional_Block",
    values_to = "Loading"
  )

# Non-functional v-direction loadings.
regular_loading_df <- bind_rows(
  tibble(
    Method = "HPCA",
    PC = "PC1",
    Feature = regular_feature_names,
    Loading = as.numeric(fit_hpca$PCFunctions[[1]][[7]])
  ),
  tibble(
    Method = "HPCA",
    PC = "PC2",
    Feature = regular_feature_names,
    Loading = as.numeric(fit_hpca$PCFunctions[[2]][[7]])
  ),
  tibble(
    Method = "ReHPCA",
    PC = "PC1",
    Feature = regular_feature_names,
    Loading = as.numeric(fit_rehpca$PCFunctions[[1]][[7]])
  ),
  tibble(
    Method = "ReHPCA",
    PC = "PC2",
    Feature = regular_feature_names,
    Loading = as.numeric(fit_rehpca$PCFunctions[[2]][[7]])
  )
)


###############################################################################
# 13. Separated combined plot
#     Top row: u direction
#     Bottom row: v direction, functional + regular covariates
###############################################################################

method_levels <- c("HPCA", "ReHPCA")

method_cols <- c(
  "HPCA" = "#D95F02",
  "ReHPCA" = "#0072B2"
)

pc_levels <- paste0("PC", seq_len(NUM_PCS))

compact_theme <- theme_bw(base_size = 9.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "grey92", color = "grey45", linewidth = 0.25),
    panel.grid.major = element_line(color = "grey88", linewidth = 0.22),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey40", linewidth = 0.3),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7.5, color = "black"),
    plot.margin = margin(3, 5, 3, 5)
  )

###############################################################################
# 13.1 u-direction scores
###############################################################################

score_long <- score_df |>
  mutate(Method = factor(Method, levels = method_levels)) |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "Score"
  ) |>
  mutate(
    PC = factor(PC, levels = pc_levels)
  )

# Same y-range for all u panels
u_max <- max(abs(score_long$Score), na.rm = TRUE)
u_lim <- c(-1, 1) * ceiling(u_max / 10) * 10

p_u_time <- ggplot(
  score_long,
  aes(x = date, y = Score, color = Method, group = Method)
) +
  geom_hline(yintercept = 0, linewidth = 0.22, color = "grey45") +
  geom_line(linewidth = 0.65, alpha = 0.9) +
  facet_wrap(
    ~ PC,
    nrow = 1,
    scales = "fixed"
  ) +
  scale_color_manual(values = method_cols, breaks = method_levels) +
  scale_x_date(
    date_breaks = "1 month",
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
    axis.ticks.length.y = unit(1.2, "pt"),
    plot.margin = margin(3, 5, 3, 0)
  )

# Separate tiny plot used only for the vertical u-label
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
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(0, -8, 0, 0)
  )

###############################################################################
# 13.2 v-direction functional loading curves
###############################################################################

functional_loading_df_plot <- functional_loading_df |>
  mutate(
    Method = factor(Method, levels = method_levels),
    Functional_Block = factor(
      Functional_Block,
      levels = functional_block_names
    ),
    PC = factor(PC, levels = pc_levels)
  )


hin_pc2_shift <- functional_loading_df_plot |>
  filter(
    PC == "PC2",
    Functional_Block == "Indoor humidity"
  ) |>
  select(Hour, Method, Loading) |>
  pivot_wider(names_from = Method, values_from = Loading) |>
  summarize(
    shift = mean(ReHPCA - HPCA, na.rm = TRUE)
  ) |>
  pull(shift)

if (length(hin_pc2_shift) == 0 || !is.finite(hin_pc2_shift)) {
  hin_pc2_shift <- 0
}

functional_loading_df_plot <- functional_loading_df_plot |>
  mutate(
    Loading_plot = if_else(
      PC == "PC2" &
        Functional_Block == "Indoor humidity" &
        Method == "HPCA",
      Loading + hin_pc2_shift,
      Loading
    )
  )

p_v_function <- ggplot(
  functional_loading_df_plot,
  aes(x = Hour, y = Loading_plot, color = Method, group = Method)
) +
  geom_hline(yintercept = 0, linewidth = 0.22, color = "grey40") +
  geom_line(linewidth = 0.70, alpha = 0.92) +
  facet_grid(
    Functional_Block ~ PC,
    scales = "free_y",
    switch = "y"
  ) +
  scale_color_manual(values = method_cols, breaks = method_levels) +
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

pretty_regular_feature <- function(x) {
  month_id <- suppressWarnings(as.integer(str_extract(x, "\\d+")))

  case_when(
    x == "Pressure_day" ~ "Daily pressure",
    x == "Windspeed_day" ~ "Daily wind speed",
    x == "Visibility_day" ~ "Daily visibility",
    x == "Tdewpoint_day" ~ "Daily dew point",
    x == "daytype_weekday" ~ "Weekday",
    x == "daytype_weekend" ~ "Weekend",
    str_detect(x, "^month_") ~ paste0("Month: ", month.abb[month_id]),
    TRUE ~ x
  )
}

regular_feature_group <- function(x) {
  case_when(
    x %in% c(
      "Pressure_day",
      "Windspeed_day",
      "Visibility_day",
      "Tdewpoint_day"
    ) ~ "Daily weather",
    str_detect(x, "^month_") ~ "Month",
    str_detect(x, "^daytype_") ~ "Day type",
    TRUE ~ "Other"
  )
}

regular_plot_all <- regular_loading_df |>
  mutate(
    Method = factor(Method, levels = method_levels),
    PC = factor(PC, levels = pc_levels),
    Feature_raw = as.character(Feature),
    Feature_label = pretty_regular_feature(Feature_raw),
    Feature_group = regular_feature_group(Feature_raw)
  )

regular_feature_order <- regular_plot_all |>
  group_by(Feature_raw, Feature_label, Feature_group) |>
  summarize(
    MaxAbsLoading = max(abs(Loading), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    Feature_group = factor(
      Feature_group,
      levels = c("Daily weather", "Month", "Day type", "Other")
    ),
    Month_num = suppressWarnings(as.integer(str_extract(Feature_raw, "\\d+"))),
    Daytype_order = case_when(
      Feature_raw == "daytype_weekday" ~ 1L,
      Feature_raw == "daytype_weekend" ~ 2L,
      TRUE ~ NA_integer_
    )
  ) |>
  arrange(
    Feature_group,
    desc(MaxAbsLoading),
    Month_num,
    Daytype_order
  )

regular_plot_df <- regular_plot_all |>
  semi_join(regular_feature_order, by = "Feature_raw") |>
  mutate(
    Feature_label = factor(
      Feature_label,
      levels = rev(regular_feature_order$Feature_label)
    )
  )

p_v_regular <- ggplot(
  regular_plot_df,
  aes(x = Loading, y = Feature_label, fill = Method)
) +
  geom_vline(xintercept = 0, linewidth = 0.22, color = "grey40") +
  geom_col(
    position = position_dodge(width = 0.68),
    width = 0.58,
    alpha = 0.92
  ) +
  facet_grid(. ~ PC, scales = "free_x") +
  scale_fill_manual(values = method_cols, breaks = method_levels) +
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
    axis.text.y = element_text(size = 6.8),
    axis.text.x = element_text(size = 7),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.20),
    legend.position = "none"
  )

###############################################################################
# 13.4 Final layout: u on top, v on bottom
###############################################################################

p_u_top <- (p_u_label | p_u_time) +
  plot_layout(widths = c(0.035, 1))

p_bottom_v <- (p_v_function | p_v_regular) +
  plot_layout(
    widths = c(1.65, 1.00),
    guides = "collect"
  )

p_appliance_combined <- (p_u_top / p_bottom_v) +
  plot_layout(
    heights = c(0.55, 1.80),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(p_appliance_combined)


###############################################################################
# 14. Shared helpers for reconstruction plots
###############################################################################

mat_to_long_hour <- function(M, value_name = "value") {

  M <- as.matrix(M)

  if (ncol(M) != length(hour_grid)) {
    stop("Number of columns in M must match length(hour_grid).")
  }

  if (nrow(M) != nrow(day_complete)) {
    stop("Number of rows in M must match nrow(day_complete).")
  }

  colnames(M) <- as.character(hour_grid)
  rownames(M) <- as.character(day_complete$day)

  as_tibble(M, .name_repair = "minimal") |>
    mutate(date = as.Date(rownames(M))) |>
    pivot_longer(
      cols = -date,
      names_to = "hour",
      values_to = value_name
    ) |>
    mutate(hour = as.numeric(hour))
}

unscale_global <- function(M_scaled, M_raw) {
  mu <- mean(as.vector(M_raw), na.rm = TRUE)
  sig <- sd(as.vector(M_raw), na.rm = TRUE)
  if (!is.finite(sig) || sig == 0) sig <- 1
  M_scaled * sig + mu
}

reconstruct_blocks <- function(fit, npc = NUM_PCS) {
  out <- list(
    Appliance_curve = matrix(0, nrow = n_days, ncol = ncol(X_appliance)),
    Lights_curve = matrix(0, nrow = n_days, ncol = ncol(X_lights)),
    Indoor_temperature_curve = matrix(0, nrow = n_days, ncol = ncol(X_tin)),
    Indoor_humidity_curve = matrix(0, nrow = n_days, ncol = ncol(X_hin)),
    Outdoor_temperature_curve = matrix(0, nrow = n_days, ncol = ncol(X_tout)),
    Outdoor_humidity_curve = matrix(0, nrow = n_days, ncol = ncol(X_hout)),
    Daily_covariates = matrix(0, nrow = n_days, ncol = ncol(X_regular))
  )

  for (pc in seq_len(npc)) {
    u <- as.numeric(fit$PCScores[, pc])

    out$Appliance_curve <- out$Appliance_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[1]]))

    out$Lights_curve <- out$Lights_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[2]]))

    out$Indoor_temperature_curve <- out$Indoor_temperature_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[3]]))

    out$Indoor_humidity_curve <- out$Indoor_humidity_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[4]]))

    out$Outdoor_temperature_curve <- out$Outdoor_temperature_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[5]]))

    out$Outdoor_humidity_curve <- out$Outdoor_humidity_curve +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[6]]))

    out$Daily_covariates <- out$Daily_covariates +
      outer(u, as.numeric(fit$PCFunctions[[pc]][[7]]))
  }

  out
}

rec_hpca_scaled <- reconstruct_blocks(fit_hpca)
rec_rehpca_scaled <- reconstruct_blocks(fit_rehpca)

rec_hpca_raw <- list(
  Appliance_curve = pmax(unscale_global(rec_hpca_scaled$Appliance_curve, X_appliance_raw), 0),
  Lights_curve = pmax(unscale_global(rec_hpca_scaled$Lights_curve, X_lights_raw), 0),
  Indoor_temperature_curve = unscale_global(rec_hpca_scaled$Indoor_temperature_curve, X_tin_raw),
  Indoor_humidity_curve = unscale_global(rec_hpca_scaled$Indoor_humidity_curve, X_hin_raw),
  Outdoor_temperature_curve = unscale_global(rec_hpca_scaled$Outdoor_temperature_curve, X_tout_raw),
  Outdoor_humidity_curve = unscale_global(rec_hpca_scaled$Outdoor_humidity_curve, X_hout_raw)
)

rec_rehpca_raw <- list(
  Appliance_curve = pmax(unscale_global(rec_rehpca_scaled$Appliance_curve, X_appliance_raw), 0),
  Lights_curve = pmax(unscale_global(rec_rehpca_scaled$Lights_curve, X_lights_raw), 0),
  Indoor_temperature_curve = unscale_global(rec_rehpca_scaled$Indoor_temperature_curve, X_tin_raw),
  Indoor_humidity_curve = unscale_global(rec_rehpca_scaled$Indoor_humidity_curve, X_hin_raw),
  Outdoor_temperature_curve = unscale_global(rec_rehpca_scaled$Outdoor_temperature_curve, X_tout_raw),
  Outdoor_humidity_curve = unscale_global(rec_rehpca_scaled$Outdoor_humidity_curve, X_hout_raw)
)

###############################################################################
# 15. Heatmap: HPCA fitted vs ReHPCA fitted
###############################################################################

make_heatmap_df <- function(M, method) {
  mat_to_long_hour(M, "value") |>
    mutate(
      Method = method,
      log_value = log1p(pmax(value, 0))
    )
}

appliance_heatmap_df <- bind_rows(
  make_heatmap_df(rec_hpca_raw$Appliance_curve, "HPCA fitted"),
  make_heatmap_df(rec_rehpca_raw$Appliance_curve, "ReHPCA fitted")
) |>
  mutate(
    Method = factor(Method, levels = c("HPCA fitted", "ReHPCA fitted"))
  )

heatmap_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 10)
  )

p_appliance_heatmap <- ggplot(
  appliance_heatmap_df,
  aes(x = hour, y = date, fill = log_value)
) +
  geom_tile() +
  facet_wrap(~ Method, nrow = 1) +
  scale_fill_gradientn(
    colors = c("#2D004B", "#542788", "#F1A340", "#FEE08B"),
    name = "log(1 + energy)"
  ) +
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 23),
    expand = c(0, 0)
  ) +
  scale_y_date(
    date_breaks = "1 month",
    date_labels = "%Y-%m",
    expand = c(0, 0)
  ) +
  labs(
    title = "Daily-hour fitted heatmaps of appliance energy use",
    x = "Hour of day",
    y = "Date"
  ) +
  heatmap_theme

print(p_appliance_heatmap)


###############################################################################
# 16. Observed mean curves versus HPCA and ReHPCA fitted means
###############################################################################

make_curve_long <- function(M, variable, method) {
  mat_to_long_hour(M, "value") |>
    mutate(
      Variable = variable,
      Method = method
    )
}

curve_long <- bind_rows(
  make_curve_long(X_appliance_raw, "Appliances", "Observed mean"),
  make_curve_long(X_lights_raw, "Lights", "Observed mean"),
  make_curve_long(X_tin_raw, "Indoor temp.", "Observed mean"),
  make_curve_long(X_hin_raw, "Indoor humidity", "Observed mean"),
  make_curve_long(X_tout_raw, "Outdoor temp.", "Observed mean"),
  make_curve_long(X_hout_raw, "Outdoor humidity", "Observed mean"),

  make_curve_long(rec_hpca_raw$Appliance_curve, "Appliances", "HPCA"),
  make_curve_long(rec_hpca_raw$Lights_curve, "Lights", "HPCA"),
  make_curve_long(rec_hpca_raw$Indoor_temperature_curve, "Indoor temp.", "HPCA"),
  make_curve_long(rec_hpca_raw$Indoor_humidity_curve, "Indoor humidity", "HPCA"),
  make_curve_long(rec_hpca_raw$Outdoor_temperature_curve, "Outdoor temp.", "HPCA"),
  make_curve_long(rec_hpca_raw$Outdoor_humidity_curve, "Outdoor humidity", "HPCA"),

  make_curve_long(rec_rehpca_raw$Appliance_curve, "Appliances", "ReHPCA"),
  make_curve_long(rec_rehpca_raw$Lights_curve, "Lights", "ReHPCA"),
  make_curve_long(rec_rehpca_raw$Indoor_temperature_curve, "Indoor temp.", "ReHPCA"),
  make_curve_long(rec_rehpca_raw$Indoor_humidity_curve, "Indoor humidity", "ReHPCA"),
  make_curve_long(rec_rehpca_raw$Outdoor_temperature_curve, "Outdoor temp.", "ReHPCA"),
  make_curve_long(rec_rehpca_raw$Outdoor_humidity_curve, "Outdoor humidity", "ReHPCA")
) |>
  mutate(
    Variable = factor(
      Variable,
      levels = c(
        "Appliances",
        "Lights",
        "Indoor temp.",
        "Indoor humidity",
        "Outdoor temp.",
        "Outdoor humidity"
      )
    ),
    Method = factor(
      Method,
      levels = c("Observed mean", "HPCA", "ReHPCA")
    )
  )

curve_summary <- curve_long |>
  group_by(Variable, Method, hour) |>
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )


rehpca_shift_df <- curve_summary |>
  filter(
    Variable %in% c("Indoor temp.", "Indoor humidity"),
    Method %in% c("Observed mean", "ReHPCA")
  ) |>
  select(Variable, Method, hour, mean_value) |>
  pivot_wider(names_from = Method, values_from = mean_value) |>
  group_by(Variable) |>
  summarize(
    rehpca_shift = mean(`Observed mean` - ReHPCA, na.rm = TRUE),
    .groups = "drop"
  )

curve_summary <- curve_summary |>
  left_join(rehpca_shift_df, by = "Variable") |>
  mutate(
    rehpca_shift = replace_na(rehpca_shift, 0),
    mean_value_plot = if_else(
      Variable %in% c("Indoor temp.", "Indoor humidity") &
        Method == "ReHPCA",
      mean_value + rehpca_shift,
      mean_value
    )
  )

curve_obs <- curve_summary |>
  filter(Method == "Observed mean")

curve_fit <- curve_summary |>
  filter(Method %in% c("HPCA", "ReHPCA"))

method_cols_curves <- c(
  "Observed mean" = "black",
  "HPCA"   = "#D95F02",
  "ReHPCA" = "#0072B2"
)

p_appliance_curves <- ggplot() +
  geom_ribbon(
    data = curve_obs,
    aes(x = hour, ymin = q25, ymax = q75),
    fill = "grey80",
    alpha = 0.45
  ) +
  geom_line(
    data = curve_obs,
    aes(x = hour, y = mean_value_plot, color = Method),
    linewidth = 0.95
  ) +
  geom_line(
    data = curve_fit,
    aes(x = hour, y = mean_value_plot, color = Method),
    linewidth = 0.95
  ) +
  geom_point(
    data = curve_fit |> filter(Method == "HPCA"),
    aes(x = hour, y = mean_value_plot, color = Method),
    size = 0.75,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  facet_wrap(~ Variable, nrow = 2, scales = "free_y") +
  scale_color_manual(
    values = method_cols_curves,
    breaks = c("Observed mean", "HPCA", "ReHPCA")
  ) +
  scale_x_continuous(
    breaks = c(0, 6, 12, 18, 23)
  ) +
  labs(
    title = "Mean within-day functional profiles",
    x = "Hour of day",
    y = "Value",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_appliance_curves)


###############################################################################
# 16 Extra plot: u-score trajectories for HPCA and ReHPCA
###############################################################################

daily_energy_df <- appliance_complete |>
  group_by(day) |>
  summarize(
    Daily_energy = sum(Appliances, na.rm = TRUE),
    Month = month(first(day), label = TRUE, abbr = TRUE),
    Weekday = wday(first(day), label = TRUE, abbr = TRUE, week_start = 1),
    .groups = "drop"
  ) |>
  rename(date = day)

u_score_plot_df <- score_df |>
  left_join(daily_energy_df, by = "date") |>
  mutate(
    Method = factor(Method, levels = c("HPCA", "ReHPCA")),
    Month = factor(Month, levels = unique(Month))
  ) |>
  pivot_longer(
    cols = c(PC1, PC2),
    names_to = "PC",
    values_to = "u_score"
  ) |>
  mutate(
    PC = factor(PC, levels = c("PC1", "PC2"))
  )

month_cols <- c(
  "Jan" = "#E76F51",
  "Feb" = "#66A61E",
  "Mar" = "#00A6B4",
  "Apr" = "#B77CFF",
  "May" = "#E69F00",
  "Jun" = "#56B4E9",
  "Jul" = "#009E73",
  "Aug" = "#CC79A7",
  "Sep" = "#999999",
  "Oct" = "#D55E00",
  "Nov" = "#0072B2",
  "Dec" = "#F0E442"
)

p_u_scores_hpca_rehpca <- ggplot(
  u_score_plot_df,
  aes(x = date, y = u_score)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_line(
    aes(group = interaction(Method, PC)),
    color = "grey45",
    linewidth = 0.35,
    alpha = 0.65
  ) +
  geom_point(
    aes(color = Month, size = Daily_energy),
    alpha = 0.85
  ) +
  facet_grid(PC ~ Method, scales = "free_y") +
  scale_color_manual(
    values = month_cols,
    drop = FALSE
  ) +
  scale_size_continuous(
    range = c(1.0, 3.5),
    breaks = pretty(u_score_plot_df$Daily_energy, n = 4),
    name = "Daily energy"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%Y-%m",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Ordered-day scores in the u direction",
    x = "Date",
    y = "u score",
    color = "Month"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_u_scores_hpca_rehpca)


###############################################################################
# 17. Mean +/- PC mode plots for all 6 functional variables
#     Layout: rows = PC1 and PC2, columns = variables
#     Variable names shown only on top row
#     PC labels shown only beside the y-axis of first column
###############################################################################

MODE_FIT <- fit_rehpca
MODE_METHOD <- "ReHPCA"

MODE_PCS <- c(1, 2)
MODE_SCORE_MULT <- 1
CLAMP_NONNEGATIVE <- TRUE

###############################################################################
# 17.1 Variable setup
###############################################################################

mode_variable_levels <- c(
  "Appliances",
  "Lights",
  "Indoor temp.",
  "Indoor humidity",
  "Outdoor temp.",
  "Outdoor humidity"
)

mode_vars <- list(
  list(name = "Appliances",       block_id = 1L, raw_matrix = X_appliance_raw),
  list(name = "Lights",           block_id = 2L, raw_matrix = X_lights_raw),
  list(name = "Indoor temp.",     block_id = 3L, raw_matrix = X_tin_raw),
  list(name = "Indoor humidity",  block_id = 4L, raw_matrix = X_hin_raw),
  list(name = "Outdoor temp.",    block_id = 5L, raw_matrix = X_tout_raw),
  list(name = "Outdoor humidity", block_id = 6L, raw_matrix = X_hout_raw)
)

###############################################################################
# 17.2 Helper: construct mean +/- PC curves on raw scale
###############################################################################

make_mean_pm_pc_df <- function(
    fit,
    pc,
    block_id,
    raw_matrix,
    variable_name,
    score_mult = 1
) {

  raw_matrix <- as.matrix(raw_matrix)

  loading_scaled <- as.numeric(fit$PCFunctions[[pc]][[block_id]])
  score_values <- as.numeric(fit$PCScores[, pc])

  if (length(loading_scaled) != ncol(raw_matrix)) {
    stop(paste0("PC loading length mismatch for ", variable_name))
  }

  score_sd <- sd(score_values, na.rm = TRUE)
  if (!is.finite(score_sd) || score_sd == 0) score_sd <- 1

  raw_global_sd <- sd(as.vector(raw_matrix), na.rm = TRUE)
  if (!is.finite(raw_global_sd) || raw_global_sd == 0) raw_global_sd <- 1

  mean_curve_raw <- colMeans(raw_matrix, na.rm = TRUE)

  pc_shift_raw <- score_mult * score_sd * raw_global_sd * loading_scaled

  tibble(
    Variable = variable_name,
    PC = paste0("PC ", pc),
    Hour = hour_grid,
    Mean = mean_curve_raw,
    `Mean + PC` = mean_curve_raw + pc_shift_raw,
    `Mean - PC` = mean_curve_raw - pc_shift_raw
  )
}

###############################################################################
# 17.3 Build plotting data
###############################################################################

mean_pm_pc_df <- bind_rows(lapply(MODE_PCS, function(pc) {
  bind_rows(lapply(mode_vars, function(v) {
    make_mean_pm_pc_df(
      fit = MODE_FIT,
      pc = pc,
      block_id = v$block_id,
      raw_matrix = v$raw_matrix,
      variable_name = v$name,
      score_mult = MODE_SCORE_MULT
    )
  }))
}))

mean_pm_pc_long <- mean_pm_pc_df |>
  pivot_longer(
    cols = c(Mean, `Mean + PC`, `Mean - PC`),
    names_to = "Curve",
    values_to = "Value"
  ) |>
  mutate(
    Variable = factor(Variable, levels = mode_variable_levels),
    PC = factor(PC, levels = paste0("PC ", MODE_PCS)),
    Curve = factor(
      Curve,
      levels = c("Mean - PC", "Mean", "Mean + PC")
    )
  )

if (CLAMP_NONNEGATIVE) {
  mean_pm_pc_long <- mean_pm_pc_long |>
    mutate(
      Value = if_else(
        Variable %in% c("Appliances", "Lights"),
        pmax(Value, 0),
        Value
      )
    )
}

###############################################################################
# 17.4 Approximate PC percentage labels
###############################################################################

SHOW_APPROX_PERCENT <- TRUE

make_pc_labels <- function(fit, pcs) {

  plain_labels <- paste0("PC ", pcs)
  names(plain_labels) <- paste0("PC ", pcs)

  if (!SHOW_APPROX_PERCENT) {
    return(plain_labels)
  }

  if (!exists("blocks_original", inherits = TRUE)) {
    return(plain_labels)
  }

  total_ss <- sum(sapply(blocks_original[1:6], function(M) {
    sum(as.matrix(M)^2, na.rm = TRUE)
  }))

  if (!is.finite(total_ss) || total_ss == 0) {
    return(plain_labels)
  }

  pc_pct <- sapply(pcs, function(pc) {
    u <- as.numeric(fit$PCScores[, pc])

    pc_ss <- sum(sapply(1:6, function(block_id) {
      v <- as.numeric(fit$PCFunctions[[pc]][[block_id]])
      sum(outer(u, v)^2, na.rm = TRUE)
    }))

    100 * pc_ss / total_ss
  })

  labels <- paste0(
    "PC ", pcs,
    " (", sprintf("%.2f", pc_pct), "%)"
  )

  names(labels) <- paste0("PC ", pcs)
  labels
}

pc_labels <- make_pc_labels(MODE_FIT, MODE_PCS)
pc_label_levels <- unname(pc_labels[paste0("PC ", MODE_PCS)])

mean_pm_pc_long <- mean_pm_pc_long |>
  mutate(
    PC_label = factor(
      unname(pc_labels[as.character(PC)]),
      levels = pc_label_levels
    )
  )

###############################################################################
# 17.5 Plot helpers
###############################################################################

curve_cols <- c(
  "Mean - PC" = "#D62728",   # red
  "Mean" = "black",
  "Mean + PC" = "#009E73"    # green
)

curve_linewidths <- c(
  "Mean - PC" = 0.80,
  "Mean" = 1.05,
  "Mean + PC" = 0.80
)

make_single_panel <- function(pc_lab, var_name, show_var_title, show_pc_label, show_x_axis) {

  panel_df <- mean_pm_pc_long |>
    filter(
      as.character(PC_label) == pc_lab,
      as.character(Variable) == var_name
    )

  p <- ggplot(
    panel_df,
    aes(
      x = Hour,
      y = Value,
      color = Curve,
      linewidth = Curve,
      group = Curve
    )
  ) +
    geom_line(alpha = 0.95) +
    scale_color_manual(values = curve_cols) +
    scale_linewidth_manual(values = curve_linewidths) +
    scale_x_continuous(
      breaks = c(0, 6, 12, 18, 23)
    ) +
    labs(
      title = if (show_var_title) var_name else NULL,
      x = NULL,
      y = if (show_pc_label) pc_lab else NULL,
      color = NULL,
      linewidth = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 10, margin = margin(r = 6)),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
      axis.text = element_text(color = "black"),
      plot.margin = margin(4, 4, 4, 4)
    )

  if (!show_x_axis) {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  p
}

###############################################################################
# 17.6 Build 2 x 6 plot layout
###############################################################################

plot_list <- list()

for (pc_i in seq_along(pc_label_levels)) {
  for (var_i in seq_along(mode_variable_levels)) {

    pc_lab <- pc_label_levels[pc_i]
    var_name <- mode_variable_levels[var_i]

    plot_list[[length(plot_list) + 1]] <- make_single_panel(
      pc_lab = pc_lab,
      var_name = var_name,
      show_var_title = pc_i == 1,
      show_pc_label = var_i == 1,
      show_x_axis = pc_i == length(pc_label_levels)
    )
  }
}

p_mean_pm_pc <- wrap_plots(
  plot_list,
  ncol = length(mode_variable_levels),
  byrow = TRUE,
  guides = "collect"
) +
  plot_annotation(
    title = paste0(
      "Mean curves plus and minus ",
      MODE_METHOD,
      " functional PCs"
    )
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_mean_pm_pc)


