########## CV Scores calculators for sparsity on rows (PC scores, u) ##########
cv_sparse_row <- function(data,
                          n_var,
                          ncol,
                          S_alpha_u,
                          K_fold,
                          cv.pick,
                          thresh,
                          maxit,
                          conditional = FALSE,
                          sparse_tuning_result_u,  # vector of candidate sparsity levels for rows (u)
                          sparse_tuning_result_v,  # fixed sparsity params (vector)
                          sparse_tuning_type) {

  set.seed(123)
  shuffled_cols <- sample(ncol(data))  # Shuffle column indices
  group_size <- ceiling(length(shuffled_cols) / K_fold)
  m <- ncol(data)

  CV_errors <- numeric(length(sparse_tuning_result_u))
  SE_errors <- numeric(length(sparse_tuning_result_u))
  fold_errors_list <- vector("list", length(sparse_tuning_result_u))

  for (j in 1:length(sparse_tuning_result_u)) {
    gamma_u <- sparsity_row_list[j]
    fold_errors <- numeric(K_fold)

    for (k in 1:K_fold) {
      # Define indices for test columns
      cols_to_remove <- shuffled_cols[((k - 1) * group_size + 1):
                                        min(k * group_size, length(shuffled_cols))]

      # Create training and test sets
      data_train <- data[, -cols_to_remove]
      data_test  <- data[, cols_to_remove]

      # Update column structure
      updated_ncol <- update_ncol(ncol, cols_to_remove)

      # Create identity smoother matrices (no smoothing)
      S_alpha_v_train <- lapply(1:n_var, function(i) diag(updated_ncol[, i]))

      # Run power algorithm
      power_result <- power_algo(data = data_train,
                                 n_var = n_var,
                                 ncol = updated_ncol,
                                 thresh = thresh,
                                 maxit = maxit,
                                 conditional = FALSE,
                                 sparse_tuning_result_u = gamma_u,
                                 sparse_tuning_result_v = sparse_tuning_result_v,
                                 S_alpha_v = S_alpha_v_train,
                                 S_alpha_u = S_alpha_u,
                                 sparse_tuning_type = sparse_tuning_type)

      u_hat <- power_result[[2]]
      v_test <- t(as.matrix(data_test)) %*% u_hat

      # Compute reconstruction error for this fold
      fold_errors[k] <- sum((data_test - u_hat %*% t(v_test))^2) / ncol(data_test)
    }

    # Compute CV error and standard error for current gamma
    CV_errors[j] <- mean(fold_errors)
    SE_errors[j] <- sd(fold_errors) / sqrt(K_fold)
    fold_errors_list[[j]] <- fold_errors
  }

  # Select gamma_u based on cv.pick
  if (cv.pick == "1se") {
    j_star <- which.min(CV_errors)
    CV_1se_threshold <- CV_errors[j_star] + SE_errors[j_star]
    j_1se <- max(which(CV_errors <= CV_1se_threshold))  # most regularized
    gamma_u <- sparse_tuning_result_u[j_1se]
    } else if (cv.pick == "min") {
      j_min <- which.min(CV_errors)
      gamma_u <- sparse_tuning_result_u[j_min]
      } else {
        stop("cv.pick must be either '1se' or 'min'")
        }

  return(list(gamma_u,
              cv_results = data.frame(sparse_tuning_result_u,
                                      CV_errors,
                                      SE_errors)))
}

########### CV Scores calculators for sparsity on columns (PCs, v) ###########
cv_sparse_col <- function(data,
                          n_var,
                          ncol,
                          S_alpha_v,
                          S_alpha_u,
                          K_fold,
                          thresh,
                          maxit,
                          conditional = FALSE,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v,
                          sparse_tuning_type) {

  set.seed(123)
  shuffled_rows <- sample(nrow(data)) # Grouping rows of data matrix
  group_size <- ifelse(round(length(shuffled_rows) / K_fold, digits = 0) == 0,
                       1,
                       round(length(shuffled_rows) / K_fold, digits = 0))

  data_tilde <- data
  fold_errors <- numeric(K_fold)  # store errors

  for (k in 1:K_fold) {
    rows_to_remove <- shuffled_rows[((k - 1) * group_size + 1):min(k * group_size, length(shuffled_rows))]
    data_train <- data.frame(data_tilde[-rows_to_remove,])
    data_test <- data.frame(matrix(data_tilde[rows_to_remove,], nrow = length(rows_to_remove)))

    # Power Algorithm
    power_train <- power_algo(data = data_train,
                              n_var = n_var,
                              ncol = ncol,
                              conditional = FALSE,
                              thresh = thresh,
                              maxit = maxit,
                              sparse_tuning_result_u = sparse_tuning_result_u,
                              sparse_tuning_result_v = sparse_tuning_result_v,
                              S_alpha_v = S_alpha_v,
                              S_alpha_u = S_alpha_u[-rows_to_remove, -rows_to_remove],
                              sparse_tuning_type = sparse_tuning_type)

    v_train <- power_train[[1]]
    u_test <- as.matrix(data_test) %*% as.matrix(v_train)

    # Fold-wise error
    fold_errors[k] <- sum((data_test - u_test %*% t(v_train))^2) / nrow(data_test)
  }
  CV_mean <- mean(fold_errors)
  CV_se <- sd(fold_errors) / sqrt(K_fold)

  return(list(CV_error = CV_mean, SE = CV_se, fold_errors = fold_errors))
}

############# Computing number of columns when eliminating some ##############
update_ncol <- function(ncol_matrix, cols_to_remove) {
  col_counts <- as.numeric(ncol_matrix[1, ])
  cumulative_cols <- c(0, cumsum(col_counts))
  updated_ncol <- col_counts
  # Loop through each dataset and count removed columns
  for (i in seq_along(updated_ncol)) {
    start_idx <- cumulative_cols[i] + 1
    end_idx <- cumulative_cols[i + 1]

    # Count how many columns to remove in this dataset
    removed_count <- sum(cols_to_remove >= start_idx & cols_to_remove <= end_idx)
    updated_ncol[i] <- updated_ncol[i] - removed_count
  }
  return(as.data.frame(t(updated_ncol)))
}
