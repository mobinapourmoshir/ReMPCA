################# CV Scores calculators for sparsity on rows (PC scores, u) #################
cv_sparse_row <- function(data, # Hybrid data
                          n_var,
                          ncol,
                          S_alpha_v,
                          S_alpha_u,
                          K_fold,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v,
                          sparse_tuning_type,
                          type) {
  set.seed(123)
  shuffled_cols <- sample(ncol(data)) # Grouping the columns of data matrix
  group_size <- ifelse(round(length(shuffled_cols)/ K_fold,
                             digits = 0) == 0, 1,
                       round(length(shuffled_cols)/ K_fold, digits = 0))

  data_tilde <- data
  error_score_sparse <- 0

  for (k in 1:K_fold) {
    cols_to_remove <- shuffled_cols[((k - 1) * group_size + 1):((k) * group_size)]
    data_train <- data.frame(data_tilde[,-cols_to_remove])  # X^{-k}
    data_test <- data.frame(matrix(data_tilde[, cols_to_remove], ncol = length(cols_to_remove))) # X^k

    updated_ncol <- update_ncol(ncol, cols_to_remove)
    S_alpha_v <- list()
    for (i in 1:n_var) {
      S_alpha_v[[i]] <- diag(updated_ncol[,i])
    }
    # Power Algorithm
    power_train <- power_algo(data = data_train,
                         n_var = n_var,
                         ncol = updated_ncol,
                         sparse_tuning_result_u = sparse_tuning_result_u, # A fixed number
                         sparse_tuning_result_v = sparse_tuning_result_v, # A vector (length = p)
                         S_alpha_v = S_alpha_v, # A list (lenght = p)
                         S_alpha_u = S_alpha_u, # A matrix
                         sparse_tuning_type = sparse_tuning_type)

    u_train <- power_train[[2]]
    v_test <- t(as.matrix(data_test)) %*% as.matrix(u_train)

    # Calculating CV Score
    error_score_sparse = error_score_sparse + sum(((
      data_test - u_train %*% t(v_test))^ 2)/ ncol(data_test))

  }
  return(error_score_sparse)  # Assuming ncol(data) is N in the formula
}



################# CV Scores calculators for sparsity on columns (PCs, v) #################
cv_sparse_col <- function(data,
                          n_var,
                          ncol,
                          S_alpha_v,
                          S_alpha_u,
                          K_fold,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v,
                          sparse_tuning_type) {

  set.seed(123)
  shuffled_rows <- sample(nrow(data)) # Grouping rows of data matrix
  group_size <- ifelse(round(length(shuffled_rows)/ K_fold,
                             digits = 0) == 0, 1,
                       round(length(shuffled_rows)/ K_fold, digits = 0))

  data_tilde <- data
  error_score_sparse <- 0

  for (k in 1:K_fold) {
    rows_to_remove <- shuffled_rows[((k - 1) * group_size + 1):((k) * group_size)]
    data_train <- data.frame(data_tilde[-rows_to_remove,])  # X^{-k}
    data_test <- data.frame(matrix(data_tilde[rows_to_remove,], nrow =length(rows_to_remove))) # X^k

    # Power Algorithm
    power_train <- power_algo(data = data_train,
                              n_var = n_var,
                              ncol = ncol,
                              sparse_tuning_result_u = sparse_tuning_result_u,
                              sparse_tuning_result_v = sparse_tuning_result_v, # A vector
                              S_alpha_v = S_alpha_v, # A list
                              S_alpha_u = S_alpha_u[-rows_to_remove,-rows_to_remove], # A matrix
                              sparse_tuning_type = sparse_tuning_type)

    v_train <- power_train[[1]]
    u_test <- as.matrix(data_test) %*% as.matrix(v_train)

    # Calculating CV Score
    error_score_sparse = error_score_sparse + sum(((
      data_test - u_test %*% v_train) ^ 2)/ nrow(data_test))
  }
  return(error_score_sparse)
}


################# Computing number of columns when eliminating some #################
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

