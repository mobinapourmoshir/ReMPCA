################# CV Scores calculators for sparsity on rows (PC scores, u) #################
cv_sparse_row <- function(data, # Hybrid data
                          S_alpha_v,
                          S_alpha_u,
                          K_fold,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v = 0,
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

    # Power Algorithm
    u_test <- power_algo(data = data_train,
                         sparse_tuning_result_u = sparse_tuning_result_u,
                         sparse_tuning_result_v = 0,
                         sparse_tuning_type = sparse_tuning_type,
                         S_alpha_v,
                         S_alpha_u,
                         type = type)

    v_test <- t(as.matrix(data_test))%*%as.matrix(u_test)

    # Calculating CV Score
    error_score_sparse = error_score_sparse + sum((
      data_test - u_test %*% t(v_test)) ^ 2)
  }
  return(error_score_sparse / nrow(data))  # Assuming nrow(data) is N in the formula
}



################# CV Scores calculators for sparsity on columns (PCs, v) #################
cv_sparse_col <- function(data, # Hybrid data
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

    # Power Algorithm
    u_test <- power_algo(data = data_train,
                         sparse_tuning_result_u = sparse_tuning_result_u,
                         sparse_tuning_result_v = 0,
                         sparse_tuning_type = sparse_tuning_type,
                         S_alpha_v,
                         S_alpha_u,
                         type = type)

    v_test <- t(as.matrix(data_test))%*%as.matrix(u_test)

    # Calculating CV Score
    error_score_sparse = error_score_sparse + sum((
      data_test - u_test %*% t(v_test)) ^ 2)
  }
  return(error_score_sparse / nrow(data))  # Assuming nrow(data) is N in the formula
}



















