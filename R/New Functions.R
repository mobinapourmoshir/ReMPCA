cv_sparse_col_SE <- function(data,
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





cv_sparse_row_SE <- function(data,
                             n_var,
                             ncol,
                             S_alpha_u,
                             K_fold,
                             sparse_tuning_result_u,  # vector of candidate sparsity levels for rows (u)
                             sparse_tuning_result_v,  # fixed right sparsity (vector)
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
      cols_to_remove <- shuffled_cols[((k - 1) * group_size + 1):min(k * group_size, length(shuffled_cols))]

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

  # One-standard-error rule selection
  j_star <- which.min(CV_errors)
  CV_1se_threshold <- CV_errors[j_star] + SE_errors[j_star]
  j_1se <- max(which(CV_errors <= CV_1se_threshold))  # most regularized (largest gamma)
  gamma_u <- sparse_tuning_result_u[j_1se]

  return(list(gamma_u, data.frame(sparse_tuning_result_u, CV_errors, SE_errors)))
}




library(Matrix)
power_algo2 <- function(data,
                        n_var,
                        ncol,
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        S_alpha_v,
                        S_alpha_u,
                        sparse_tuning_type,
                        conditional = FALSE,
                        alpha_Omega_v,
                        alpha_Omega_u) {
  rownames(data) <- NULL; colnames(data) <- NULL
  data <- as.matrix(data)

  start_col <- 1; splitted_data <- list()
  for (i in 1:n_var) {
    num_cols <- ncol[[i]]
    splitted_data[[i]] <- data[, start_col:(start_col + num_cols - 1)]
    start_col <- start_col + num_cols
  }

  v_old <- svd(data)$v[, 1]
  errors <- Inf; thresh <- 1e-10

  while (errors > thresh) {
    # Compute u
    Xv <- data %*% v_old
    u_raw <- sparse_pen_fun(y = Xv,
                            tuning_parameter = sparse_tuning_result_u,
                            type = sparse_tuning_type)
    if (all(u_raw == 0)) u_raw <- Xv

    if (conditional) {
      v_old <- as.vector(v_old)
      alpha_Omega_v <- as.matrix(alpha_Omega_v)
      storage.mode(alpha_Omega_v) <- "numeric"
      denom_u <- as.numeric(t(v_old) %*% (diag(nrow(alpha_Omega_v)) + alpha_Omega_v) %*% v_old)
      u_new <- S_alpha_u %*% u_raw / denom_u
    } else {
      u_new <- S_alpha_u %*% u_raw
    }

    # Compute v
    v_new <- c()
    for (i in 1:n_var) {
      sparse_param <- as.integer(sparse_tuning_result_v[i])
      Xu_i <- t(splitted_data[[i]]) %*% u_new
      v_part_raw <- sparse_pen_fun(y = Xu_i,
                                   tuning_parameter = sparse_param,
                                   type = sparse_tuning_type)
      if (all(v_part_raw == 0)) v_part_raw <- Xu_i

      if (conditional) {
        norm_u2 <- as.numeric(norm_vec(u_new)^2)
        R_u <- as.numeric(t(u_new) %*% alpha_Omega_u %*% u_new / norm_u2)
        denom_v <- 1 + R_u
        v_part <- S_alpha_v[[i]] %*% v_part_raw / denom_v
      } else {
        v_part <- S_alpha_v[[i]] %*% v_part_raw
      }

      v_new <- c(v_new, v_part)
    }

    v_new <- v_new / norm_vec(v_new)
    errors <- sum((v_new - v_old)^2)
    v_old <- v_new
  }

  if (!conditional) u_new <- u_new / norm_vec(u_new)
  return(list(v_new, u_new))
}



opt_alpha_v <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v,   # List of length = n_iter, each element is list of length = n_var
                        S_alphas_u,   # Matrix (n × n)
                        alphas_v,     # Matrix of alpha_v values (n_iter × n_var)
                        alpha_u,      # scalar
                        Omega_u,      # A matrix associated with alpha_u
                        Omegas_v,     # List of length = n_iter, each element is list of length = n_var
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- nrow(alphas_v)
  n <- nrow(X)
  m <- ncol(X)
  GCV <- numeric(n_iter)

  if (all(alphas_v == 0)) {
    return(list(GCV_v = Inf,
                opt.alpha_v = rep(0, n_var),
                opt_s.alpha_v = replicate(n_var, diag(m / n_var), simplify = FALSE),
                GCVdf_v = data.frame(alphas_v, rep(Inf, n_iter))))
  } else {
    alpha_Omega_u <- alpha_u * Omega_u
    for (i in 1:n_iter) {
      S_alpha_v <- S_alphas_v[[i]]
      Omega_v <- Omegas_v[[i]]
      alpha_Omega_v <- as.matrix(bdiag(lapply(1:n_var, function(j) alphas_v[i, j] * Omega_v[[j]])))

      # Power Algorithm
      power_result <- power_algo2(data = X,
                                  n_var = n_var,
                                  ncol = ncol,
                                  sparse_tuning_result_u = sparse_tuning_result_u,
                                  sparse_tuning_result_v = sparse_tuning_result_v,
                                  S_alpha_v = S_alpha_v,
                                  S_alpha_u = S_alphas_u,
                                  sparse_tuning_type = sparse_tuning_type,
                                  conditional = TRUE,
                                  alpha_Omega_v = alpha_Omega_v,
                                  alpha_Omega_u = alpha_Omega_u)

      v <- power_result[[1]]
      u <- power_result[[2]]

      norm_u2 <- as.numeric(norm_vec(u)^2)
      Xu_proj <- as.vector(t(X) %*% u / norm_u2)
      R_u <- as.numeric(t(u) %*% Omega_u %*% u / norm_u2)

      GCV_alpha <- 0
      idx <- 0
      for (k in 1:n_var) {
        m_k <- as.integer(ncol[k])
        v_k <- v[(idx + 1):(idx + m_k)]
        xuk <- Xu_proj[(idx + 1):(idx + m_k)]
        S_k <- S_alpha_v[[k]]

        trace_k <- sum(diag(S_k))
        denom_k <- (1 - (1 / m_k) * (trace_k / (1 + alpha_u * R_u)))^2
        if(denom_k == 0){denom_k<-1}
        numer_k <- norm_vec(xuk - v_k)^2

        GCV_alpha <- GCV_alpha + numer_k / denom_k
        idx <- idx + m_k
      }

      GCV[i] <- GCV_alpha / m
    }

    opt.alpha <- alphas_v[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]

    return(list(GCV_v = GCV,
                opt.alpha_v = opt.alpha,
                opt_s.alpha_v = opt_s.alpha,
                GCVdf_v = data.frame(alphas_v, GCV)))
  }
}



opt_alpha_u <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas_u,   # A vector of alpha_u
                        alpha_v,    # scalar
                        Omega_v,    # A matrix associated with alpha_v
                        Omegas_u,   # List of length, each element is matrix
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- length(alphas_u)
  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas_u == 0)) {
    return(list(GCV_u = Inf, opt.alpha_u = 0, opt_s.alpha_u = diag(n), GCVdf_u = data.frame(alphas_u, rep(Inf, n_iter))))
  } else {
    alpha_Omega_v <- alpha_v * Omega_v
    for (i in 1:n_iter) {
      GCV_alpha <- 0
      S <- S_alphas_u[[i]]
      Omega_u <- Omegas_u[[i]]
      alpha_Omega_u <- as.matrix(alphas_u[i] * Omega_u)

      power_result <- power_algo2(data = X,
                                  n_var,
                                  ncol,
                                  sparse_tuning_result_u, # A fixed number
                                  sparse_tuning_result_v, # A vector (length = p)
                                  S_alpha_v = list(S_alphas_v), # A list (lenght = p)
                                  S_alpha_u = S, # A matrix
                                  conditional = TRUE,
                                  alpha_Omega_v = alpha_Omega_v,
                                  alpha_Omega_u = alpha_Omega_u,
                                  sparse_tuning_type = sparse_tuning_type)

      u <- power_result[[2]]
      v <- t(X) %*% u
      v <- as.vector(v)
      Omega_v <- as.matrix(Omega_v)
      storage.mode(Omega_v) <- "numeric"

      norm_v2 <- as.numeric(norm_vec(v)^2)
      R_v <- as.numeric(t(v) %*% Omega_v %*% v / norm_v2)

      GCV_alpha <- ( ((1/n) * ((X %*% v) / norm_vec(v)) - u )^2 ) / ( 1 -  (1/n)* (sum(diag(S)))/(1 + alpha_v*R_v))^2
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas_u[which.min(GCV)]
    opt_s.alpha <- S_alphas_u[[which.min(GCV)]]
    Omega_u <- Omega_u[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV_u = GCV, opt.alpha_u = opt.alpha, opt_s.alpha_u = opt_s.alpha, GCVdf = data.frame(alphas_u, GCV), Omega_u = Omega_u)

    return(result)
  }
}
