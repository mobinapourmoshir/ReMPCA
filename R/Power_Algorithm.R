################ Sparse penalty for coefficients ################
# Lemma 2 (Sparse PCA via regularized low rank matrix approximation by Huang)
# y is either coefficients (u's) or PCs (v's)
sparse_pen_fun <- function(y,
                           tuning_parameter,
                           type,alpha = 3.7) {

  y_sorted <- sort(abs(y))
  lambda <- y_sorted[tuning_parameter]
  if (tuning_parameter == 0 ||
      tuning_parameter > length(y)) {
    return(y)
  }
  if (type == "soft") {
    return(sign(y) * pmax(abs(y) - lambda, 0))
  }
  else if (type == "hard") {
    return(ifelse(abs(y) > lambda, y, 0))
  }
  else if (type == "SCAD") {
    res <- ifelse(
      abs(y) <= 2 * lambda,
      sign(y) * pmax(abs(y) - lambda, 0),
      ifelse(
        abs(y) <= alpha * lambda,
        ((alpha - 1) * y - sign(y) * alpha * lambda) / (alpha - 2),
        y
      )
    )
    return(res)
  }
}

##################### Calculating the norm of a vector #####################
norm_vec <- function(x) sqrt(sum(x^2))

############################### Power Algorithm ###############################
power_algo <- function(data,
                       n_var,
                       ncol,
                       sparse_tuning_result_u,
                       sparse_tuning_result_v,
                       S_alpha_v,
                       S_alpha_u,
                       sparse_tuning_type,
                       conditional = FALSE,
                       alpha_Omega_v,
                       alpha_Omega_u,
                       thresh,
                       maxit) {

  rownames(data) <- NULL; colnames(data) <- NULL
  data <- as.matrix(data)

  start_col <- 1; splitted_data <- list()
  for (i in 1:n_var) {
    num_cols <- ncol[[i]]
    splitted_data[[i]] <- data[, start_col:(start_col + num_cols - 1)]
    start_col <- start_col + num_cols
  }

  v_old <- svd(data)$v[, 1]
  errors <- Inf
  iter <- 0  # initialize iteration counter

  while (errors > thresh && iter < maxit) {
    iter <- iter + 1

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
      denom_u <- as.numeric(t(v_old) %*%
                              (diag(nrow(alpha_Omega_v)) + alpha_Omega_v) %*%
                              v_old)
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
        v_part <- S_alpha_v[[i]] %*% as.vector(v_part_raw) / denom_v
      } else {
        v_part <- S_alpha_v[[i]] %*% as.vector(v_part_raw)
      }
      v_new <- c(v_new, v_part)
    }
    v_new <- v_new / norm_vec(v_new)
    errors <- sum((v_new - v_old)^2)
    v_old <- v_new
  }

  if (iter == maxit && errors > thresh) {
    warning("Algorithm did not converge within the maximum number of iterations.")
  }

  if (!conditional) u_new <- u_new / norm_vec(u_new)

  return(list(v_new = v_new,
              u_new = u_new,
              iterations = iter,
              error = errors))
}

