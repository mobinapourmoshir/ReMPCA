############################### Sparse penalty for coefficients ###############################
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

############################### Calculating the norm of a vector ###############################
norm_vec <- function(x) sqrt(sum(x^2))

############################### Power Algorithm ###############################
power_algo = function(data,
                      n_var,
                      ncol,
                      sparse_tuning_result_u, # A fixed number
                      sparse_tuning_result_v, # A vector (length = p)
                      S_alpha_v, # A list (lenght = p)
                      S_alpha_u, # A matrix
                      sparse_tuning_type){

  rownames(data) <- NULL; colnames(data) <- NULL
  data <- as.matrix(data)

  # Splitting data
  start_col <- 1; splitted_data <- list()
  for (i in 1:n_var) {
    num_cols <- ncol[[i]]  # Get number of columns for the ith matrix
    splitted_data[[i]] <- data[, start_col:(start_col + num_cols - 1)]  # Extract matrix
    start_col <- start_col + num_cols
  }

  v_old <- svd(as.matrix(data))$v[, 1]
  errors <- Inf; thresh <- 1e-10

  # Power Algorithm
  while (errors > thresh) {
    u_new <- S_alpha_u %*% sparse_pen_fun(y = data%*%as.matrix(v_old),
                                          tuning_parameter = sparse_tuning_result_u,
                                          type = sparse_tuning_type) # u = h_{gamma} Xv
    if(all(u_new == 0)){ # To avoid 0 in the denominator
      u_new <- S_alpha_u %*% data%*%v_old
      }
    v_new <- c()
    for(i in 1:n_var){
      sparse_param <- as.integer(sparse_tuning_result_v[i])
      v <- S_alpha_v[[i]] %*% sparse_pen_fun(y = t(splitted_data[[i]]) %*% as.matrix(u_new),
                                             tuning_parameter = sparse_param,
                                             type = sparse_tuning_type) # v = S_{alpha} h_{gamma} t(X)u

      if(all(v == 0)){ # To avoid 0 in the denominator
        v <- S_alpha_v[[i]] %*% t(splitted_data[[i]]) %*% as.matrix(u_new)}
      v_new <- c(v_new, v)
      }

    # Adjust the sign of v based on the direction of maximum variance in the original data
    #max_var_index = which.max(apply(data, 2, var))
    #v_new = v_new * sign(v_new[max_var_index])
    v_new = v_new / norm_vec(v_new) # v/||v||

    # Convergence condition
    errors = sum((v_new - v_old)^2)
    v_old <- v_new
  }
  u_new = u_new/norm_vec(u_new) # u/||u||
  return(list(v_new,u_new))
}

