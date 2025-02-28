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
                      sparse_tuning_result_u,
                      sparse_tuning_result_v,
                      S_alpha_v,
                      S_alpha_u,
                      sparse_tuning_type,
                      type = "real"){

  rownames(data) <- NULL; colnames(data) <- NULL
  data <- as.matrix(data)
  v_old <- svd(as.matrix(data))$v[, 1]
  errors <- 10^60; thresh <- 1e-10

  # Power Algorithm
  while (errors > thresh) {
    if (type == "CV") {
      u_new <- sparse_pen_fun(y = data%*%as.matrix(v_old),
                              tuning_parameter = sparse_tuning_result_u,
                              type = sparse_tuning_type) # u = h_{gamma} Xv

      if(all(u_new == 0)){u_new <- data%*%v_old}
      v_new = sparse_pen_fun(y = t(data)%*%u_new,
                             tuning_parameter = sparse_tuning_result_v,
                             type = sparse_tuning_type)
      if(all(v_new == 0)){v_new <- t(data)%*%u_new} # To avoid 0 in the denominator

    }else{
      u_new <- S_alpha_u %*% sparse_pen_fun(y = data%*%v_old,
                              tuning_parameter = sparse_tuning_result_u,
                              type = sparse_tuning_type) # u = h_{gamma} Xv
      if(all(u_new == 0)){u_new <- data%*%v_old}

      v_new = S_alpha_v %*% sparse_pen_fun(y = t(data)%*%u_new,
                                         tuning_parameter = sparse_tuning_result_v,
                                         type = sparse_tuning_type) # v = S_{alpha} h_{gamma} t(X)u
      if(all(v_new == 0)){v_new <- t(data)%*%u_new}
      # Adjust the sign of v based on the direction of maximum variance in the original data
      #max_var_index = which.max(apply(data, 2, var))
      #v_new = v_new * sign(v_new[max_var_index])
    }
    v_new = v_new / norm_vec(v_new) # v/||v||



    # Convergence condition
    errors = sum((v_new - v_old)^2)
    v_old <- v_new
  }

  u_new = u_new/norm_vec(u_new) # u/||u||
  if (type == "CV") {
    return(u_new)
  }
  #if (type == "CV-two-way") {
  #  return(v_new)
  #}
  else{
    # v_new = v_new %*% solve(sqrt(t(v_new) %*% solve(S_alpha) %*% v_new))
    return(list(v_new,u_new))
  }
}

