############################### Power Algorithm ###############################
power_algo = function(data,sparse_tuning_result,sparse_tuning_type,S_alpha = NULL,type = "real"){

  v_old = svd(data)$v[,1]
  errors = 10^60; thresh <- 1e-10

  # Power Algorithm
  while (errors > thresh) {
    u_new = csparse_pen_fun(y = as.vector(data%*%v_old),tuning_parameter = sparse_tuning_result,sparse_tuning_type) # u = h_{gamma} Xv
    if (type == "CV") {
      v_new = t(data)%*%u_new
    } else{
      v_new = S_alpha %*% t(data) %*% u_new # v = S_{alpha}t(X)u
    }
    v_new = v_new / norm_vec(v_new) # v/||v||

    # Adjust the sign of v based on the direction of maximum variance in the original data
    max_var_index = which.max(apply(data, 2, var))
    v_new = v_new * sign(v_new[max_var_index])

    # Convergence condition
    errors = sum((v_new - v_old)^2)
    v_old = v_new
  }

  u_new = u_new/norm_vec(u_new) # u/||u||
  if (type == "CV") {
    return(u_new)
  }
  else{
    # v_new = v_new %*% solve(sqrt(t(v_new) %*% solve(S_alpha) %*% v_new))
    return(list(v_new,u_new))
  }
}
