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



############################### Power Algorithm with tuning Parameters ###############################

Tuning_Power <- function(X_temp =  X_temp, Y_temp = Y_temp , nvar = n_var, ncol = n_cols, n = n, smooth_tuning = smooth_tuning,
                         sparse_tuning_u_fd  = sparse_tuning_u_fd, sparse_tuning_v_fd = sparse_tuning_v_fd,
                         sparse_tuning_u_nfd  = sparse_tuning_u_nfd, sparse_tuning_v_nfd = sparse_tuning_nv_fd,
                         sparse_tuning_type = sparse_tuning_type, K_fold = K_fold,
                         S_alpha_List_v = S_alpha_list_v,S_alpha_list_u = S_alpha_list_u ,
                         two_way_smoothness = two_way_smoothness , two_way_sparsity = two_way_sparsity, j = j){


  # Functional data only!
  if(!(is.null(X_temp))){
    smooth_tuning_result_fd  <- sparse_tuning_result_fd <- list()
    gcv_fd <- opt_S_fd  <- funcs_fd <- GCVdf_fd <- list()


    lsv_fd <- lsu_fd <- c() # List for storing v's  and u's
    variance_fd <- vector() # % of variability explained by PC


    # Tuning Parameters ########### Functional data ###########
    opt_parameters_result <- opt_alpha_result <- list()
    opt_parameters_result <- parameter_selection_conditional(data = X_temp, nvar = n_var, ncol = n_cols, n = n, smooth_tuning = smooth_tuning,
                                                             sparse_tuning_u  = sparse_tuning_u_fd, sparse_tuning_v = sparse_tuning_v_fd,
                                                             sparse_tuning_type = sparse_tuning_type, K_fold = K_fold,
                                                             S_alpha_List_v = S_alpha_list_v,S_alpha_list_u = S_alpha_list_u ,
                                                             two_way_smoothness = two_way_smoothness , two_way_sparsity = two_way_sparsity)



    ################# Take care of sparse_tuning = 0 #################
    sparse_tuning_result_fd[[j]] <- opt_parameters_result[[1]] # Optimal level of sparsity (CV)
    opt_alpha_result <- opt_parameters_result[[2]] # Optimal Smoothness (GCV)

    opt_S_fd[[j]] <- opt_alpha_result$opt_s.alpha
    smooth_tuning_result_fd[[j]] <- opt_alpha_result$opt.alpha
    gcv_fd[[j]] <- opt_alpha_result$GCV
    GCVdf_fd[[j]] <- opt_alpha_result$GCVdf



    # Extracting v and u having the optimal parameters
    test_result <- power_algo(data = X_temp, sparse_tuning_result_fd = sparse_tuning_result_fd[[j]] ,
                              sparse_tuning_type = sparse_tuning_type, S_alpha = opt_alpha_result$opt_s.alpha, type = "real")

    v_fd <- test_result[[1]]
    u_fd <- test_result[[2]]
    lsv_fd <- cbind(lsv_fd, v)
    lsu_fd <- cbind(lsu_fd, u)
    funcs[[j]] <- u_fd%*%t(v_fd)



    result_fd <- list(Estimated = funcs, PC_Scores = lsu_fd, lsv_fd = lsv_fd, opt_alpha_for_PC = smooth_tuning_result_fd,
                   opt_gamma_for_PC = sparse_tuning_result_fd, GCV = gcv_fd, GCV_df = GCVdf_fd)


    # Non-Functional data only!
  }else if(!(is.null(Y_temp))){

    smooth_tuning_result_nfd  <- sparse_tuning_result_nfd <- list()
    gcv_nfd <- opt_S_nfd  <- funcs_nfd <- GCVdf_nfd <- list()


    lsv_nfd <- lsu_nfd <- c() # List for storing v's  and u's
    variance_nfd <- vector() # % of variability explained by PC


    # Tuning Parameters ########### Functional data ###########
    opt_parameters_result <- opt_alpha_result <- list()
    opt_parameters_result <- parameter_selection_conditional(data = Y_temp, nvar = ncol(Y_temp) , ncol = ncol(Y_temp), n = n, smooth_tuning = smooth_tuning,
                                                             sparse_tuning_u  = sparse_tuning_u_nfd, sparse_tuning_v = sparse_tuning_v_nfd,
                                                             sparse_tuning_type = sparse_tuning_type, K_fold = K_fold,
                                                             S_alpha_List_v = S_alpha_list_v,S_alpha_list_u = S_alpha_list_u ,
                                                             two_way_smoothness = two_way_smoothness , two_way_sparsity = two_way_sparsity)



    ################# Take care of sparse_tuning = 0 #################
    sparse_tuning_result_nfd[[j]] <- opt_parameters_result[[1]] # Optimal level of sparsity (CV)
    opt_alpha_result <- opt_parameters_result[[2]] # Optimal Smoothness (GCV)

    opt_S_nfd[[j]] <- opt_alpha_result$opt_s.alpha
    smooth_tuning_result_nfd[[j]] <- opt_alpha_result$opt.alpha
    gcv_nfd[[j]] <- opt_alpha_result$GCV
    GCVdf_nfd[[j]] <- opt_alpha_result$GCVdf



    # Extracting v and u having the optimal parameters
    test_result_nfd <- power_algo(data = Y_temp, sparse_tuning_result_nfd = sparse_tuning_result_nfd[[j]] ,
                              sparse_tuning_type = sparse_tuning_type, S_alpha = opt_alpha_result$opt_s.alpha, type = "real")

    v_nfd <- test_result_nfd[[1]]
    u_nfd <- test_result_nfd[[2]]
    lsv_nfd <- cbind(lsv_nfd, v)
    lsu_nfd <- cbind(lsu_nfd, u)
    funcs[[j]] <- u_nfd%*%t(v_nfd)



    result_nfd <- list(Estimated = funcs, PC_Scores = lsu_nfd, lsv_nfd = lsv_nfd , opt_alpha_for_PC = smooth_tuning_result_nfd,
                   opt_gamma_for_PC = sparse_tuning_result_nfd, GCV = gcv_nfd, GCV_df = GCVdf_nfd)

  }


  return(list(result_fd, result_nfd))


}



