############################ Conditional Tuning Parameters  - CV and GCV ############################
parameter_selection_conditional <- function(X_temp,
                                            n_var,
                                            ncol,
                                            n,
                                            smooth_tuning_v, # expand.grid matrix
                                            smooth_tuning_u, # Vector
                                            sparse_tuning_u, # vector
                                            sparse_tuning_v, # list
                                            sparse_tuning_type,
                                            K_fold,
                                            S_alpha_list_v ,
                                            S_alpha_list_u,
                                            Omegas_u,
                                            tuning_order){



  # tuning_order = 'Sparsity'
  if(tuning_order == "Sparsity"){
    # No smoothness list:
    S_alpha_v <- list()
    for (i in 1:n_var) {
      S_alpha_v[[i]] <- diag(ncol[,i])
    }

    CV_score_sparse_u <- CV_score_sparse_v <- GCV_score_smooth_u <- GCV_score_smooth_v <- Inf
    CV_scores_result_u <- CV_scores_result_v <- c()
    result = c()

    # Step 1: Sparsity on u  with no smoothness
    for (sparse_tuning_single in sparse_tuning_u) {
      sparse_score = cv_sparse_row(data=X_temp,
                                   ncol = ncol,
                                   n_var = n_var,
                                   S_alpha_v = S_alpha_v,
                                   S_alpha_u = diag(nrow(X_temp)),
                                   K_fold = K_fold,
                                   sparse_tuning_result_u = sparse_tuning_single,
                                   sparse_tuning_result_v = rep(0, n_var),
                                   sparse_tuning_type = sparse_tuning_type)
      CV_scores_result_u <- c(CV_scores_result_u, sparse_score)
      if (sparse_score <= CV_score_sparse_u) {
        CV_score_sparse_u = sparse_score
        gamma_u = sparse_tuning_single
      }
    }

    # Step 2: Conditional sparsity on v having sparsity on u with no smoothness
    gamma_v <- rep(0, n_var)
    for (i in 1:n_var) {
      gamma_Xi <- sparsity_col_list[[i]]
      CV_scores_result_v <- c()
      for (sparse_tuning_single in gamma_Xi) {
        gamma_v[i] <- sparse_tuning_single
        sparse_score = cv_sparse_col(data = X_temp,
                                     n_var = n_var,
                                     ncol = ncol,
                                     S_alpha_v = S_alpha_v,
                                     S_alpha_u = diag(nrow(X_temp)),
                                     K_fold = K_fold,
                                     sparse_tuning_result_u = gamma_u,
                                     sparse_tuning_result_v = gamma_v,
                                     sparse_tuning_type = sparse_tuning_type)

        CV_scores_result_v <- c(CV_scores_result_v, sparse_score)
        if (sparse_score <= CV_score_sparse_v) {
          CV_score_sparse_v = sparse_score
          sparse_tuning_selection_v = sparse_tuning_single
        }
        gamma_v[i] <- gamma_Xi[which.min(CV_scores_result_v)]
      }
    }
  }


  # Step 3: Conditional smoothness on u having sparsity on u and v and no smoothness on v
  GCV_result_u <- opt_alpha_u(X = X_temp,
                            n_var = n_var,
                            ncol = ncol,
                            S_alphas_v = S_alpha_v,
                            S_alphas_u = S_alpha_list_u,
                            alphas = smooth_tuning_u,
                            sparse_tuning_result_u = gamma_u,
                            sparse_tuning_result_v = gamma_v,
                            sparse_tuning_type = sparse_tuning_type)











  ###### Smoothing tuning parameter using GCV ######
  GCV_score_smooth = opt_alpha(X = X_temp,
                                 n_var = n_var,
                                 ncol = ncol,
                                 S_alphas_v = S_alpha_list_v,
                                 S_alphas_u = S_alpha_list_u,
                                 alphas = smooth_tuning,
                                 sparse_tuning_result_u = sparse_tuning_selection_u,
                                 sparse_tuning_result_v = sparse_tuning_selection_v,
                                 sparse_tuning_type = sparse_tuning_type,
                                 two_way_smoothness)



  close(pb) # Close the connection
  result = list(sparse_tuning_selection_u = sparse_tuning_selection_u,
                sparse_tuning_selection_v = sparse_tuning_selection_v,
                GCV_score_smooth = GCV_score_smooth)
  return(result)
}



############################### Process bar indexing ###############################
ordinal <- function(i) {
  if (i == 1) {
    return(paste0(i, "st"))

  } else if (i == 2) {
    return(paste0(i, "nd"))

  } else if (i == 3) {
    return(paste0(i, "rd"))

  }
  else {
    return(paste0(i, "th"))

  }
}

