############################### Considering some values for alpha ###############################
get.Alphas = function(n=103,a=2,s=-100) {return(a^seq(s,s+n))}


############################### Calculating the optimal alpha using CV ###############################
opt_alpha <- function(X,
                      n_var,
                      ncol,
                      S_alphas_v,
                      S_alphas_u,
                      alphas,
                      sparse_tuning_result_u,
                      sparse_tuning_result_v,
                      sparse_tuning_type,
                      two_way_smoothness) {

  n_iter <- nrow(alphas)  # Update to get the number of iterations
  #pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  #                     max = n_iter, # Maximum value of the progress bar
  #                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
  #                     width = 50,   # Progress bar width. Defaults to getOption("width")
  #                     char = "=")   # Character used to create the bar

  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    #close(pb)
    return(list(GCV = Inf, opt.alpha = 0, opt_s.alpha = diag(n), GCVdf = data.frame(alphas, rep(Inf, n_iter))))
  } else {
    for (i in 1:n_iter) {
      S <- S_alphas_v[[i]]
      GCV_alpha <- 0
      m <- ncol(X)

      power_result <- power_algo(data = X,
                                 S_alpha = S,
                                 sparse_tuning_result_u = sparse_tuning_result_u,
                                 sparse_tuning_result_v = sparse_tuning_result_v,
                                 sparse_tuning_type = sparse_tuning_type)
      u <- power_result[[2]]
      GCV_alpha <- (1/m) * (norm_vec((diag(m) - S) %*% (t(X) %*% u))^2 / (norm_vec(u)^2 * (1 - (1/m) * sum(diag(S)))^2))

      #setTxtProgressBar(pb, i + (k - 1) / n_var)
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV))

    ########### Smoothness on u ###########
    if(two_way_smoothness != 0){
      GCV_u <- numeric(length(two_way_smoothness))
      for(j in 1:length(two_way_smoothness)){
        alpha_u <-  two_way_smoothness[j]
        S_u <-  S_alphas_u[[j]]
        GCV_alpha_u <- 0
        m <- ncol(t(X))

        power_result <- power_algo(data = t(X),
                                   S_alpha = S_u,
                                   sparse_tuning_result_u = sparse_tuning_result_u,
                                   sparse_tuning_result_v = sparse_tuning_result_v,
                                   sparse_tuning_type = sparse_tuning_type)
        u_u <- power_result[[2]]
        GCV_alpha_u <- (1/m) * (norm_vec((diag(m) - S_u) %*% (X %*% u_u))^2 / (norm_vec(u_u)^2 * (1 - (1/m) * sum(diag(S_u)))^2))
        GCV_u[j] <- GCV_alpha_u
      }
      opt.alpha_u <- two_way_smoothness[which.min(GCV_u), ]
      opt_s.alpha_u <- S_alphas_u[[which.min(GCV_u)]]

      result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV),
                     GCV_u = GCV_u, opt.alpha_u = opt.alpha_u, opt_s.alpha_u = opt_s.alpha_u, GCVdf_u = data.frame(two_way_smoothness, GCV_u))
    }

    return(result)
  }
}


############################ Conditional Tuning Parameters  - CV and GCV ############################
parameter_selection_conditional <- function(X_temp =  X_temp,
                                            n_var = n_var,
                                            ncol = ncol,
                                            n = n,
                                            smooth_tuning = smooth_tuning,
                                            sparse_tuning_u = sparsity_row_list, # vector
                                            sparse_tuning_v = sparsity_col_list, # list
                                            sparse_tuning_type = sparse_tuning_type,
                                            K_fold,
                                            S_alpha_list_v ,
                                            S_alpha_list_u,
                                            tuning_order){



  # tuning_order = 'Sparsity'
  if(tuning_order == "Sparsity"){
    CV_score_sparse_u <- CV_score_sparse_v <- GCV_score_smooth_u <- GCV_score_smooth_v <- Inf
    CV_scores_result_u <- CV_scores_result_u <- c()
    result = c()

    # Step 1: Sparsity on u  with no smoothness
    for (sparse_tuning_single in sparse_tuning_u) {
      sparse_score = cv_sparse_row(data=X_temp,
                                   S_alpha_v = diag(ncol(X_temp)),
                                   S_alpha_u = diag(nrow(X_temp)),
                                   K_fold = K_fold,
                                   sparse_tuning_result_u = sparse_tuning_single,
                                   sparse_tuning_result_v = 0,
                                   sparse_tuning_type = sparse_tuning_type,
                                   type = "CV")
      CV_scores_result_u <- c(CV_scores_result_u, sparse_score)
      if (sparse_score <= CV_score_sparse_u) {
        CV_score_sparse_u = sparse_score
        sparse_tuning_selection_u = sparse_tuning_single
      }
    }


    # Step 2: Conditional sparsity on v having sparsity on u with no smoothness
    for (sparse_tuning_single in sparse_tuning_v) {
      sparse_score = cv_score_sparse(data=X_temp,
                                     K_fold,
                                     sparse_tuning_result_v = sparse_tuning_single,
                                     sparse_tuning_result_u = sparse_tuning_selection_u,
                                     sparse_tuning_type,
                                     S_alpha_v = diag(ncol(X_temp)),
                                     S_alpha_u = diag(nrow(X_temp)),
                                     type = "CV") # Returns u only in the power func!

      if (sparse_score <= CV_score_sparse_v) {
        CV_score_sparse_v = sparse_score
        sparse_tuning_selection_v = sparse_tuning_single
        sparse_tuning_selection_u = 0
      }
    }








  }


















  # Step 1: Sparsity (No smoothness)
  CV_score_sparse_u <- CV_score_sparse_v <- GCV_score_smooth_u <- GCV_score_smooth_v <- Inf
  result = c()

  count = 0

  n_iter <- nrow(smooth_tuning) + length(sparse_tuning_u) + length(sparse_tuning_v)
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar


  ######  Sparsity on v (default)  ######
  # Sparsity tuning parameter using CV
  for (sparse_tuning_single in sparse_tuning_v) {
    count = count +1
    setTxtProgressBar(pb, count)


    sparse_score = cv_score_sparse(data=X_temp,
                                   K_fold,
                                   sparse_tuning_result_v = sparse_tuning_single,
                                   sparse_tuning_result_u = 0,
                                   sparse_tuning_type,
                                   S_alpha_v = diag(ncol(X_temp)),
                                   S_alpha_u = diag(nrow(X_temp)),
                                   type = "CV") # Returns u only in the power func!

    if (sparse_score <= CV_score_sparse_v) {
      CV_score_sparse_v = sparse_score
      sparse_tuning_selection_v = sparse_tuning_single
      sparse_tuning_selection_u = 0
    }
  }


  ###### Sparsity on u  ######
  if (two_way_sparsity == TRUE){

    # Sparsity tuning parameter using CV
    for (sparse_tuning_single in sparse_tuning_u) {
      count = count +1
      setTxtProgressBar(pb, count)

      sparse_score = cv_score_sparse(data=t(hd),
                                     K_fold,
                                     sparse_tuning_result_u = sparse_tuning_single,
                                     sparse_tuning_result_v = 0,
                                     sparse_tuning_type,
                                     S = diag(ncol(data)), # No Smoothness
                                     type = "CV")

      if (sparse_score <= CV_score_sparse_u) {
        CV_score_sparse_u = sparse_score
        sparse_tuning_selection_u = sparse_tuning_single
      }
    }
  }


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






