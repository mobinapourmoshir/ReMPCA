############################ CV Scores calculators for sparsity and smoothness ############################
cv_score_sparse <- function(data,
                            S,
                            K_fold,
                            sparse_tuning_result_u,
                            sparse_tuning_result_v,
                            sparse_tuning_type,
                            type) {
  set.seed(123)
  shuffled_row <- sample(ncol(data)) # Grouping the rows of data matrix
  group_size <- ifelse(round(length(shuffled_row)/ K_fold,
                             digits = 0) == 0, 1,
                       round(length(shuffled_row)/ K_fold, digits = 0))

  data_tilde <- t(data)  # Group the rows of data
  error_score_sparse <- 0



  for (k in 1:K_fold) {
    rows_to_remove <- shuffled_row[((k - 1) * group_size + 1):((k) * group_size)]
    data_train <- data.frame(data_tilde[-rows_to_remove, ])  # X^{-k}
    data_test <- data.frame(matrix(data_tilde[rows_to_remove, ], nrow = length(rows_to_remove)))    # X^k

    # Ensure u_test is a column vector with the same number of rows as columns in data_train
    u_test <- power_algo(data = t(data_train),
                         sparse_tuning_result_u = sparse_tuning_result_u,
                         sparse_tuning_result_v = sparse_tuning_result_v,
                         sparse_tuning_type = sparse_tuning_type,
                         S_alpha = S,
                         type = type)

    rownames(u_test) <- NULL
    colnames(u_test) <- NULL
    # Ensure data_test has the same number of columns as the length of u_test
    v_test <- as.matrix(data_test)%*%as.matrix(u_test)


    data_test_back = t(data_tilde)[, rows_to_remove]
    error_score_sparse = error_score_sparse + sum((
      t(data_test_back) - v_test %*% t(u_test)) ^ 2)
  }

  return(error_score_sparse / ncol(data))  # Assuming ncol(data) is the number of grid points N
}

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
                                            n_var,
                                            ncol,
                                            n,
                                            smooth_tuning,
                                            sparse_tuning_u,
                                            sparse_tuning_v,
                                            sparse_tuning_type,
                                            K_fold,
                                            S_alpha_list_v ,
                                            S_alpha_list_u ,
                                            two_way_smoothness,
                                            two_way_sparsity){


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



  # Splitting X_temp
  start_col <- 1
  for (i in 1:n_var) {
    num_cols <- ncol[[i]]  # Get number of columns for the ith matrix
    sub_matrix <- X_temp[, start_col:(start_col + num_cols - 1)]  # Extract matrix
    start_col <- start_col + num_cols

  }


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
                                   S = diag(nrow(X_temp)), # No Smoothness
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

