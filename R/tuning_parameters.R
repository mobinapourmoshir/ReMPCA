############################ CV Scores calculators for sparsity and smoothness ############################
cv_score_sparse <- function(data, S, K_fold, sparse_tuning_single, sparse_tuning_type, shuffled_row, group_size) {
  data_tilde <- (data)  # Group the rows of data
  error_score_sparse <- 0

  for (k in 1:K_fold) {
    rows_to_remove <- shuffled_row[((k - 1) * group_size + 1):(k * group_size)]
    data_train <- data_tilde[-rows_to_remove, ]  # X^{-k}
    data_test <- data_tilde[rows_to_remove, ]    # X^k

    # Ensure u_test is a column vector with the same number of rows as columns in data_train
    u_test <- power_algo(t(data_train), sparse_tuning_result = sparse_tuning_single, sparse_tuning_type, S_alpha = S, type = "CV") # Returns u only!

    # Ensure data_test has the same number of columns as the length of u_test
    v_test <- data_test %*% as.matrix(u_test)

    # Ensure the dimensions of data_test and the reconstructed data match
    reconstructed <- as.matrix(data_test - t(as.matrix(u_test) %*% as.matrix(t(v_test))))
    error_score_sparse <- error_score_sparse + (norm_vec(reconstructed)^2)
  }

  return(error_score_sparse / ncol(data))  # Assuming ncol(data) is the number of grid points N
}

############################### Considering some values for alpha ###############################
get.Alphas = function(n=103,a=2,s=-100) {return(a^seq(s,s+n))}


############################### Calculating the optimal alpha using CV ###############################
opt_alpha <- function(X, nvar, ncol, S_alphas, alphas, CV_sparse_tuning_result ,sparse_tuning_type) {

  n_iter <- nrow(alphas)  # Update to get the number of iterations
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    close(pb)
    return(list(GCV = Inf, opt.alpha = 0, opt_s.alpha = diag(n), GCVdf = data.frame(alphas, GCV)))
  } else {
    for (i in 1:n_iter) {
      S <- S_alphas[[i]]
      GCV_alpha <- 0
      df <- X

      for (k in 1:nvar) {
        m <- as.integer(ncol[k]) # m_k
        Xk <- df[, 1:m]
        df <- df[, -(1:m), drop = FALSE]
        Sk <- S[1:m, 1:m]
        S <- S[-(1:m), -(1:m), drop = FALSE]

        power_result <- power_algo(data = Xk, S_alpha = Sk, sparse_tuning_result = CV_sparse_tuning_result ,sparse_tuning_type = sparse_tuning_type)
        u <- power_result[[2]]
        GCV_alpha <- GCV_alpha + (1/m) * (norm_vec((diag(m) - Sk) %*% (t(Xk) %*% u))^2 / (norm_vec(u)^2 * (1 - (1/m) * sum(diag(Sk)))^2))

        # Update progress bar within the inner loop
        setTxtProgressBar(pb, i + (k - 1) / nvar)
      }

      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas[[which.min(GCV)]]
    close(pb)

    return(list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV)))
  }
}


############################ Conditional Tuning Parameters  - CV and GCV ############################
parameter_selection_conditional <- function(data, nvar, ncol, n , smooth_tuning, sparse_tuning, sparse_tuning_type, K_fold = 5,
                                            S_alpha_List_v ,S_alpha_list_u , two_way_smoothness, two_way_sparsity,
                                            sparse_tuning_u, sparse_tuning_v){

  CV_score_sparse_u <- CV_score_sparse_v <- GCV_score_smooth_u <- GCV_score_smooth_v <- 10^60
  result = c()

  count = 0

  n_iter <- nrow(smooth_tuning) + length(sparse_tuning)
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar


  ######  Sparsity on u  ######
  shuffled_row = sample(nrow(data)) # Grouping the rows of data matrix
  group_size = length(shuffled_row) / K_fold


  # Sparsity tuning parameter using CV
  for (sparse_tuning_single in sparse_tuning_u) {
    count = count +1
    setTxtProgressBar(pb, count)
    if (sparse_tuning_single == 0) {
      sparse_score = 0
    } else{
      sparse_score = cv_score_sparse(data=data,K_fold,sparse_tuning_single,sparse_tuning_type,shuffled_row,group_size,S = diag(ncol(data)))
    }
    if (sparse_score <= CV_score_sparse_u) {
      CV_score_sparse_u = sparse_score
      sparse_tuning_selection_u = sparse_tuning_single
    }
  }


  ###### Sparsity on v  ######
  if (two_way_sparsity == TRUE){

    shuffled_col = sample(ncol(data)) # Grouping the columns of data matrix
    group_size_v = length(shuffled_col) / K_fold


    # Sparsity tuning parameter using CV
    for (sparse_tuning_single in sparse_tuning_v) {
      count = count +1
      setTxtProgressBar(pb, count)
      if (sparse_tuning_single == 0) {
        sparse_score = 0
      } else{
        sparse_score = cv_score_sparse(data=data,K_fold,sparse_tuning_single,sparse_tuning_type,shuffled_col,group_size,S = diag(nrow(data)))
      }
      if (sparse_score <= CV_score_sparse_v) {
        CV_score_sparse_v = sparse_score
        sparse_tuning_selection_v = sparse_tuning_single
      }
    }
  }


  ############## For Functional data only! ##############
  ######  Smoothness on v  ######
  # Smoothing tuning parameter using GCV
  GCV_score_smooth_v = opt_alpha(X = data , nvar = nvar, ncol = ncol, S_alphas = S_alpha_list_v , alphas = smooth_tuning ,
                                CV_sparse_tuning_result = sparse_tuning_selection_u, sparse_tuning_type = sparse_tuning_type)


  ######  Smoothness on u  ######
  if(two_way_smoothness == TRUE){
    GCV_score_smooth_u = opt_alpha(X = t(data) , nvar = nvar, ncol = n, S_alphas = S_alpha_list_u , alphas = smooth_tuning ,
                                   CV_sparse_tuning_result = sparse_tuning_selection_u, sparse_tuning_type = sparse_tuning_type)
  }


  close(pb) # Close the connection

  if(two_way_smoothness == FALSE & two_way_sparsity == FALSE){
    result = list(sparse_tuning_selection_u, GCV_score_smooth_v)
  }else if(two_way_smoothness == TRUE & two_way_sparsity == FALSE){
    result = list(sparse_tuning_selection_u , GCV_score_smooth_u , GCV_score_smooth_v)
  }else if(two_way_smoothness == FALSE & two_way_sparsity == TRUE){
    result = list(sparse_tuning_selection_u  ,sparse_tuning_selection_v , GCV_score_smooth_v)
  }else if(two_way_smoothness == TRUE & two_way_sparsity == TRUE){
    result = list(sparse_tuning_selection_u  ,sparse_tuning_selection_v , GCV_score_smooth_v, GCV_score_smooth_u)
  }

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


