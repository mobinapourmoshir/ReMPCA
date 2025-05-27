################## Conditional Tuning Parameters  - CV and GCV #################
parameter_selection <- function(X_temp,
                                n_var,
                                ncol,
                                n,
                                GridPoints_u,
                                GridPoints_v,
                                smooth_tuning_v,  # matrix
                                smooth_tuning_u,  # vector
                                sparse_tuning_u,  # vector
                                sparse_tuning_v,  # list
                                sparse_tuning_type,
                                nfolds_u,
                                nfolds_v,
                                S_alpha_list_v,
                                S_alpha_list_u,
                                Omegas_u,
                                Omegas_v,
                                tuning_iter,
                                tuning_order,
                                thresh,
                                maxit,
                                cv.pick = "1se",
                                smoothness_type) {

  alpha_u <- 0
  alpha_v <- rep(0, n_var)
  gamma_u <- 0
  gamma_v <- rep(0, n_var)

  all_alphas_u <- numeric(tuning_iter)
  all_gamma_u <- numeric(tuning_iter)
  all_alphas_v <- matrix(0, nrow = tuning_iter, ncol = n_var)
  all_gamma_v <- matrix(0, nrow = tuning_iter, ncol = n_var)

  for (iter in 1:tuning_iter) {
    if (tuning_order == "Sparsity") {

      # Progress bar
      total_steps <- nrow(smooth_tuning_v) + length(smooth_tuning_u) +
        length(sparse_tuning_u) + sum(lengths(sparse_tuning_v))
      pb <- txtProgressBar(min = 0,
                           max = total_steps,
                           width = 50,
                           style = 3)
      step <- 0

      # S_alpha for v
      S_alpha_v0 <- Omega_v0 <- list()
      for (i in 1:n_var) {
        tds <- GridPoints_v[[i]]
        if(is.null(tds)){
          S_alpha_v0[[i]] <- Omega_v0[[i]] <- diag(ncol[,i])
        }else{
          getpenresult <- get.pen(td = tds,
                                  alpha = as.numeric(alpha_v[i]),
                                  type = smoothness_type)
          S_alpha_v0[[i]] <- getpenresult$S.alpha
          Omega_v0[[i]] <- getpenresult$Omega
        }
      }

      S_alpha_u0 <- Omega_u0 <- matrix()
      tds <- GridPoints_u
      if(is.null(tds)){
        S_alpha_u0 <- Omega_u0 <- diag(n)
      } else{
        getpenresult <- get.pen(td = tds,
                                alpha = as.numeric(alpha_u),
                                type = smoothness_type)
        S_alpha_u0 <- getpenresult$S.alpha
        Omega_u0 <- getpenresult$Omega
      }

      # Sparsity on u
      cv_row_result <- cv_sparse_row(data = X_temp,
                                     n_var = n_var,
                                     ncol = ncol,
                                     S_alpha_u = S_alpha_u0,
                                     S_alpha_v = S_alpha_v0,
                                     K_fold = nfolds_u,
                                     thresh = thresh,
                                     maxit = maxit,
                                     conditional = FALSE,
                                     cv.pick = cv.pick,
                                     sparse_tuning_result_u = sparse_tuning_u,
                                     sparse_tuning_result_v = gamma_v,
                                     sparse_tuning_type = sparse_tuning_type)
      gamma_u <- cv_row_result[[1]]
      last_cv_result_u <- cv_row_result[[2]]
      step <- step + 1; setTxtProgressBar(pb, step)

      CV_results <- CV_scores_result_v <- list()
      gamma_v <- rep(0, n_var)

      # Sparsity on v
      for (i in 1:n_var) {
        K_fold <- nfolds_v[i]
        gamma_Xi <- sparsity_col_list[[i]]
        cv_means <- cv_ses  <- numeric(length(gamma_Xi))

        for (j in seq_along(gamma_Xi)) {
          gamma_v[i] <- gamma_Xi[j]

          cv_result <- cv_sparse_col(data = X_temp,
                                     n_var = n_var,
                                     ncol = ncol,
                                     thresh = thresh,
                                     maxit = maxit,
                                     conditional = FALSE,
                                     S_alpha_v = S_alpha_v0,
                                     S_alpha_u = S_alpha_u0,
                                     K_fold = K_fold,
                                     sparse_tuning_result_u = gamma_u,
                                     sparse_tuning_result_v = gamma_v,
                                     sparse_tuning_type = sparse_tuning_type)

          cv_means[j] <- cv_result$CV_error
          cv_ses[j] <- cv_result$SE *sqrt(K_fold)
          step <- step + 1
          setTxtProgressBar(pb, step)
        }

        CV_results[[i]] <- data.frame(gamma_Xi, cv_means, cv_ses)

        CV_scores_result_v[[i]] <- list(errors = cv_means, SEs = cv_ses)

        if (cv.pick == "1se") {
          j_min <- which.min(cv_means)
          threshold <- cv_means[j_min] + cv_ses[j_min]
          j_1se <- max(which(cv_means <= threshold))
          sparse_tuning_selection_v <- gamma_Xi[j_1se]
        } else if (cv.pick == "min") {
          j_min <- which.min(cv_means)
          sparse_tuning_selection_v <- gamma_Xi[j_min]
        } else {
          stop("cv.pick must be either '1se' or 'min'")
        }

        gamma_v[i] <- sparse_tuning_selection_v
      }
      last_cv_result_v <- CV_results
      step <- step + 1; setTxtProgressBar(pb, step)

      # Smoothness on u
      opt_u <- opt_alpha_u(X = X_temp,
                           n_var = n_var,
                           ncol = ncol,
                           thresh = thresh,
                           maxit = maxit,
                           conditional = TRUE,
                           S_alphas_v = S_alpha_v0,
                           S_alphas_u = S_alpha_list_u,
                           alphas_u = smooth_tuning_u,
                           alpha_v = alpha_v,
                           Omega_v = Omega_v0,
                           Omegas_u = Omegas_u,
                           sparse_tuning_result_u = gamma_u,
                           sparse_tuning_result_v = gamma_v,
                           sparse_tuning_type = sparse_tuning_type)
      alpha_u <- opt_u$opt.alpha_u
      last_gcv_result_u <- opt_u$GCVdf
      last_opt_S_u <- opt_u$opt_s.alpha_u
      last_opt_alpha_omega_u <- opt_u$opt_alpha_Omega_u
      step <- step + 1; setTxtProgressBar(pb, step)

      # Smoothness on v
      opt_v <- opt_alpha_v(X = X_temp,
                           n_var = n_var,
                           ncol = ncol,
                           thresh = thresh,
                           maxit = maxit,
                           conditional = TRUE,
                           S_alphas_v = S_alpha_list_v,
                           S_alphas_u = opt_u$opt_s.alpha_u,
                           alphas_v = smooth_tuning_v,
                           alpha_u = alpha_u,
                           Omega_u = opt_u$Omega_u,
                           Omegas_v = Omegas_v,
                           sparse_tuning_result_u = gamma_u,
                           sparse_tuning_result_v = gamma_v,
                           sparse_tuning_type = sparse_tuning_type)
      alpha_v <- opt_v$opt.alpha_v
      last_gcv_result_v <- opt_v$GCVdf
      last_opt_S_v <- opt_v$opt_s.alpha_v
      last_opt_alpha_omega_v <- opt_u$opt_alpha_Omega_v
      step <- step + 1; setTxtProgressBar(pb, step)

      close(pb)
      cat("Progress completed:", step, "of", total_steps, "\n")
      stopifnot(step == total_steps)

    } else if (tuning_order == "Smoothness") {

      # Progress bar
      total_steps <- nrow(smooth_tuning_v) + length(smooth_tuning_u) +
        length(sparse_tuning_u) + sum(lengths(sparse_tuning_v))
      pb <- txtProgressBar(min = 0,
                           max = total_steps,
                           width = 50,
                           style = 3)
      step <- 0

      # Current S_alpha_v
      S_alpha_v0 <- Omega_v0 <- list()
      for (i in 1:n_var) {
        tds <- GridPoints_v[[i]]
        if(is.null(tds)){
          S_alpha_v0[[i]] <- Omega_v0[[i]] <- diag(ncol[,i])
        }else{
          getpenresult <- get.pen(td = tds,
                                  alpha = as.numeric(alpha_v[i]),
                                  type = smoothness_type)
          S_alpha_v0[[i]] <- getpenresult$S.alpha
          Omega_v0[[i]] <- getpenresult$Omega
        }
      }

      # Current S_alpha_u
      S_alpha_u0 <- Omega_u0 <- matrix()
      tds <- GridPoints_u
      if(is.null(tds)){
        S_alpha_u0 <- Omega_u0 <- diag(n)
      }else{
        getpenresult <- get.pen(td = tds,
                                alpha = as.numeric(alpha_u),
                                type = smoothness_type)
        S_alpha_u0 <- getpenresult$S.alpha
        Omega_u0 <- getpenresult$Omega
      }

      # Smoothness on u
      opt_u <- opt_alpha_u(X = X_temp,
                           n_var = n_var,
                           ncol = ncol,
                           thresh = thresh,
                           maxit = maxit,
                           conditional = TRUE,
                           S_alphas_v = S_alpha_v0,
                           S_alphas_u = S_alpha_list_u,
                           alphas_u = smooth_tuning_u,
                           alpha_v = alpha_v,
                           Omega_v = Omega_v0,
                           Omegas_u = Omegas_u,
                           sparse_tuning_result_u = gamma_u,
                           sparse_tuning_result_v = gamma_v,
                           sparse_tuning_type = sparse_tuning_type)

      alpha_u <- opt_u$opt.alpha_u
      last_gcv_result_u <- opt_u$GCVdf
      last_opt_S_u <- opt_u$opt_s.alpha_u
      last_opt_alpha_omega_u <- opt_u$opt_alpha_Omega_u
      step <- step + 1; setTxtProgressBar(pb, step)

      # Smoothness on v
      opt_v <- opt_alpha_v(X = X_temp,
                           n_var = n_var,
                           ncol = ncol,
                           thresh = thresh,
                           maxit = maxit,
                           conditional = TRUE,
                           S_alphas_v = S_alpha_list_v,
                           S_alphas_u = opt_u$opt_s.alpha_u,
                           alphas_v = smooth_tuning_v,
                           alpha_u = alpha_u,
                           Omega_u = opt_u$Omega_u,
                           Omegas_v = Omegas_v,
                           sparse_tuning_result_u = gamma_u,
                           sparse_tuning_result_v = gamma_v,
                           sparse_tuning_type = sparse_tuning_type)
      alpha_v <- opt_v$opt.alpha_v
      last_gcv_result_v <- opt_v$GCVdf
      last_opt_S_v <- opt_v$opt_s.alpha_v
      last_opt_alpha_omega_v <- opt_v$opt_alpha_Omega_v
      step <- step + 1; setTxtProgressBar(pb, step)

      # Sparsity on u
      cv_row_result <- cv_sparse_row(data = X_temp,
                                     n_var = n_var,
                                     ncol = ncol,
                                     S_alpha_u = opt_u$opt_s.alpha_u,
                                     S_alpha_v = opt_v$opt_s.alpha_v,
                                     K_fold = nfolds_u,
                                     thresh = thresh,
                                     maxit = maxit,
                                     conditional = FALSE,
                                     cv.pick = cv.pick,
                                     sparse_tuning_result_u = sparse_tuning_u,
                                     sparse_tuning_result_v = gamma_v,
                                     sparse_tuning_type = sparse_tuning_type)
      gamma_u <- cv_row_result[[1]]
      last_cv_result_u <- cv_row_result[[2]]
      step <- step + 1; setTxtProgressBar(pb, step)

      # Sparsity on v
      CV_scores_result_v <- CV_results <- list()
      gamma_v <- rep(0, n_var)

      for (i in 1:n_var) {
        K_fold <- nfolds_v[i]
        gamma_Xi <- sparsity_col_list[[i]]
        cv_means <- cv_ses  <- numeric(length(gamma_Xi))

        for (j in seq_along(gamma_Xi)) {
          gamma_v[i] <- gamma_Xi[j]

          cv_result <- cv_sparse_col(data = X_temp,
                                     n_var = n_var,
                                     ncol = ncol,
                                     thresh = thresh,
                                     maxit = maxit,
                                     S_alpha_v = opt_v$opt_s.alpha_v,
                                     S_alpha_u = opt_u$opt_s.alpha_u,
                                     K_fold = K_fold,
                                     sparse_tuning_result_u = gamma_u,
                                     sparse_tuning_result_v = gamma_v,
                                     sparse_tuning_type = sparse_tuning_type)

          cv_means[j] <- cv_result$CV_error
          cv_ses[j] <- cv_result$SE *sqrt(K_fold)
          step <- step + 1
          setTxtProgressBar(pb, step)
        }

        CV_results[[i]] <- data.frame(gamma_Xi, cv_means, cv_ses)

        # Store results
        CV_scores_result_v[[i]] <- list(errors = cv_means, SEs = cv_ses)

        if (cv.pick == "1se") {
          j_min <- which.min(cv_means)
          threshold <- cv_means[j_min] + cv_ses[j_min]
          j_1se <- max(which(cv_means <= threshold))  # Most sparse within 1-SE
          sparse_tuning_selection_v <- gamma_Xi[j_1se]
        } else if (cv.pick == "min") {
          j_min <- which.min(cv_means)
          sparse_tuning_selection_v <- gamma_Xi[j_min]
        } else {
          stop("cv.pick must be either '1se' or 'min'")
        }

        gamma_v[i] <- sparse_tuning_selection_v
      }
      last_cv_result_v <- CV_results
      step <- step + 1; setTxtProgressBar(pb, step)

      close(pb)
      cat("Progress completed:", step, "of", total_steps, "\n")
      stopifnot(step == total_steps)
    }

    # Store all iterations
    all_alphas_u[iter] <- alpha_u
    all_gamma_u[iter] <- gamma_u
    all_alphas_v[iter, ] <- alpha_v
    all_gamma_v[iter, ] <- gamma_v
  }

  # Output of the function
  final_results <- list(
    sparse_tuning_selection_u = gamma_u,
    sparse_tuning_selection_v = gamma_v,
    smooth_tuning_selection_u = alpha_u,
    smooth_tuning_selection_v = as.numeric(alpha_v),
    all_alphas_u = all_alphas_u,
    all_gamma_u = all_gamma_u,
    all_alphas_v = all_alphas_v,
    all_gamma_v = all_gamma_v,
    last_cv_result_u = last_cv_result_u,
    last_cv_result_v = last_cv_result_v,
    last_gcv_result_u = last_gcv_result_u,
    last_gcv_result_v = last_gcv_result_v,
    last_opt_S_u = last_opt_S_u,
    last_opt_S_v = last_opt_S_v,
    last_opt_alpha_omega_u = last_opt_alpha_omega_u,
    last_opt_alpha_omega_v =last_opt_alpha_omega_v
  )

  return(final_results)
}

########################### Process bar indexing ###########################
ordinal <- function(i) {
  if (i == 1) {
    return(paste0(i, "st"))
  } else if (i == 2) {
    return(paste0(i, "nd"))
  } else if (i == 3) {
    return(paste0(i, "rd"))
  } else {
    return(paste0(i, "th"))
  }
}
