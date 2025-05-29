########## CV Scores calculators for sparsity on rows (PC scores, u) ##########
cv_sparse_row <- function(data,
                          n_var,
                          ncol,
                          S_alpha_u,
                          S_alpha_v,
                          K_fold,
                          cv.pick,
                          thresh,
                          maxit,
                          parallel = FALSE,
                          cl       = NULL,
                          conditional = FALSE,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v,
                          sparse_tuning_type) {

  # 1) pre‐compute things that don't depend on j:
  set.seed(123)
  shuffled_cols <- sample(ncol(data))
  group_size    <- ceiling(length(shuffled_cols) / K_fold)

  # 2) helper that does the exact same inner logic for a single j:
  compute_for_j <- function(j) {
    gamma_u    <- sparse_tuning_result_u[j]
    fold_errs  <- numeric(K_fold)

    for (k in seq_len(K_fold)) {
      cols_out   <- shuffled_cols[((k-1)*group_size + 1):
                                    min(k*group_size, length(shuffled_cols))]
      data_train <- data[, -cols_out, drop=FALSE]
      data_test  <- data[,  cols_out, drop=FALSE]
      up_ncol    <- update_ncol(ncol, cols_out)

      # rebuild S_alpha_v for train
      ncol_vec        <- as.numeric(ncol[1,])
      cum_ncol        <- c(0, cumsum(ncol_vec))
      S_alpha_v_train <- lapply(1:n_var, function(i) {
        cols_i    <- (cum_ncol[i]+1):cum_ncol[i+1]
        keep      <- setdiff(cols_i, cols_out)
        idx_local <- match(keep, cols_i)
        S_alpha_v[[i]][ idx_local, idx_local, drop=FALSE ]
      })

      pr <- power_algo(
        data                    = data_train,
        n_var                   = n_var,
        ncol                    = up_ncol,
        thresh                  = thresh,
        maxit                   = maxit,
        conditional             = conditional,
        sparse_tuning_result_u  = gamma_u,
        sparse_tuning_result_v  = sparse_tuning_result_v,
        S_alpha_v               = S_alpha_v_train,
        S_alpha_u               = S_alpha_u,
        sparse_tuning_type      = sparse_tuning_type
      )

      u_hat      <- pr[[2]]
      v_pred     <- t(as.matrix(data_test)) %*% u_hat
      fold_errs[k] <- sum((data_test - u_hat %*% t(v_pred))^2) / ncol(data_test)
    }

    # return named vector with error & se
    err <- mean(fold_errs)
    se  <- sd(fold_errs) / sqrt(K_fold)
    list(err = err, se = se, folds = fold_errs)
  }

  # 3) pick the right apply‐fun
  if (parallel) {
    if (is.null(cl))
      stop("When parallel=TRUE you must pass a cluster 'cl'")
    # export only the binding names your compute_for_j needs:
    parallel::clusterExport(cl,
                            varlist = c("data","n_var","ncol","S_alpha_u","S_alpha_v",
                                        "K_fold","thresh","maxit","conditional",
                                        "sparse_pen_fun", "norm_vec",
                                        "sparse_tuning_result_v","sparse_tuning_type",
                                        "update_ncol","power_algo"),
                            envir = environment()
    )
    applyFun <- function(X, FUN) parallel::parLapplyLB(cl, X, FUN)
  } else {
    applyFun <- lapply
  }

  # 4) run over all j in parallel or serial
  out_list <- applyFun(seq_along(sparse_tuning_result_u), compute_for_j)

  # 5) pull out CV_errors, SE_errors, and fold_errors_list
  CV_errors       <- sapply(out_list, `[[`, "err")
  SE_errors       <- sapply(out_list, `[[`, "se")
  fold_errors_list<- lapply(out_list, `[[`, "folds")

  # 6) pick gamma_u as before
  if (cv.pick == "1se") {
    j0      <- which.min(CV_errors)
    cutoff  <- CV_errors[j0] + SE_errors[j0]
    j_star  <- max(which(CV_errors <= cutoff))
  } else if (cv.pick == "min") {
    j_star <- which.min(CV_errors)
  } else {
    stop("cv.pick must be either '1se' or 'min'")
  }
  gamma_u <- sparse_tuning_result_u[j_star]

  # 7) return exactly what you did before
  list(
    gamma_u,
    cv_results = data.frame(
      sparse_tuning_result_u,
      CV_errors,
      SE_errors
    )
  )
}

########### CV Scores calculators for sparsity on columns (PCs, v) ###########
cv_sparse_col <- function(data,
                          n_var,
                          ncol,
                          S_alpha_v,
                          S_alpha_u,
                          K_fold,
                          thresh,
                          maxit,
                          parallel = FALSE,
                          cl       = NULL,
                          conditional = FALSE,
                          sparse_tuning_result_u,
                          sparse_tuning_result_v,
                          sparse_tuning_type) {

  # 1) prepare fold assignments
  set.seed(123)
  shuffled_rows <- sample(nrow(data))
  group_size    <- ifelse(
    round(length(shuffled_rows) / K_fold, 0) == 0,
    1,
    round(length(shuffled_rows) / K_fold, 0)
  )

  # 2) helper: compute the error for a single fold k
  compute_fold <- function(k) {
    # which rows go to test
    rows_out   <- shuffled_rows[((k - 1) * group_size + 1):
                                  min(k * group_size, length(shuffled_rows))]
    # train / test split
    train_df   <- data.frame(data[-rows_out, , drop = FALSE])
    test_df    <- data.frame(data[ rows_out, , drop = FALSE])

    # run your exact same power‐algorithm call
    pr <- power_algo(
      data                    = train_df,
      n_var                   = n_var,
      ncol                    = ncol,
      conditional             = conditional,
      thresh                  = thresh,
      maxit                   = maxit,
      sparse_tuning_result_u  = sparse_tuning_result_u,
      sparse_tuning_result_v  = sparse_tuning_result_v,
      S_alpha_v               = S_alpha_v,
      S_alpha_u               = S_alpha_u[-rows_out, -rows_out],
      sparse_tuning_type      = sparse_tuning_type
    )

    v_train <- pr[[1]]
    u_test  <- as.matrix(test_df) %*% as.matrix(v_train)

    # exactly the same fold‐error
    sum((test_df - u_test %*% t(v_train))^2) / nrow(test_df)
  }

  # 3) pick apply‐function
  if (parallel) {
    if (is.null(cl)) stop("Must supply a cluster 'cl' when parallel=TRUE")
    # export needed objects
    parallel::clusterExport(cl,
                            varlist = c("data","n_var","ncol","S_alpha_v","S_alpha_u",
                                        "thresh","maxit","conditional",
                                        "sparse_tuning_result_u","sparse_tuning_result_v",
                                        "sparse_tuning_type","power_algo"),
                            envir = environment()
    )
    applyFun <- function(X, FUN) parallel::parLapplyLB(cl, X, FUN)
  } else {
    applyFun <- lapply
  }

  # 4) run folds (in parallel or serial)
  fold_errors_list <- applyFun(seq_len(K_fold), compute_fold)
  fold_errors      <- unlist(fold_errors_list)

  # 5) compute CV summary
  CV_mean <- mean(fold_errors)
  CV_se   <- sd(fold_errors) / sqrt(K_fold)

  # 6) return exactly as before
  list(
    CV_error   = CV_mean,
    SE         = CV_se,
    fold_errors = fold_errors
  )
}


############# Computing number of columns when eliminating some ##############
update_ncol <- function(ncol_matrix, cols_to_remove) {
  col_counts <- as.numeric(ncol_matrix[1, ])
  cumulative_cols <- c(0, cumsum(col_counts))
  updated_ncol <- col_counts
  # Loop through each dataset and count removed columns
  for (i in seq_along(updated_ncol)) {
    start_idx <- cumulative_cols[i] + 1
    end_idx <- cumulative_cols[i + 1]

    # Count how many columns to remove in this dataset
    removed_count <- sum(cols_to_remove >= start_idx & cols_to_remove <= end_idx)
    updated_ncol[i] <- updated_ncol[i] - removed_count
  }
  return(as.data.frame(t(updated_ncol)))
}
