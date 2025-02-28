#' ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis
#'
#' @param object_list An hdClass object.
#' @param centerfns A logical; if True, it demeans the data before calculating the principal components.
#' @param num_pcs An integer. The number of principal components.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v, must be one of "Second_order" (default), "First_order" or "Indicator".
#' @param K_fold An integer. It's used in cross validation approach for tuning the level of sparsity.
#' @param tuning_order A character string representing the tuning order. If set to 'Sparsity', sparsity parameters are tuned first, followed by smoothness. If set to 'Smoothness', the order is reversed.
#'
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return PC scores, PC functions, ...
#' @export
#'


############################ Smooth and Sparse Multivariate PCA ############################

ReMPCA <- function(object_list,
                   centerfns = TRUE,
                   num_pcs = 1,
                   smoothness_type = "Second_order",
                   sparse_tuning_type = "soft",
                   K_fold = 5,
                   tuning_order = "Sparsity") {

  # Combine matrices side by side
  hd <- object_list
  n <- nrow(hd)
  n_var <- attr(object_list, "n_var") # Number of variables (# of matrices in object_list)
  ncol <- attr(object_list, "ncol") # Number of columns of each matrix

  ####### Smoothing Parameter ##########
  # Generate all combinations alphas (one row per combination)
  smooth_tuning_col <- expand.grid(attr(object_list, "Smoothing_parameter_col"))
  smooth_tuning_row <- attr(object_list, "Smoothing_parameter")

  ####### level of sparsity (for both functional and non-functional data) #######
  sparsity_row_list <- attr(object_list, "Sparsity_parameter")
  sparsity_col_list <- attr(object_list, "Sparsity_parameter_col")

  ####### Pre-processing: Centralizing the data #######
  X <- data.frame()
  if (centerfns) {
    X <- apply(hd, 2, function(x) x - mean(x))
  }else{
    X <- hd
  }

  ####### Grid Points #######
  GridPoints_u <- attr(object_list, "GridPoints_u")
  GridPoints_v <- attr(object_list, "GridPoints_v")

  ####### S_alpha for all alphas #######
  S_alpha_list_u <- S_alpha_list_v <- list()
  index <- 0
  cat("Preprocessing ...\n")
  n_iter1 <- nrow(smooth_tuning)     # The number of alphas
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter1,# Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  # S_alpha for v
  for (alpha_index in 1:nrow(smooth_tuning)) {
    index <- index + 1
    S_u <- S_v <- list()
    for (i in 1:n_var) {
      alpha <- as.numeric(smooth_tuning[alpha_index,i])
      if(is.null(GridPoints_v[[i]])){
        S_v[[i]] <- diag(ncol[,i])
      }else{
        S_v[[i]] <- get.pen(td = GridPoints_v[[i]],
                            alpha = alpha,
                            type = smoothness_type)
      }
    }
    S_alpha_list_v[[index]] <- as.matrix(bdiag(S_v))
    setTxtProgressBar(pb, index)
  }

  # S_alpha for u
  if(!is.null(Smoothing_parameter)){
    for (alpha_index in 1:length(smooth_tuning_row)) {
      alpha <- as.numeric(smooth_tuning_row[alpha_index])
      S_alpha_list_u[[alpha_index]] <- get.pen(td = GridPoints_u,
                                               alpha = alpha,
                                               type = smoothness_type)
    }
  }
  close(pb)

  # Initializing the lists
  opt_S_v <- opt_S_u <- list()
  smooth_tuning_result_v <- smooth_tuning_result_u <- list()
  GCV_v <- GCV_u <- list()
  GCVdf_v <- GCVdf_u <- list()
  sparse_tuning_result_u <- sparse_tuning_result_v <- list()
  lsv <- lsu <- c()
  funcs <- list()

  ####### ReMPCA Implementation #######
  for (j in 1:num_pcs) {
    cat(sprintf("Computing the %s PC ...\n", ordinal(j)))
    if (j == 1) {
      X_temp = X
    } else{
      SVD_result = svd(X_temp)
      v_original = SVD_result$v[,1]
      u_original = SVD_result$u[,1]
      sigma = SVD_result$d[1]
      X_temp = X_temp - sigma * u_original%*%t(v_original)
    }

    # Tuning Parameters
    opt_parameters_result <- opt_alpha_result <- list()
    opt_parameters_result <- parameter_selection_conditional(X_temp =  X_temp,
                                                             n_var = n_var,
                                                             ncol = ncol,
                                                             n = n,
                                                             smooth_tuning = smooth_tuning,
                                                             sparse_tuning_u = sparsity_row_list,
                                                             sparse_tuning_v = sparsity_col_list,
                                                             sparse_tuning_type = sparse_tuning_type,
                                                             K_fold,
                                                             S_alpha_list_v ,
                                                             S_alpha_list_u,
                                                             tuning_order)


    # Optimal parameters
    sparse_tuning_result_u[[j]] <- opt_parameters_result$sparse_tuning_selection_u # Optimal level of sparsity for u (CV)
    sparse_tuning_result_v[[j]] <- opt_parameters_result$sparse_tuning_selection_v # Optimal level of sparsity for v (CV)
    opt_alpha_result <- opt_parameters_result$GCV_score_smooth # Optimal Smoothness (GCV)

    opt_S_v[[j]] <- opt_alpha_result$opt_s.alpha
    smooth_tuning_result_v[[j]] <- opt_alpha_result$opt.alpha
    GCV_v[[j]] <- opt_alpha_result$GCV
    GCVdf_v[[j]] <- opt_alpha_result$GCVdf

    # Handling two-way smoothness
    if(two_way_smoothness != 0){
      opt_S_u[[j]] <- opt_alpha_result$opt_s.alpha_u
      smooth_tuning_result_u[[j]] <- opt_alpha_result$opt.alpha_u
      GCV_u[[j]] <- opt_alpha_result$GCV_u
      GCVdf_u[[j]] <- opt_alpha_result$GCVdf_u
    }else{
      opt_S_u[[j]] <- diag(n)
      smooth_tuning_result_u[[j]] <- 0
      GCV_u[[j]] <- Inf
      GCVdf_u[[j]] <- data.frame(0,Inf)
    }

    # Extracting v and u having the optimal parameters
    test_result <- power_algo(data = X_temp,
                              sparse_tuning_result_u = sparse_tuning_result_u[[j]],
                              sparse_tuning_result_v = sparse_tuning_result_v[[j]],
                              S_alpha_v = opt_S_v[[j]],
                              S_alpha_u = opt_S_u[[j]],
                              sparse_tuning_type = sparse_tuning_type,
                              type = "real")

    v <- test_result[[1]]
    u <- test_result[[2]]
    lsv <- cbind(lsv, v)
    lsu <- cbind(lsu, u)
    funcs[[j]] <- u%*%t(v)
  }


  # Splitting v for variables
  PCs <- list()
  for (i in 1:n_var) {
    lsv <- data.frame(lsv)
    rows_to_extract <- 1:as.integer(ncol[i])
    PCs[[i]] <- lsv[rows_to_extract,]
    lsv <- lsv[-rows_to_extract,]
  }



  return(list(Estimated = funcs, PC_functions = PCs, PC_Scores = lsu,
              opt_alpha_for_PC = smooth_tuning_result_v, opt_alpha_for_u = smooth_tuning_result_u,
              opt_gamma_for_PC = sparse_tuning_result_v, opt_gamma_for_u = sparse_tuning_result_u,
              GCV_v = GCV_v, GCVdf_v = GCVdf_v, GCV_u = GCV_u, GCVdf_u = GCVdf_u))
}

