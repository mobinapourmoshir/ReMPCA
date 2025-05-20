#' ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis
#'
#' @param object_list An hdClass object.
#' @param centerhds A logical; if True, it demeans the data before calculating the principal components.
#' @param num_pcs An integer. The number of principal components.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v,
#' must be one of "Second_order" (default), "First_order" or "Indicator".
#' @param nfolds_u  An integer. It's used in cross validation approach for tuning the level of sparsity of rows.
#' @param nfolds_v An integer vector of length \code{p} (the number of variables), where
#' the \code{i}-th element specifies the number of cross-validation folds to use for the
#' \code{i}-th variable. If not provided (\code{NULL}), a value of 5 will be assigned to
#' all variables by default.
#' @param tuning_order A character string representing the tuning order.
#' If set to 'Sparsity', sparsity parameters are tuned first, followed by smoothness.
#'  If set to 'Smoothness', the order is reversed.
#' @param cv.pick A character string specifying the rule used to select the optimal
#' tuning parameter during cross-validation for sparsity. If set to `'min'`, the tuning
#' parameter corresponding to the minimum cross-validation error is chosen. If set to `'1se'`
#' (the default), the 1-standard-error rule is applied, selecting the most regularized
#' model whose error is within one standard error of the minimum.
#' @param parallel Logical; if \code{TRUE}, parallel computation is used to fit models
#' across different sparsity parameter values. Users must register a parallel backend
#' beforehand using packages such as \pkg{doParallel}, \pkg{doMC}, or similar.
#' @param tuning_iter Integer specifying the number of iterations to perform during the tuning process for conditional smoothing and sparsity parameters.
#' @param sparse_tuning_u Optional. Specifies the sparsity level(s) for rows:
#' \itemize{
#'   \item A single non-negative integer for fixed sparsity.
#'   \item A numeric vector of candidate values, to be selected via cross-validation (CV).
#'   \item Set to \code{0} for no sparsity.
#'   \item If \code{NULL}, the function defaults to the \code{Sparsity_parameter} attribute from the \code{hdClass} object.
#' }
#' @param sparse_tuning_v Optional. A list of length \code{p} (number of variables), where the \code{i}-th element is either a numeric value or a vector specifying candidate sparsity levels for the \code{i}-th variable.
#' If set to \code{NULL}, the function will default to using the \code{Sparsity_parameter_col} attribute from the input object of class \code{hdClass}.
#' @param smooth_tuning_u Optional. Specifies smoothing parameter for rows:
#' \itemize{
#'   \item A single number for fixed smoothness.
#'   \item A numeric vector of candidate values, to be selected via generalized cross-validation (GCV).
#'   \item Set to \code{0} for no smoothness
#'   \item If \code{NULL}, the function defaults to the \code{Smoothing_parameter} attribute from the \code{hdClass} object.
#' }
#' @param smooth_tuning_v Optional. A list of length \code{p} (number of variables), where the \code{i}-th element is either a numeric value or a vector specifying candidate smoothing parameters for the \code{i}-th variable.
#' If set to \code{NULL}, the function will default to using the \code{Smoothing_parameter_col} attribute from the input object of class \code{hdClass}.
#'
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return PC scores, PC functions, ...
#' @export
#'


################### Smooth and Sparse Multivariate PCA ###################
ReMPCA <- function(hd,
                   centerhds = TRUE,
                   num_pcs = 1,
                   smoothness_type = "Second_order",
                   sparse_tuning_type = "soft",
                   nfolds_u = 5,
                   nfolds_v = NULL,
                   thresh = 1e-10,
                   maxit = 100,
                   tuning_iter = 1,
                   parallel = FALSE,
                   tuning_order = "Sparsity",
                   cv.pick = "1se",
                   sparse_tuning_u = NULL,
                   sparse_tuning_v = NULL,
                   smooth_tuning_u = NULL,
                   smooth_tuning_v = NULL) {

  # Check if hd is a hd object
  if (!inherits(hd, "hdClass")) {
    stop("hd must be of class 'hdClass'!")
  }

  # Combine matrices side by side
  n <- nrow(hd)
  n_var <- attr(object_list, "n_var") # Number of variables (# of matrices in object_list)
  ncol <- data.frame(attr(object_list, "ncol")) # Number of columns of each matrix
  # nfolds_v with no default
  if (is.null(nfolds_v)) {
    nfolds_v <- rep(5, attr(hd, "n_var"))
  }

  ####### Smoothing Parameter ##########
  # Smoothing parameters for column
  # Generate all combinations alphas (one row per combination): A matrix
  if(!is.null(smooth_tuning_v)){
    smooth_tuning_col <- expand.grid(smooth_tuning_v)
    }else{
      smooth_tuning_col <- expand.grid(attr(object_list, "Smoothing_parameter_col"))
    }

  # Smoothing parameters for row: A vector
  if(!is.null(smooth_tuning_u)){
    smooth_tuning_row <- smooth_tuning_u
  }else{
    smooth_tuning_row <- attr(object_list, "Smoothing_parameter")
  }

  ####### level of sparsity #######
  if(!is.null(sparse_tuning_u)){
    sparsity_row_list <- sparse_tuning_u
  }else{
    sparsity_row_list <- attr(object_list, "Sparsity_parameter")
  }
  if(!is.null(sparse_tuning_v)){
    sparsity_col_list <- sparse_tuning_v
  }else{
    sparsity_col_list <- attr(object_list, "Sparsity_parameter_col")
  }

  ####### Pre-processing: Centralizing the data #######
  X <- data.frame()
  if (centerhds) {
    X <- apply(hd, 2, function(x) x - mean(x))
  }else{
    X <- hd
  }

  ####### Grid Points #######
  GridPoints_u <- attr(object_list, "GridPoints_u")
  GridPoints_v <- attr(object_list, "GridPoints_v")

  ####### S_alpha for all alphas #######
  S_alpha_list_u <- S_alpha_list_v <- Omegas_v <- list()
  index <- 0
  cat("Preprocessing ...\n")
  n_iter1 <- nrow(smooth_tuning_col)     # The number of alphas
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter1,# Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  # S_alpha for v
  for (alpha_index in 1:nrow(smooth_tuning_col)) {
    index <- index + 1
    S_u <- S_v <- Omega_v <- Omegas_u <- list()
    for (i in 1:n_var) {
      alpha <- as.numeric(smooth_tuning_col[alpha_index,i])
      if(is.null(GridPoints_v[[i]])){
        S_v[[i]] <- diag(ncol[,i])
        Omega_v[[i]] <- diag(ncol[,i])
      }else{
        get.pen.result <- get.pen(td = GridPoints_v[[i]],
                                  alpha = alpha,
                                  type = smoothness_type)
        S_v[[i]] <- get.pen.result$S.alpha
        Omega_v[[i]] <- get.pen.result$ Omega
      }
    }
    S_alpha_list_v[[index]] <- S_v #as.matrix(bdiag(S_v))
    Omegas_v[[index]] <- Omega_v
    setTxtProgressBar(pb, index)
  }

  # S_alpha for u
  if(!is.null(smooth_tuning_row)){
    for (alpha_index in 1:length(smooth_tuning_row)) {
      alpha <- as.numeric(smooth_tuning_row[alpha_index])
      get.pen.result <- get.pen(td = GridPoints_u,
                                alpha = alpha,
                                type = smoothness_type)

      S_alpha_list_u[[alpha_index]] <- get.pen.result$S.alpha
      Omegas_u[[alpha_index]] <- get.pen.result$Omega
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
                                                             smooth_tuning_v = smooth_tuning_col,
                                                             smooth_tuning_u = smooth_tuning_row,
                                                             sparse_tuning_u = sparsity_row_list,
                                                             sparse_tuning_v = sparsity_col_list,
                                                             sparse_tuning_type,
                                                             nfolds_u,
                                                             nfolds_v,
                                                             S_alpha_list_v ,
                                                             S_alpha_list_u,
                                                             Omegas_u = Omegas_u,
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

