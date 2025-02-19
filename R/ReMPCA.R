#' ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis
#'
#' @param mhd_obj A list of of all hdClass objects. Each might represent either functional data (including grid points in the columns)
#' or regular data (characterized by a smoothness attribute of zero). The order of putting functional and regular data does not matter.
#' @param centerfns A logical; if True, it demeans the data before calculating the principal components.
#' @param num_pcs An integer. The number of principal components.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param sparse_tuning A number that shows the level of sparsity.
#'  Set to 0 to have no sparsity (default). Tune it automatically by setting it to NULL.
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v, must be one of "Second_order" (default), "First_order" or "Indicator".
#' @param K_fold An integer. It's used in cross validation approach for tuning the level of sparsity.
#'
#' @param two_way_smoothness A fix number that represents the smoothness parameter (alpha) for u. Or a vector of different alphas to be tuned.
#' Set to 0 (default) to have no smoothness. If NULL, it analyzes a sequence of 2^seq(-30, 5, length.out = 10) and attempts to tune it.
#' @param two_way_sparsity A logical; if True, the function implements the two-way sparsity on both u and v and
#' if False (default) it only implement the sparsity on the principal components (v).
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
                   sparse_tuning = 0,
                   two_way_smoothness = 0,
                   K_fold = 5,
                   two_way_sparsity = FALSE) {

  if(!(is.list(object_list))){
    object_list <- list(object_list)
  }

  # Validate input: Ensure all elements are of class "hdClass"
  if (!all(sapply(object_list, function(obj) inherits(obj, "hd")))) {
    stop("All elements in the list must be of class 'hdClass'.")
  }

  # Make sure that all matrices have the same number of observations
  nrows <- sapply(object_list, function(obj) nrow(obj))
  if (!all(nrows == nrows[1])) {
    stop("Error: Not all matrices have the same number of rows!")
  }

  # Combine matrices side by side
  hd <- do.call(cbind, object_list)
  n <- nrow(hd)
  n_var <- length(object_list) # Number of variables (# of matrices in object_list)
  ncol <- as.data.frame(sapply(object_list, dim))[2,] # Number of columns of each matrix

  ####### Smoothing Parameter ##########
  # Generate all combinations alphas (one row per combination)
  smooth_tuning <- expand.grid(lapply(object_list, function(obj) attr(obj, "Smoothing_parameter")))

  ####### level of sparsity (for both functional and non-functional data) #######
  sparsity_row_list <- lapply(object_list, function(obj) attr(obj, "Sparsity_parameter_row"))
  sparsity_col_list <- lapply(object_list, function(obj) attr(obj, "Sparsity_parameter_col"))

  ####### Pre-processing: Centralizing the data #######
  X <- data.frame()
  if (centerfns) {
    X <- apply(hd, 2, function(x) x - mean(x))
  }else{
    X <- hd
  }

  ####### Grid Points #######
  GridPoints_u <- lapply(object_list, function(obj) attr(obj, "GridPoints_u"))[[1]] # A vector
  GridPoints_v <- lapply(object_list, function(obj) attr(obj, "GridPoints_v")) # A list

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

      S_v[[i]] <- get.pen(td = GridPoints_v[[i]],
                          alpha = alpha,
                          type = smoothness_type)
    }
    S_alpha_list_v[[index]] <- as.matrix(bdiag(S_v))
    setTxtProgressBar(pb, index)
  }

  # S_alpha for u
  if(is.null(two_way_smoothness)){
    two_way_smoothness <- 2^seq(-30,5, length.out = 10)
  }
  for (alpha_index in 1:length(two_way_smoothness)) {
    alpha <- as.numeric(two_way_smoothness[alpha_index])
    S_alpha_list_u[[alpha_index]] <- get.pen(td = GridPoints_u,
                                             alpha = alpha,
                                             type = smoothness_type)
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
                                                             S_alpha_list_u ,
                                                             two_way_smoothness,
                                                             two_way_sparsity)


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

