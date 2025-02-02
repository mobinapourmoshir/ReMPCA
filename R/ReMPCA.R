#' ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis
#'
#' @param mhd_obj Two lists of data matrices exist: one for functional data (fd_matrices) and another for non-functional data (nfd_matrices).
#' Each matrix is regarded as a variable, with observations organized in the rows and grid points in the columns for functional data.
#'  The timeline or grid points may also incorporate the column name.
#' @param argval A list of grid points for functional data corresponding to each variable, where the length of each component matches the number of columns in the related data matrix.
#' @param centerfns A logical; if True, it demeans the data before calculating the principal components.
#' @param num_pcs An integer. The number of principal components.
#' @param smooth_tuning A vector with p elements that each represent a fixed smoothing parameter alpha for all p variables,
#' OR A matrix with different combinations of alphas for all variables, it should have p columns for p variables!
#' By default, it is null, and it looks at a matrix of all the possible alphas in 2^seq(-30,5, length.out = 10).
#' Set to 0 to have no smoothness.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param sparse_tuning A number that shows the level of sparsity.
#'  Set to 0 to have no sparsity (default). Tune it automatically by setting it to NULL.
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v, must be one of "Second_order" (default), "First_order" or "Indicator".
#' @param K_fold An integer. It's used in cross validation approach for tuning the level of sparsity.
#'
#' @param two_way_smoothness A fix number that represents the smoothness parameter (alpha) for u. Or a vector of different alphas to be tuned.
#' Set to 0 (default) to have no smoothness.
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
                   argval = NULL,
                   centerfns = TRUE,
                   num_pcs = 1,
                   smoothness_type = "Second_order",
                   sparse_tuning_type = "soft",
                   sparse_tuning = 0,
                   K_fold = 5,
                   two_way_smoothness = 0,
                   two_way_sparsity = FALSE) {

  # Validate input: Ensure all elements are of class "hdClass"
  if (!all(sapply(object_list, function(obj) inherits(obj, "hd")))) {
    stop("All elements in the list must be of class 'hdClass'.")
  }


  # Combine matrices side by side
  hd <- do.call(cbind, object_list)
  n <- nrow(hd)

  n_var <- length(object_list) # Number of variables (# of matrices in object_list)
  ncol <- as.data.frame(sapply(object_list, dim))[2,]

  ####### Smoothing Parameter ##########
  # Generate all combinations alphas (one row per combination)
  smooth_tuning <- expand.grid(lapply(object_list, function(obj) attr(obj, "custom_attr")))


  ####### level of sparsity (for both functional and non-functional data) #######
  # level of sparsity for u can be either 0 or any number between 1 through the length of u (Coefficients - # of observations)
  # level of sparsity for v can be  either 0 or any number between 1 through the length of v (# of columns)

  if (is.null(sparse_tuning) ||
      sparse_tuning > n) {
    sparse_tuning_u <- seq(from = 0, to = n-1, by = 1)
    } else{
      sparse_tuning_u <- sparse_tuning
    }
  if(is.null(sparse_tuning) ||
     sparse_tuning > ncol(hd)){

    sparse_tuning_v <-  seq(from = 0, to = ncol(hd)-1, by = 1)

  }else{
    sparse_tuning_v <- sparse_tuning}


  ####### Pre-processing: Centralizing the data #######
  X <- data.frame()
  if (centerfns) {
    X <- apply(hd, 2, function(x) x - mean(x))
  }else{
    X <- hd
  }


  ####### Grid Points (input or assigning) - Smoothness for functional data only #######
  if(is.null(argval) == FALSE){
    if(length(argval) != n_var ||
       sum(as.data.frame(sapply(argval, length))) != sum(n_cols_fd)){ # if argval is not defined appropriately
      warning("'argval' is not assigned appropriately!")
      argval <- NULL
    }
  }

  GridPoints_v <- GridPoints_u <- list()

  if (!is.null(argval)) { # For the given argval for v
    GridPoints_v <- argval
  } else {
    for (i in 1:n_var) {
      cycle_v <- seq(1:ncol(object_list[[i]])) / ncol(object_list[[i]])
      GridPoints_v[[i]] <- cycle_v
    }
  }
  GridPoints_u <- seq(1:nrow(object_list[[1]])) / nrow(object_list[[1]])

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

  for (alpha_index in 1:length(two_way_smoothness)) {
    alpha <- as.numeric(two_way_smoothness[alpha_index])
    S_alpha_list_u[[alpha_index]] <- get.pen(td = GridPoints_u,
                                   alpha = two_way_smoothness,
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
                                                             sparse_tuning_u = sparse_tuning_u,
                                                             sparse_tuning_v = sparse_tuning_v,
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

