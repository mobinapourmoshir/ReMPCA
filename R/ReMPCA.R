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
#' By default, it is null, and it looks at a matrix of all the possible alphas in 2^seq(-30,5, length.out = 15).
#' Set to 0 to have no smoothness.
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param sparse_tuning_fd A number that shows the level of sparsity for functional data.
#'  Set to 0 to have no sparsity (default). Tune it automatically by setting it to NULL.
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v, must be one of "Second_order" (default), "First_order" or "Indicator".
#'
#' @param two_way_smoothness A logical; if True, the function implements the two-way smoothness on both u and v and
#' if False (default) it only implement the smoothness on the functions (v).
#' @param two_way_sparsity A logical; if True, the function implements the two-way sparsity on both u and v and
#' if False (default) it only implement the sparsity on the coefficients (u).
#' @param sparse_tuning_nfd A number that shows the level of sparsity for non-functional data.
#'  Set to 0 to have no sparsity (default). Tune it automatically by setting it to NULL.
#'
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return PC scores, PC functions, ...
#' @export
#'


############################ Smooth and Sparse Multivariate PCA ############################

ReMPCA <- function(mhd_obj, argval = NULL, centerfns = TRUE, num_pcs = 1,
                      smooth_tuning = NULL, smoothness_type = "Second_order",
                      sparse_tuning_type = "soft", sparse_tuning_fd = 0,
                      two_way_smoothness = FALSE, two_way_sparsity = FALSE,
                      sparse_tuning_nfd = 0) {


  if(smooth_tuning = 0){ # No smoothness penalty
    mhd_obj <- c(mhd_obj$fd, mhd_obj$nfd) # All the f and nf matrices in one list
    n_var <- length(mhd_obj) # Number of variables
    n <- nrow(mhd_obj[[1]]) # Number of observations

    # Pre-processing: Centralizing the data
    Y <- c()
    if (centerfns) {
      for (p in 1:n_var) {
        c <-  apply(mhd_obj[[p]], 2, function(x) x - mean(x))
        Y <- cbind(Y,c)
      }
    }else{Y <- do.call(cbind, mhd_obj)} # Side_by_side data






  }else{ # Penalty is added

    # If penalty is added, just implement it for the functional part of data!
    mfd <- mhd_obj$fd
    mnfd <- mhd_obj$nfd
    n_var <- length(mfd) # Number of variables
    n <- nrow(mfd[[1]]) # Number of observations
    n_cols_fd <- as.data.frame(sapply(mfd, dim))[2,]
    n_cols_nfd <- as.data.frame(sapply(mnfd, dim))[2,]

    ####### Smoothing Parameter ##########

    # If smooth_tuning is a vector of p alphas (fixed and pre-defined)
    if (length(smooth_tuning) != n_var) {
      warning("The length of 'smooth_tuning' does not match 'p'. Setting 'smooth_tuning' to NULL.")
      smooth_tuning <- NULL
    } else {
      smooth_tuning <- data.frame(matrix(smooth_tuning, nrow = 1))
      colnames(smooth_tuning) <- paste0("var", seq_along(smooth_tuning))
    }


    if (is.null(smooth_tuning)) {
      for (i in 1:n_var) {
        smooth_tuning <- c(smooth_tuning, list(2^seq(-30,5, length.out = 15)))# 15 alphas for each variable
      }
    }
    smooth_tuning <- expand.grid(smooth_tuning) # Matrix of all possible alphas for p variables



    # level of sparsity
    # level of sparsity for u can be either 0 or any number between 1 through the length of u (Coefficients - # of observations)
    if (is.null(sparse_tuning_fd) | is.null(sparse_tuning_nfd) | sparse_tuning_fd > n | sparse_tuning_nfd > n) {
      sparse_tuning_u_fd <- sparse_tuning_u_nfd <- seq(0:floor(n-1))} else{
        sparse_tuning_u_fd <- sparse_tuning_fd
        sparse_tuning_u_nfd <- sparse_tuning_nfd}

    # level of sparsity for v can be  either 0 or any number between 1 through the length of v (# of columns)
    if(two_way_sparsity == TRUE){
      if(is.null(sparse_tuning_fd) | is.null(sparse_tuning_nfd) | sparse_tuning_fd > sum(n_cols_fd) | sparse_tuning_nfd > sum(n_cols_nfd)){
        sparse_tuning_v_fd <- seq(0:floor(sum(n_cols_fd)-1))
        sparse_tuning_v_nfd <- seq(0:floor(sum(n_cols_nfd)-1))
      }else{
        sparse_tuning_v_fd <- sparse_tuning_fd
        sparse_tuning_v_nfd <- sparse_tuning_nfd}
    }





    # Pre-processing: Centralizing the data
    X <- c()
    if (centerfns) {
      for (p in 1:n_var) {
        c <-  apply(mfd[[p]], 2, function(x) x - mean(x))
        X <- cbind(X,c)
      }
    }else{X <- do.call(cbind, mfd)}

    Y <- c()
    n_var_nfd <- length(mnfd)
    if (centerfns) {
      for (p in 1:n_var_nfd) {
        c <-  apply(mnfd[[p]], 2, function(x) x - mean(x))
        Y <- cbind(Y,c)
      }
    }else{Y <- do.call(cbind, mnfd)}




    # Grid Points (input or assigning)
    GridPoints_v <- GridPoints_u <- list()
    if (!is.null(argval)) {
      GridPoints <- argval
    } else {
      for (i in 1:n_var) {
        cycle_v <- seq(1:ncol(mfd[[i]])) / ncol(mfd[[i]])
        cycle_u <- seq(1:nrow(mfd[[i]])) / nrow(mfd[[i]])

        GridPoints_v[[i]] <- cycle_v
        GridPoints_u[[i]] <- cycle_u
      }
    }


    # S_alpha for all alphas
    alphas <- smooth_tuning
    S_alpha_list_u <- S_alpha_list_v <- list()
    index <- 0
    cat("Preprocessing ...\n")
    n_iter1 <- dim(smooth_tuning)[1]
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = n_iter1, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar

    for (alpha_index in 1:nrow(smooth_tuning)) {
      index <- index + 1
      S <- list()
      for (i in 1:n_var) {
        alpha <- as.numeric(smooth_tuning[alpha_index,i])
        S_v[[i]] <- get.pen(td = GridPoints_v[[i]], alpha = alpha, type = smoothness_type)
        S_u[[i]] <- get.pen(td = GridPoints_u[[i]], alpha = alpha, type = smoothness_type)
      }
      S_alpha_list_v[[index]] <- as.matrix(bdiag(S_v))
      S_alpha_list_u[[index]] <- as.matrix(bdiag(S_u))
      setTxtProgressBar(pb, index)
    }
    close(pb)


    smooth_tuning_result  <- sparse_tuning_result <- list()
    gcv <- opt_S  <- funcs <- GCVdf <- list()

    } # Penalty is added




  lsv <- lsu <- c() # List for storing v's  and u's
  variance <- vector() # % of variability explained by PC



  for (j in 1:num_pcs) {
    cat(sprintf("Computing the %s PC ...\n", ordinal(j)))
    if (j == 1) {
      X_temp = X # Penalty added

    } else{
      SVD_result = svd(X_temp)
      v_original = SVD_result$v[,1]
      u_original = SVD_result$u[,1]
      sigma = SVD_result$d[1]
      X_temp = X_temp - sigma * u_original%*%t(v_original)
    }


    # Tuning Parameters
    opt_parameters_result <- opt_alpha_result <- list()
    opt_parameters_result <- parameter_selection_conditional(data = X_temp, nvar = n_var, ncol = n_cols, smooth_tuning = smooth_tuning,
                                                             sparse_tuning = sparse_tuning, sparse_tuning_type = sparse_tuning_type, K_fold = 5, S_alpha_List = S_alpha_list)



    ################# Take care of sparse_tuning = 0 #################
    sparse_tuning_result[[j]] <- opt_parameters_result[[1]] # Optimal level of sparsity (CV)
    opt_alpha_result <- opt_parameters_result[[2]] # Optimal Smoothness (GCV)

    opt_S[[j]] <- opt_alpha_result$opt_s.alpha
    smooth_tuning_result[[j]] <- opt_alpha_result$opt.alpha
    gcv[[j]] <- opt_alpha_result$GCV
    GCVdf[[j]] <- opt_alpha_result$GCVdf



    # Extracting v and u having the optimal parameters
    test_result <- power_algo(data = X_temp, sparse_tuning_result = sparse_tuning_result[[j]] ,
                              sparse_tuning_type = sparse_tuning_type, S_alpha = opt_alpha_result$opt_s.alpha, type = "real")

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
    rows_to_extract <- 1:as.integer(n_cols[i])
    PCs[[i]] <- lsv[rows_to_extract,]
    lsv <- lsv[-rows_to_extract,]
  }


  return(list(Estimated = funcs, PC_functions = PCs, PC_Scores = lsu, opt_alpha_for_PC = smooth_tuning_result,
              opt_gamma_for_PC = sparse_tuning_result, GCV = gcv, GCV_df = GCVdf))
}

