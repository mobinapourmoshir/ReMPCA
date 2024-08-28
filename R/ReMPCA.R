#' ReMPCA Smooth and Sparse Multivariate Functional Principal Component Analysis
#'
#' @param mvfd_obj A list of data matrices, where each one is considered of as a variable and observations are stored in the rows, and grid points are in the columns. It is also possible for the timeline or grid points to include the column name.
#' @param argval A list of grid points corresponding to each variable, where the length of each component matches the number of columns in the related data matrix.
#' @param centerfns A logical
#' @param num_pcs Logical: if True, it demeans the data before calculating the principal components.
#' @param smooth_tuning A vector with p elements that each represent a fixed smoothing parameter alpha for all p variables, OR A matrix with different combinations of alphas for all variables, OR A list of two vectors, one for each variable, that each represents possible alphas. By default, it is null, and it looks at a matrix of all the possible alphas in 2^seq(-30,5, length.out = 15).
#' @param sparse_tuning_type A character string specifying the sparse calculation method. Must be one of "soft" (default), "hard", or "SCAD".
#' @param sparse_tuning A number that shows the level of sparsity. Set to 0 to have no sparsity (default).
#' @param smoothness_type A character string specifying the method used in smoothing u and/or v, must be one of "Second_order" (default), "First_order" or "Indicator".
#'
#' @importFrom utils  txtProgressBar setTxtProgressBar
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return PC scores, PC functions, ...
#' @export
#'


############################ Smooth and Sparse Multivariate PCA ############################
ReMPCA <- function(mvfd_obj, argval = NULL, centerfns = TRUE, num_pcs = 1,
                       smooth_tuning = NULL, sparse_tuning_type = "soft",
                   sparse_tuning = 0, smoothness_type = "Second_order") {

  n_var <- length(mvfd_obj) # Number of variables
  n <- nrow(mvfd_obj[[1]]) # Number of observations
  n_cols <- as.vector(as.data.frame(sapply(mvfd_obj, dim))[2,])


  ####### Smoothing Parameter ##########
  if (is.null(smooth_tuning)) {
    for (i in 1:n_var) {
      smooth_tuning <- c(smooth_tuning, list(2^seq(-30,5, length.out = 15)))
    }
  }
  smooth_tuning <- expand.grid(smooth_tuning)

  # sparse_tuning is the level of sparsity
  # It can be either 0 or any number between 1 through the length of u (Coefficients)
  if (is.null(sparse_tuning)) {
    sparse_tuning <- seq(0:floor(n-1))
  }


  # Pre-processing: Centralizing the data
  X <- c()
  if (centerfns) {
    for (p in 1:n_var) {
      c <-  apply(mvfd_obj[[p]], 2, function(x) x - mean(x))
      X <- cbind(X,c)
    }
  }else{X <- do.call(cbind, mvfd_obj)}


  # Grid Points (input or assigning)
  GridPoints <- list()
  if (!is.null(argval)) {
    GridPoints <- argval
  } else {
    for (i in 1:n_var) {
      cycle <- seq(1:ncol(mvfd_obj[[i]])) / ncol(mvfd_obj[[i]])

      GridPoints[[i]] <- cycle
    }
  }


  # S_alpha for all alphas
  alphas <- smooth_tuning
  S_alpha_list <- list()
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
      S[[i]] <- get.pen(td = GridPoints[[i]], alpha = alpha)
    }
    S_alpha_list[[index]] <- as.matrix(bdiag(S))
    setTxtProgressBar(pb, index)
  }
  close(pb)



  lsv <- lsu <- c() # List for storing v's  and u's
  variance <- vector() # % of variability explained by PC
  smooth_tuning_result  <- sparse_tuning_result <- list()
  gcv <- opt_S  <- funcs <- GCVdf <- list()


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
    opt_parameters_result <- parameter_selection_conditional(data = X_temp, nvar = n_var, ncol = n_cols, smooth_tuning = smooth_tuning,
                                                             sparse_tuning = sparse_tuning, sparse_tuning_type = sparse_tuning_type, K_fold = 5, S_alpha_List = S_alpha_list)


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

