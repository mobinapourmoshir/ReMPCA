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
#' @param two_way_smoothness A logical; if True, the function implements the two-way smoothness on both u and v and
#' if False (default) it only implement the smoothness on principal components (v).
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

ReMPCA <- function(mhd_obj,
                   argval = NULL,
                   centerfns = TRUE,
                   num_pcs = 1,
                   smooth_tuning = NULL,
                   smoothness_type = "Second_order",
                   sparse_tuning_type = "soft",
                   sparse_tuning = 0,
                   K_fold = 5,
                   two_way_smoothness = FALSE,
                   two_way_sparsity = FALSE) {

  fd <- mhd_obj$fd
  nfd <- mhd_obj$nfd
  n <- nrow(fd[[1]]) # Number of observations (same for both fd, nfd)
  fd_n_var <- length(fd) # Number of variables (functional data)

  n_cols_fd <- as.data.frame(sapply(fd, dim))[2,] # Number of columns of each variable in fd
  n_cols_nfd <- sum(as.data.frame(sapply(nfd, dim))[2,]) # Number of columns of each variable in nfd

  ####### Smoothing Parameter (for functional data only) ##########
  if (all(is.vector(smooth_tuning)) & all(smooth_tuning !=0)) { # For a given vector (fixed and pre-defined)
    if(length(smooth_tuning) != fd_n_var){
      warning("The length of 'smooth_tuning' does not match 'p'. Setting 'smooth_tuning' to NULL!")
      smooth_tuning <- NULL
    }
    smooth_tuning <- matrix(smooth_tuning, nrow = 1)

  } else if(!(is.null(smooth_tuning))){
    if((smooth_tuning == 0 ||
        all(smooth_tuning == 0))){
      smooth_tuning <- matrix(rep(0,fd_n_var), nrow = 1)

  }}else if(is.matrix(smooth_tuning)){

    if(ncol(smooth_tuning) > fd_n_var){
      warning("'smooth_tuning' matrix should have p columns for p variables. Considering the first p columns!")
      smooth_tuning <- smooth_tuning[,1:p]

    }

    smooth_tuning <- data.frame(matrix(smooth_tuning)) # For a given matrix
    colnames(smooth_tuning) <- paste0("var", seq_along(smooth_tuning))

  }else{
    smooth_tuning <- NULL
  }


  # smooth_tuning = NUll -> Assigning different combinations of alphas for p variables
  if (is.null(smooth_tuning)) {

    for (i in 1:fd_n_var) {
      smooth_tuning <- c(smooth_tuning, list(2^seq(-30,5, length.out = 10)))# 10 alphas for each variable
    }
    smooth_tuning <- expand.grid(smooth_tuning) # Matrix of all possible alphas for p variables
  }



  ####### level of sparsity (for both functional and non-functional data) #######
  # level of sparsity for u can be either 0 or any number between 1 through the length of u (Coefficients - # of observations)
  if (is.null(sparse_tuning) ||
      sparse_tuning > n) {

    sparse_tuning_u <- seq(from = 0, to = n-1, by = 1)

    } else{
      sparse_tuning_u <- sparse_tuning
      }

  # level of sparsity for v can be  either 0 or any number between 1 through the length of v (# of columns)
  if(two_way_sparsity == TRUE){

    if(is.null(sparse_tuning) ||
       sparse_tuning > sum(n_cols_fd) + n_cols_nfd){

      sparse_tuning_v <-  seq(from = 0, to = sum(n_cols_fd) + sum(n_cols_nfd) -1, by = 1)

    }else{
      sparse_tuning_v <- sparse_tuning}
  }



  ####### Pre-processing: Centralizing the data #######
  # Functional data
  X <- c()
  if (centerfns) {
    for (p in 1:fd_n_var) {
      c <-  apply(fd[[p]], 2, function(x) x - mean(x))
      X <- cbind(X,c) # Demeaned Side by side functional data
    }
  }else{X <- do.call(cbind, fd)}

  # non-functional data
  Y <- c()
  if (centerfns) {
    c <-  apply(nfd[[1]], 2, function(x) x - mean(x)) # Demeaned non-functional data
    Y <- cbind(Y,c)
    }else{Y <- do.call(cbind, nfd)}


  ####### Grid Points (input or assigning) - Smoothness for functional data only #######
  if(is.null(argval) == FALSE){
    if(length(argval) != fd_n_var ||
       sum(as.data.frame(sapply(argval, length))) != sum(n_cols_fd)){ # if argval is not defined appropriately
      warning("'argval' is not assigned appropriately!")
      argval <- NULL
    }
  }


  GridPoints_v <- GridPoints_u <- list()

  if (!is.null(argval)) {
    GridPoints_v <- argval
    for (i in 1:fd_n_var) {
      cycle_u <- seq(1:nrow(fd[[i]])) / nrow(fd[[i]])
      GridPoints_u[[i]] <- cycle_u
    }

  } else {
    for (i in 1:fd_n_var) {
      cycle_v <- seq(1:ncol(fd[[i]])) / ncol(fd[[i]])
      GridPoints_v[[i]] <- cycle_v
    }
    cycle_u <- seq(1:nrow(fd[[1]])) / nrow(fd[[1]])
    GridPoints_u <- cycle_u
  }


  ####### S_alpha for all alphas #######
  alphas <- smooth_tuning
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
    for (i in 1:fd_n_var) {
      alpha <- as.numeric(smooth_tuning[alpha_index,i])

      S_v[[i]] <- get.pen(td = GridPoints_v[[i]],
                          alpha = alpha,
                          type = smoothness_type)
    }

    S_u <- get.pen(td = GridPoints_u,
                        alpha = alpha,
                        type = smoothness_type)

    S_alpha_list_v[[index]] <- as.matrix(bdiag(S_v))
    S_alpha_list_u[[index]] <- as.matrix(bdiag(S_u))
    setTxtProgressBar(pb, index)
  }
  close(pb)



  ####### ReMPCA Implementation #######
  for (j in 1:num_pcs) {
    cat(sprintf("Computing the %s PC ...\n", ordinal(j)))
    if (j == 1) {
      X_temp = X
      Y_temp = Y

    } else{

      # Functional data
      SVD_result = svd(X_temp)
      v_original = SVD_result$v[,1]
      u_original = SVD_result$u[,1]
      sigma = SVD_result$d[1]
      X_temp = X_temp - sigma * u_original%*%t(v_original)


      # non-Functional data
      SVD_result_nfd = svd(Y_temp)
      v_original_nfd = SVD_result_nfd$v[,1]
      u_original_nfd = SVD_result_nfd$u[,1]
      sigma_nfd = SVD_result_nfd$d[1]
      Y_temp = Y_temp - sigma_nfd * u_original_nfd%*%t(v_original_nfd)
    }

    results <- Tuning_Power(X_temp =  X_temp,
                            Y_temp = Y_temp,
                            n_var = fd_n_var,
                            n_cols_fd = n_cols_fd,
                            n = n,
                            smooth_tuning = smooth_tuning,
                            sparse_tuning_u  = sparse_tuning_u,
                            sparse_tuning_v = sparse_tuning_v,
                            sparse_tuning_type = sparse_tuning_type,
                            K_fold = K_fold,
                            S_alpha_List_v = S_alpha_list_v,
                            S_alpha_list_u = S_alpha_list_u,
                            two_way_smoothness = two_way_smoothness ,
                            two_way_sparsity = two_way_sparsity,
                            j = j)

    result_fd <- results[[1]]  # Functional data results
    result_nfd <- results[[2]] # Non-Functional data results


  }

  # Splitting v for variables
  PCs <- list()
  for (i in 1:n_var) {
    lsv_fd <- data.frame(result_fd$lsv_fd)
    rows_to_extract <- 1:as.integer(n_cols[i])
    PCs[[i]] <- lsv_fd[rows_to_extract,]
    lsv_fd <- lsv_fd[-rows_to_extract,]
  }

  result_fd <- c(result_fd, PC_functions = PCs)

  return(list(result_fd, result_nfd))
}

