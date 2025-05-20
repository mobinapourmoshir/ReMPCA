#' @title Hybrid Data Class
#'
#' @description Constructs an object of class `hdClass`, representing hybrid data
#' composed of a list of `fdClass` and/or `rdClass` objects.
#'
#' @details
#' The `hdClass` class represents hybrid data.
#' - It is a list containing any number of `fdClass` and/or `rdClass` objects in no specific order.
#' - Users can assign smoothing and sparsity parameters for the rows, as well as grid points for the rows.
#' - All `fdClass` and `rdClass` objects within the `hdClass` list must have the same number of observations (rows).
#'
#' @param hdlist A list of `rdClass` and/ or `hdClass` objects.
#'
#' @param argval A vector of grid points. For hybrid data (`hdClass`), it assigns grid points to
#' the rows, with a length equal to the number of rows in the data.
#' - If `NULL`, grid points are automatically assigned from 0 to 1.
#'
#' @param Smoothing_parameter : Smoothing parameter for rows It can be:
#'   - A fixed number representing the smoothing parameter.
#'   - A vector of numerical values, which will undergo generalized cross-validation (GCV) to determine the optimal value.
#'   - Set to 0 for no smoothing.
#'   - If `NULL`, it analyzes a sequence of `2^seq(-30, 5, length.out = 10)` and attempts to tune it.
#'
#' @param Sparsity_parameter
#' - A fixed number representing the level of sparsity for rows, or a vector of
#' numerical values that will undergo cross-validation (CV) to determine the optimal value.
#' - For no sparsity, set it to 0.
#' - If `NULL`, the row sparsity parameter will be tuned automatically.
#' - This parameter is defined only for `hdClass` objects.
#'
#' @example
#' Example for Hybrid Data (hd)
#' fd_object2 <- fdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another fd object
#' rd_object2 <- rdClass(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another rd object
#'
#' hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
#' hd_object <- hdClass(hdlist = hd_list,
#'                      argval = seq(0, 1, length.out = 10),  # Grid points for rows
#'                      Smoothing_parameter = 0.5,  # Custom smoothing parameter for rows
#'                      Sparsity_parameter = 2)  # Custom sparsity parameter for rows
#'
#' # Display the created hd object
#' print(hd_object)
#' print(attr(hd_object, "GridPoints_u"))  # Display grid points for rows
#' print(attr(hd_object, "Smoothing_parameter"))  # Display row smoothing parameter
#' print(attr(hd_object, "Sparsity_parameter"))  # Display row sparsity parameter
#'
#'
#' is.hdClass(hd_object)
#'
#' @export

####################### Define an S3 Class for Hybrid Data #######################
hdClass <- function(hdlist,
                    argval = NULL,
                    Smoothing_parameter = 0,
                    Sparsity_parameter = 0) {

  if(!(is.list(hdlist))){
    hdlist <- list(hdlist)
  }

  # Validate input: Ensure all elements are of class "fdClass" or "rdClass"
  if (!all(sapply(hdlist, function(obj) any(class(obj) %in% c("rdClass", "fdClass"))))) {
    stop("All elements in the list must be of class 'fdClass' or 'rdClass'.")
  }

  # Make sure that all matrices have the same number of observations
  nrows <- sapply(hdlist, function(obj) nrow(obj))
  if (!all(nrows == nrows[1])) {
    stop("Error: Not all matrices have the same number of rows!")
  }

  hd <- do.call(cbind, hdlist)

  ####### Grid Points for rows #######
  if (!is.null(argval)) {
    if (length(argval) != nrow(hd)) {
      warning("There should be an equal number of grid points and rows for hybrid data. 'argval' is set to NULL!")
      argval <- NULL
    }
  }

  # Assigning GridPoints for u
  GridPoints_u <- vector()

  if (!is.null(argval)) {
    GridPoints_u <- argval  # Use provided argval
  } else {
    GridPoints_u <- seq(from = 1/nrow(hd), to = 1 , length.out =nrow(hd))
  }
  attr(hd, "GridPoints_u") <- GridPoints_u

  ####### Smoothing_parameter #######
  # Validate 'Smoothing_parameter'
  if (!is.null(Smoothing_parameter) &&
      !is.numeric(Smoothing_parameter) &&
      !identical(Smoothing_parameter, 0)) {
    stop("Smoothing_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (is.null(Smoothing_parameter)) {
    Smoothing_parameter <- 2^seq(-30, 5, length.out = 10)
  }
  attr(hd, "Smoothing_parameter") <- Smoothing_parameter

  ####### Sparsity_parameter #######
  # Validate 'Sparsity_parameter'
  if (!is.null(Sparsity_parameter)) {
    # Check it's numeric and non-negative integers
    if (!is.numeric(Sparsity_parameter) ||
        any(Sparsity_parameter < 0) ||
        any(Sparsity_parameter != floor(Sparsity_parameter))) {
      stop("Sparsity_parameter must be a vector of non-negative integers.")
    }

    # Check that all values are within range
    if (any(Sparsity_parameter > nrow(hd) - 1)) {
      stop("All elements of Sparsity_parameter must be between 0 and nrow(hd) - 1.")
    }
  } else {
    # If NULL, generate sequence
    if (nrow(hd) <= 15) {
      Sparsity_parameter <- 0:(nrow(hd) - 1)
    } else {
      extra_vals <- unique(c(0:3, 2^(0:floor(log2(nrow(hd) - 1))), nrow(hd) - 1))
      Sparsity_parameter <- sort(unique(extra_vals[extra_vals <= (nrow(hd) - 1)]))
    }
  }

  attr(hd, "Sparsity_parameter") <- Sparsity_parameter
  attr(hd, "n_var") <- length(hdlist) # Number of variables (# of matrices in object_list)
  attr(hd, "ncol") <- as.data.frame(sapply(hdlist, dim))[2,] # Number of columns of each matrix

  ####### Smoothing parameters on columns #######
  Smoothing_parameter_col <- lapply(hdlist, function(obj) attr(obj, "Smoothing_parameter"))
  attr(hd, "Smoothing_parameter_col") <- Smoothing_parameter_col
  GridPoints_v <- lapply(hdlist, function(obj) attr(obj, "GridPoints_v"))
  attr(hd, "GridPoints_v") <- GridPoints_v

  ####### Sparsity parameters on columns #######
  Sparsity_parameter_col <- lapply(hdlist, function(obj) attr(obj, "Sparsity_parameter"))
  attr(hd, "Sparsity_parameter_col") <- Sparsity_parameter_col

  # Set the class of the object
  class(hd) <- "hdClass"
  return(hd)
}
