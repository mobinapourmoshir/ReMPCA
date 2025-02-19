#' @title  Functional, Regular, and Hybrid Structures
#'
#' @description
#' The `fd` class represents an element of functional data.
#' - Data objects are constructed using matrices, where columns represent grid points (user-defined or NULL) and rows represent observations.
#' - Users can assign both smoothing and sparsity parameters.
#' - If the smoothing parameter is set to zero, no smoothing is applied. Otherwise, users can specify a fixed smoothing value or provide a vector of values, which will be optimized using generalized cross-validation (GCV).
#'
#' The `rd` class denotes an element of regular data.
#' - Data objects are constructed with observations in rows and variables in columns.
#' - The smoothing parameter is automatically set to zero (no smoothing).
#' - Users can assign sparsity parameters only!
#' - No grid points are involved in `rd` objects.
#'
#' The `hd` class represents hybrid data.
#' - It is a list containing any number of `fd` and/or `rd` objects in no specific order.
#' - Users can assign smoothing and sparsity parameters for the rows, as well as grid points for the rows.
#' - All `fd` and `rd` objects within the `hd` list must have the same number of observations (rows).

#'
#' @param data A matrix representing the data, with rows indicating observations and columns representing variables (for `rd`) or grid points (for `fd`).
#'
#' @param argval A vector of grid points.
#' - For functional data (`fd`), it assigns grid points to the columns, with a length equal to the number of columns in the data.
#' - For hybrid data (`hd`), it assigns grid points to the rows, with a length equal to the number of rows in the data.
#' - If `NULL`, grid points are automatically assigned from 0 to 1.
#'
#' @param Smoothing_parameter
#' - In `fd`: Smoothing parameter for columns. It can be:
#'   - A fixed number representing the smoothing parameter.
#'   - A vector of numerical values, which will undergo generalized cross-validation (GCV) to determine the optimal value.
#'   - Set to 0 for no smoothing.
#'   - If `NULL`, it analyzes a sequence of `2^seq(-30, 5, length.out = 10)` and attempts to tune it.
#' - In `hd`: Smoothing parameter for rows, can be defined in the same way as for the columns in `fd`.
#'
#' @param Sparsity_parameter
#' - A fixed number representing the level of sparsity for columns, or a vector of numerical values that will undergo cross-validation (CV) to determine the optimal value.
#' - For no sparsity, set it to 0.
#' - If `NULL`, the sparsity parameter will be tuned automatically.
#' - Can be defined for both `fd` and `rd` objects.
#'
#' @param row_sparsity_parameter
#' - A fixed number representing the level of sparsity for rows, or a vector of numerical values that will undergo cross-validation (CV) to determine the optimal value.
#' - For no sparsity, set it to 0.
#' - If `NULL`, the row sparsity parameter will be tuned automatically.
#' - This parameter is defined only for `hd` objects.
#'
#'
#'
#' @example
#' # Example for Functional Data (fd)
#' fd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example functional data matrix (10 rows, 10 columns)
#' fd_object <- fdClaa(data = fd_data,
#'                     argval = seq(0, 1, length.out = 10),  # Grid points for columns
#'                     Smoothing_parameter = 0.5,  # Custom smoothing parameter
#'                     Sparsity_parameter = 2)  # Custom sparsity parameter
#'
#' # Display the created fd object
#' print(fd_object)
#' print(attr(fd_object, "GridPoints_v"))  # Display grid points for columns
#' print(attr(fd_object, "Smoothing_parameter"))  # Display smoothing parameter
#'
#' # Example for Regular Data (rd)
#' rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # Example regular data matrix (10 rows, 10 columns)
#' rd_object <- rdClaa(data = rd_data,
#'                     Sparsity_parameter = 3)  # Custom sparsity parameter
#'
#' # Display the created rd object
#' print(rd_object)
#' print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter
#'
#' # Example for Hybrid Data (hd)
#' fd_object2 <- fdClaa(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another fd object
#' rd_object2 <- rdClaa(data = matrix(rnorm(100), nrow = 10, ncol = 10))  # Another rd object
#'
#' hd_list <- list(fd_object, rd_object)  # List of fd and rd objects
#' hd_object <- hdClass(hdlist = hd_list,
#'                      row_argval = seq(0, 1, length.out = 10),  # Grid points for rows
#'                      row_smoothing_parameter = 0.5,  # Custom smoothing parameter for rows
#'                      row_sparsity_parameter = 2)  # Custom sparsity parameter for rows
#'
#' # Display the created hd object
#' print(hd_object)
#' print(attr(hd_object, "GridPoints_u"))  # Display grid points for rows
#' print(attr(hd_object, "row_smoothing_parameter"))  # Display row smoothing parameter
#' print(attr(hd_object, "row_sparsity_parameter"))  # Display row sparsity parameter
#'
#' @export

####################### Define an S3 Class for Hybrid Data #######################
############ Functional Data ############
fdClaa <- function(fdmatrix,
                   argval = NULL,
                   Smoothing_parameter = 0,
                   Sparsity_parameter = 0){

  # Validation on the fdmatrix
  if (!is.matrix(fdmatrix)) {
    stop("Input 'data' must be a matrix.")
  }

  ####### Grid Points #######
  if (!is.null(argval)) {
    if (length(argval) != ncol(fdmatrix)) {
      warning("There should be an equal number of grid points and columns. 'argval' is set to NULL!")
      argval <- NULL
    }
  }

  # Assigning GridPoints for u and v
  GridPoints_v <- vector()

  if (!is.null(argval)) {
    GridPoints_v <- argval  # Use provided argval
  } else {
    GridPoints_v<- seq(from = 1/ncol(fdmatrix), to = 1 , length.out =ncol(fdmatrix)) # Default sequence
  }
  attr(fdmatrix, "GridPoints_v") <- GridPoints_v

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
  attr(fdmatrix, "Smoothing_parameter") <- Smoothing_parameter

  ####### Sparsity_parameter #######
  # Validate 'Sparsity_parameter'
  if (!is.null(Sparsity_parameter) &&
      !is.numeric(Sparsity_parameter) &&
      !identical(Sparsity_parameter, 0)) {
    stop("Sparsity_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (any(Sparsity_parameter > ncol(fdmatrix))) {
    warning("An integer between 0 and ncol(fdmatrix) must be used to represent the level of sparsity for columns. Setting 'Sparsity_parameter' to NULL!")
    Sparsity_parameter <- NULL
  }

  if (is.null(Sparsity_parameter)) {
    Sparsity_parameter <- seq(from = 0, to = ncol(fdmatrix) - 1, by = 1)
  }
  attr(fdmatrix, "Sparsity_parameter") <- Sparsity_parameter

  # Set the class of the object
  class(fdmatrix) <- "fd"
  return(fdmatrix)

}


############ Regular Data ############
rdClaa <- function(data,
                   Sparsity_parameter = 0){

  # Validation on the data
  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }


  ####### Smoothing_parameter #######
  attr(data, "Smoothing_parameter") <- 0

  ####### Sparsity_parameter #######
  # Validate 'Sparsity_parameter'
  if (!is.null(Sparsity_parameter) &&
      !is.numeric(Sparsity_parameter) &&
      !identical(Sparsity_parameter, 0)) {
    stop("Sparsity_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (any(Sparsity_parameter > ncol(data))) {
    warning("An integer between 0 and ncol(data) must be used to represent the level of sparsity for columns. Setting 'Sparsity_parameter' to NULL!")
    Sparsity_parameter <- NULL
  }

  if (is.null(Sparsity_parameter)) {
    Sparsity_parameter <- seq(from = 0, to = ncol(data) - 1, by = 1)
  }
  attr(data, "Sparsity_parameter") <- Sparsity_parameter

  # Set the class of the object
  class(data) <- "rd"

  return(data)

}


############ Hybrid Data ############
hdClass <- function(hdlist,
                    row_argval = NULL,
                    row_smoothing_parameter = 0,
                    row_sparsity_parameter = 0) {

  if(!(is.list(hdlist))){
    hdlist <- list(hdlist)
  }

  # Validate input: Ensure all elements are of class "fdClass" or "rdClass"
  if (!all(sapply(hdlist, function(obj) any(class(obj) %in% c("fd", "rd"))))) {
    stop("All elements in the list must be of class 'fd' or 'rd'.")
  }

  # Make sure that all matrices have the same number of observations
  nrows <- sapply(hdlist, function(obj) nrow(obj))
  if (!all(nrows == nrows[1])) {
    stop("Error: Not all matrices have the same number of rows!")
  }

  ####### Grid Points for rows #######
  if (!is.null(row_argval)) {
    if (length(row_argval) != nrow(data)) {
      warning("There should be an equal number of grid points and rows for hybrid data. 'row_argval' is set to NULL!")
      row_argval <- NULL
    }
  }

  # Assigning GridPoints for u
  GridPoints_u <- vector()

  if (!is.null(row_argval)) {
    GridPoints_u <- row_argval  # Use provided row_argval
  } else {
    GridPoints_u <- seq(from = 1/nrow(data), to = 1 , length.out =nrow(data))
  }
  attr(hdlist, "GridPoints_u") <- GridPoints_u

  ####### row_smoothing_parameter #######
  # Validate 'row_smoothing_parameter'
  if (!is.null(row_smoothing_parameter) &&
      !is.numeric(row_smoothing_parameter) &&
      !identical(row_smoothing_parameter, 0)) {
    stop("row_smoothing_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (is.null(row_smoothing_parameter)) {
    row_smoothing_parameter <- 2^seq(-30, 5, length.out = 10)
  }
  attr(hdlist, "row_smoothing_parameter") <- row_smoothing_parameter

  ####### row_sparsity_parameter #######
  # Validate 'row_sparsity_parameter'
  if (!is.null(row_sparsity_parameter) &&
      !is.numeric(row_sparsity_parameter) &&
      !identical(row_sparsity_parameter, 0)) {
    stop("row_sparsity_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (any(row_sparsity_parameter > nrow(data))) {
    warning("An integer between 0 and nrow(data) must be used to represent the level of sparsity for rows Setting 'row_sparsity_parameter' to NULL!")
    row_sparsity_parameter <- NULL
  }

  if (is.null(row_sparsity_parameter)) {
    row_sparsity_parameter <- seq(from = 0, to = nrow(data) - 1, by = 1)
  }
  attr(hdlist, "row_sparsity_parameter") <- row_sparsity_parameter


  # Set the class of the object
  class(hdlist) <- "hd"

  return(hdlist)
}
