#' @title  Functional Data Class
#'
#' @description
#' The `fd` class represents an element of functional data.
#' - Data objects are constructed using matrices, where columns represent grid points (user-defined or NULL) and rows represent observations.
#' - Users can assign both smoothing and sparsity parameters.
#' - If the smoothing parameter is set to zero, no smoothing is applied. Otherwise, users can specify a fixed smoothing value or provide a vector of values, which will be optimized using generalized cross-validation (GCV).
#'
#' @param data A matrix representing the data, with rows indicating observations and columns representing grid points.
#'
#' @param argval A vector of grid points, it assigns grid points to the columns, with a length equal to the number of columns in the data.
#' If `NULL`, grid points are automatically assigned from 0 to 1.
#'
#' @param Smoothing_parameter : Smoothing parameter for columns. It can be:
#'   - A fixed number representing the smoothing parameter.
#'   - A vector of numerical values, which will undergo generalized cross-validation (GCV) to determine the optimal value.
#'   - Set to 0 for no smoothing.
#'   - If `NULL`, it analyzes a sequence of `2^seq(-30, 5, length.out = 10)` and attempts to tune it.
#'
#' @param Sparsity_parameter
#' - A fixed number representing the level of sparsity for columns, or a vector of numerical values that will undergo cross-validation (CV) to determine the optimal value.
#' - For no sparsity, set it to 0.
#' - If `NULL`, the sparsity parameter will be tuned automatically.
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

####################### Define an S3 Class for Functional Data #######################
fdClass <- function(data,
                    argval = NULL,
                    Smoothing_parameter = 0,
                    Sparsity_parameter = 0){

  # Validation on the data
  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }

  #obj <- list(data = data)
  obj = data

  ####### Grid Points #######
  if (!is.null(argval)) {
    if (length(argval) != ncol(data)) {
      warning("There should be an equal number of grid points and columns. 'argval' is set to NULL!")
      argval <- NULL
    }
  }

  # Assigning GridPoints for u and v
  GridPoints_v <- vector()

  if (!is.null(argval)) {
    GridPoints_v <- argval  # Use provided argval
  } else {
    GridPoints_v <- seq(from = 1/ncol(data), to = 1 , length.out =ncol(data)) # Default sequence
  }
  attr(obj, "GridPoints_v") <- c(GridPoints_v)

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
  attr(obj, "Smoothing_parameter") <- as.vector(Smoothing_parameter)

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
  attr(obj, "Sparsity_parameter") <- as.vector(Sparsity_parameter)

  # Set the class of the object
  class(obj) <- "fdClass"
  return(obj)

}
