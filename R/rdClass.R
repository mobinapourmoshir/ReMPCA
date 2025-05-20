#' @title  Regular Data Class
#'
#' @description
#' The `rdClass` class denotes an element of regular data.
#' - Data objects are constructed with observations in rows and variables in columns.
#' - The smoothing parameter is automatically set to zero (no smoothing).
#' - Users can assign sparsity parameters only!
#' - No grid points are involved in `rdClass` objects.
#'
#' @param data A matrix representing the data, with rows indicating observations and columns representing variables.
#' @param Sparsity_parameter
#' - A fixed number representing the level of sparsity for columns, or a vector of numerical values that will undergo cross-validation (CV) to determine the optimal value.
#' - For no sparsity, set it to 0.
#' - If `NULL`, the sparsity parameter will be tuned automatically.
#'
#' @example
#' # Example for Regular Data (rd)
#' rd_data <- matrix(rnorm(100), nrow = 10, ncol = 10)  # 10 rows, 10 columns
#' rd_object <- rdClass(data = rd_data,
#'                      Sparsity_parameter = 3)  # Custom sparsity parameter
#'
#' # Display the created rd object
#' print(rd_object)
#' print(attr(rd_object, "Sparsity_parameter"))  # Display sparsity parameter
#'
#' is.rd(rd_object)
#' is.fd(rd_object)
#'
#' convert2fd <- as.fdClass(rd_object)
#' convert2rd <- as.rdClass(fd_object)
#'
#' @export

####################### Define an S3 Class for Regular Data #######################
####################### Define an S3 Class for Regular Data #######################
rdClass <- function(data,
                    Sparsity_parameter = 0) {

  # Validation on the data
  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }

  ####### Smoothing_parameter (always 0 for rdClass) #######
  attr(data, "Smoothing_parameter") <- 0

  ####### Sparsity_parameter #######
  if (!is.null(Sparsity_parameter)) {
    # Must be numeric, non-negative integers
    if (!is.numeric(Sparsity_parameter) ||
        any(Sparsity_parameter < 0) ||
        any(Sparsity_parameter != floor(Sparsity_parameter))) {
      stop("Sparsity_parameter must be a vector of non-negative integers.")
    }

    # Ensure values are in valid range
    if (any(Sparsity_parameter > ncol(data) - 1)) {
      stop("All elements of Sparsity_parameter must be between 0 and ncol(data) - 1.")
    }
  } else {
    # Generate default sequence
    if (ncol(data) <= 15) {
      Sparsity_parameter <- 0:(ncol(data) - 1)
    } else {
      extra_vals <- unique(c(0:3, 2^(0:floor(log2(ncol(data) - 1))), ncol(data) - 1))
      Sparsity_parameter <- sort(unique(extra_vals[extra_vals <= (ncol(data) - 1)]))
    }
  }

  attr(data, "Sparsity_parameter") <- Sparsity_parameter

  # Set the class
  class(data) <- "rdClass"
  return(data)
}
