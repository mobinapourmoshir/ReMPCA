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
rdClass <- function(data,
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
  attr(data, "Smoothing_parameter") <- 0
  # Set the class of the object
  class(data) <- "rdClass"
  return(data)
}
