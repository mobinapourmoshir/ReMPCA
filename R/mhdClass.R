#' @title  Define an object for the elements of the hybrid data.
#'
#' #' @description
#' The `hd` class denotes an element of hybrid data, which can be either functional data or regular data.
#' Data that is functional Objects are constructed using matrices, with columns representing grid points and rows indicating observations.
#' The Smoothing_parameter attributes are the distinguishing factor between the functional and regular data.
#' The smoothness of the data will not be implemented if the value is zero, or user can choose a fix number for smoothing parameter.
#' Alternatively, a vector of numerical values that will undergo generalized cross validation (GCV) to determine the most optimal one.
#'
#'
#' @param data A matrix which represents the data with rows indicating observations.
#' @param argval A vector of grid points for functional data, with a length equal to the number of columns in the data.
#'
#' @param Smoothing_parameter A fix number representing the smoothing parameter,
#' or a vector of numerical values that will undergo generalized cross validation (GCV) to determine the most optimal one.
#' Set it to 0 for no smoothing.
#' @param Sparsity_parameter A fix number representing the level of sparsity,
#' or a vector of numerical values that will undergo cross validation (CV) to determine the most optimal one.
#' For no sparsity, set it to 0 and tune it automatically by setting it to NULL.
#'
#' @description Check for validity of the data in the `hd` object
#'
#' @examples
#'
#' # Create some example matrices
#' mat1 <- matrix(c(1:30), nrow = 5)
#' mat2 <- matrix(c(1:25), nrow = 5)
#'
#' obj1 <- hdClass(mat1,
#'                 argval = NULL,
#'                 Smoothing_parameter = 0,
#'                 Sparsity_parameter = NULL)
#'
#' obj2 <- hdClass(mat2, attr = 0)
#'
#' hybrid.data <- list(obj1, obj2)
#'
#' # Print the object
#' print(obj1)
#'
#'
#' print(obj2)
#'
#' @export

####################### Define an S3 Class for Hybrid Data #######################
hdClass <- function(data,
                    argval = NULL,
                    Smoothing_parameter = 0,
                    Sparsity_parameter = 0) {

  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }

  # Validate 'Smoothing_parameter' to be a number, vector, 0, or NULL
  if (!is.null(Smoothing_parameter) &&
      !is.numeric(Smoothing_parameter) &&
      !identical(Smoothing_parameter, 0)) {
    stop("Smoothing_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  # Assign the attribute to the matrix
  if(is.null(Smoothing_parameter)){
    Smoothing_parameter <- 2^seq(-30,5, length.out = 10)}
  attr(data, "Smoothing_parameter") <- Smoothing_parameter


  # Validate 'Sparsity_parameter' to be a number, vector, 0, or NULL
  if (!is.null(Sparsity_parameter) &&
      !is.numeric(Sparsity_parameter) &&
      !identical(Sparsity_parameter, 0)) {
    stop("Sparsity_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if(any(Sparsity_parameter > ncol(data))){
    warning("An integer between 0 and ncol(data) must be used to represent the level of sparsity. Setting the 'Sparsity_parameter' to NULL!")
    Sparsity_parameter <- NULL
  }

  # Assign the attribute to the matrix
  if(is.null(Sparsity_parameter)){
    Sparsity_parameter <- seq(from = 0, to = ncol(data)-1, by = 1)}
  attr(data, "Sparsity_parameter") <- Sparsity_parameter


  # Set the class of the object
  class(data) <- "hd"

  return(data)
}
