#' @title  Define an object for the elements of the hybrid data.
#'
#' #' @description
#' The `hd` class denotes an element of hybrid data, which can be either functional data or regular data.
#' Data that is functional Objects are constructed using matrices, with columns representing grid points and rows indicating observations.
#' The Smoothing_parameter attributes are the distinguishing factor between the functional and regular data.
#' The smoothness of the data will not be implemented if the value is zero, or user can choose a fix number for smoothing parameter.
#' Alternatively, a vector of numerical values that will undergo generalized cross validation (GCV) to determine the most optimal one.
#'
#' @param data A matrix which represents the data with rows indicating observations.
#' @param argval A vector of grid points for functional data, with a length equal to the number of columns in the data.
#' Grid poits from 0 to 1 will be assigned if NULL.
#' @param Smoothing_parameter A fix number representing the smoothing parameter,
#' or a vector of numerical values that will undergo generalized cross validation (GCV) to determine the most optimal one.
#' Set it to 0 for no smoothing. If NULL, it analyzes a sequence of 2^seq(-30, 5, length.out = 10) and attempts to tune it.
#' @param Sparsity_parameter_col A fix number representing the level of sparsity for column,
#' or a vector of numerical values that will undergo cross validation (CV) to determine the most optimal one.
#' For no sparsity, set it to 0 and tune it automatically by setting it to NULL.
#' @param Sparsity_parameter_row A fix number representing the level of sparsity for rows,
#' or a vector of numerical values that will undergo cross validation (CV) to determine the most optimal one.
#' For no sparsity, set it to 0 and tune it automatically by setting it to NULL.
#'
#'
#' @description Check for validity of the data in the `hd` object
#'
#' @examples
#'
#' # Create some example matrices
#' mat1 <- matrix(c(1:30), nrow = 5)
#' mat2 <- matrix(c(31:55), nrow = 5)
#'
#' obj1 <- hdClass(mat1,
#'                 argval = NULL,
#'                 Smoothing_parameter = 0,
#'                 Sparsity_parameter_col = NULL,
#'                 Sparsity_parameter_row = NULL)
#'
#' obj2 <- hdClass(mat2,
#'                argval = NULL,
#'                Smoothing_parameter = NULL,
#'                Sparsity_parameter_col = NULL,
#'                Sparsity_parameter_row = NULL)
#'
#' hybrid.data <- list(obj1, obj2)
#'
#' # Print the object
#' print(obj1)
#'
#' print(obj2)
#'
#' @export

####################### Define an S3 Class for Hybrid Data #######################
hdClass <- function(data,
                    argval = NULL,
                    Smoothing_parameter = 0,
                    two_way_smoothness = 0,
                    Sparsity_parameter_col = 0,
                    Sparsity_parameter_row = 0) {

  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }

  ####### Grid Points (input or assigning) - Smoothness for functional data only #######
  if (!is.null(argval)) {
    if (length(argval) != ncol(data)) {
      warning("There should be an equal number of grid points and columns. 'argval' is set to NULL!")
      argval <- NULL
    }
  }

  # Assigning GridPoints for u and v
  GridPoints_v <- GridPoints_u <- vector()

  if (!is.null(argval)) {
    GridPoints_v <- argval  # Use provided argval
  } else {
    GridPoints_v<- seq(from = 1/ncol(data), to = 1 , length.out =ncol(data)) # Default sequence
  }
  GridPoints_u <- seq(from = 1/nrow(data), to = 1 , length.out =nrow(data))


  # Validate 'Smoothing_parameter'
  if (!is.null(Smoothing_parameter) &&
      !is.numeric(Smoothing_parameter) &&
      !identical(Smoothing_parameter, 0)) {
    stop("Smoothing_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (is.null(Smoothing_parameter)) {
    Smoothing_parameter <- 2^seq(-30, 5, length.out = 10)
  }
  attr(data, "Smoothing_parameter") <- Smoothing_parameter

  # Validate 'Sparsity_parameter_col'
  if (!is.null(Sparsity_parameter_col) &&
      !is.numeric(Sparsity_parameter_col) &&
      !identical(Sparsity_parameter_col, 0)) {
    stop("Sparsity_parameter_col must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (any(Sparsity_parameter_col > ncol(data))) {
    warning("An integer between 0 and ncol(data) must be used to represent the level of sparsity for columns. Setting 'Sparsity_parameter_col' to NULL!")
    Sparsity_parameter_col <- NULL
  }

  if (is.null(Sparsity_parameter_col)) {
    Sparsity_parameter_col <- seq(from = 0, to = ncol(data) - 1, by = 1)
  }
  attr(data, "Sparsity_parameter_col") <- Sparsity_parameter_col

  # Validate 'Sparsity_parameter_row'
  if (!is.null(Sparsity_parameter_row) &&
      !is.numeric(Sparsity_parameter_row) &&
      !identical(Sparsity_parameter_row, 0)) {
    stop("Sparsity_parameter_row must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (any(Sparsity_parameter_row > nrow(data))) {
    warning("An integer between 0 and nrow(data) must be used to represent the level of sparsity for rows. Setting 'Sparsity_parameter_row' to NULL!")
    Sparsity_parameter_row <- NULL
  }

  if (is.null(Sparsity_parameter_row)) {
    Sparsity_parameter_row <- seq(from = 0, to = ncol(data) - 1, by = 1)
  }
  attr(data, "Sparsity_parameter_row") <- Sparsity_parameter_row

  # Assign GridPoints as attributes
  attr(data, "GridPoints_u") <- GridPoints_u
  attr(data, "GridPoints_v") <- GridPoints_v

  # Set the class of the object
  class(data) <- "hd"

  return(data)
}
