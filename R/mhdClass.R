#' @title  Define a Set of Multivariate Hybrid Data objects
#'
#' #' @description
#' The `mhd` class represents a set of multivariate hybrid data.
#' Functional data Objects are constructed using matrices, with columns representing grid points and rows indicating observations.
#'
#' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
#'
#' @param mvfd_obj List of matrices or arrays (Multivariate Grid Functional Data)
#' @param centerfns logical. If TRUE  the input data undergoes centralization or demeaning prior to being processed by the function.
#' @param num_pcs  The number of PCs. The default is one (The first principal component only). But the user is able to see higher PCs
#' @description Check for validity of the data in the `mvgfd` object
#'
#' @examples
#'
#' # Create some example matrices
#' obj1 <- hdClass(mat, attr = NULL)
#' obj2 <- hdClass(mat, attr = 0)
#'
#' # Print the object
#' print(obj1)
#' print(obj2)
#'
#' @export


####################### Define an S3 Class for Hybrid Data #######################
hdClass <- function(data, attr = 0) {
  if (!is.matrix(data)) {
    stop("Input 'data' must be a matrix.")
  }

  # Validate 'attr' to be a number, vector, 0, or NULL
  if (!is.null(attr) && !is.numeric(attr) && !identical(attr, 0)) {
    stop("Attribute must be a numeric value, numeric vector, 0, or NULL.")
  }

  # Assign the attribute to the matrix
  if(is.null(attr)){
    attr <- 2^seq(-30,5, length.out = 10)}
  attr(data, "custom_attr") <- attr

  # Set the class of the object
  class(data) <- "hd"

  return(data)
}
