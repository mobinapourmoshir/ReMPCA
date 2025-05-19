##################### Print method for regular data #####################
print.rdClass <- function(x, ...) {

  sparsity_param <- attributes(x)["Sparsity_parameter"][[1]]
  cat("Regular Data (rdClass) Object\n")
  cat("-----------------------------------\n")
  cat("Dimensions: ", dim(x)[1], " x ", dim(x)[2], "\n", sep = "")

  # Print Sparsity Parameters
  cat("Sparsity Parameter: ")
  if (!is.null(sparsity_param)) {
    if (length(sparsity_param) > 5) {
      cat(head(sparsity_param, 3), "...", tail(sparsity_param, 2), "\n")
    } else {
      cat(sparsity_param, "\n")
    }
  } else {
    cat("NULL\n")
  }

  cat("-----------------------------------\n")
  cat("First few rows and columns of the data:\n")
  # Extract and print only the first 5 rows and 10 columns
  rows_to_show <- min(5, dim(x)[1])
  cols_to_show <- min(5, dim(x)[2])
  print(as.matrix(x[1:rows_to_show, 1:cols_to_show]))

}

#' Check if an object is of class 'rdClass'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'rdClass', FALSE otherwise.
#' @export
is.rdClass <- function(x) {
  inherits(x, "rdClass")
}

#' Coerce an object of class 'fdClass', 'rdClass', or 'hdClass' to class 'rdClass'
#'
#' @param x An object of class 'fdClass', 'rdClass', or 'hdClass'
#' @return An object of class 'rdClass' with only Sparsity_parameter(s) preserved.
#'         All smoothing-related attributes and grid points are removed.
#' @export
as.rdClass <- function(x) {
  # Validate input class
  if (!(inherits(x, "fdClass") ||
        inherits(x, "rdClass") ||
        inherits(x, "hdClass"))) {
    stop("Input must be of class 'fdClass', 'rdClass', or 'hdClass'")
  }

  # Convert to matrix
  data <- as.matrix(x)

  # Extract and preserve sparsity parameters
  sparsity_col <- attr(x, "Sparsity_parameter_col")
  sparsity <- attr(x, "Sparsity_parameter")

  # Strip all smoothing and grid-related attributes
  attr(data, "Smoothing_parameter") <- NULL
  attr(data, "Smoothing_parameter_col") <- NULL
  attr(data, "GridPoints_v") <- NULL
  attr(data, "GridPoints_u") <- NULL

  # Retain sparsity attribute(s)
  if (!is.null(sparsity)) {
    attr(data, "Sparsity_parameter") <- sparsity
  }
  if (!is.null(sparsity_col)) {
    attr(data, "Sparsity_parameter_col") <- sparsity_col
  }

  # Set class
  class(data) <- "rdClass"
  return(data)
}
