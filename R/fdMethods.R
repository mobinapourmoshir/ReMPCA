##################### Print method for functional data #####################
print.fdClass <- function(x, ...) {
  cat("Functional Data (fdClass) Object\n")
  cat("-----------------------------------\n")
  cat("Dimensions: ", dim(x)[1], " x ", dim(x)[2], "\n", sep = "")

  # Get attributes
  smoothing_param <- attr(x, "Smoothing_parameter")
  sparsity_param <- attr(x, "Sparsity_parameter")
  grid_points_v <- attr(x, "GridPoints_v")

  # Print Smoothing Parameter
  cat("Smoothing Parameter: ")
  if (!is.null(smoothing_param)) {
    if (length(smoothing_param) > 5) {
      cat(head(smoothing_param, 3), "...", tail(smoothing_param, 2), "\n")
    } else {
      cat(smoothing_param, "\n")
    }
  } else {
    cat("NULL\n")
  }

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

  # Print GridPoints_v
  cat("GridPoints_v: ")
  if (!is.null(grid_points_v)) {
    if (length(grid_points_v) > 5) {
      cat(head(grid_points_v, 3), "...", tail(grid_points_v, 2), "\n")
    } else {
      cat(grid_points_v, "\n")
    }
  } else {
    cat("NULL\n")
  }

  cat("-----------------------------------\n")
  cat("First few rows and columns of the data:\n")
  rows_to_show <- min(5, dim(x)[1])
  cols_to_show <- min(5, dim(x)[2])
  print(x[1:rows_to_show, 1:cols_to_show])

  invisible(x)  # Return the object invisibly
}

#' Custom `$` operator for fdClass
#' Allows access to the underlying data matrix via `fd$matrix`
#'
#' @param x An object of class 'fdClass'
#' @param name The name of the element to extract
#' @export
`$.fdClass` <- function(x, name) {
  if (name == "matrix") {
    return(as.data.frame(unclass(x)))
  } else {
    stop(sprintf("Unknown field '%s'. Only 'matrix' is supported for fdClass."), call. = FALSE)
  }
}


#' Coerce an object of class 'rdClass', 'fdClass', or 'hdClass' to class 'fdClass'
#'
#' @param x An object of class 'rdClass', 'fdClass', or 'hdClass'.
#' @param Smoothing_parameter Optional smoothing parameter to assign. If \code{NULL}, defaults are used.
#' @param argval Optional vector of grid points for columns. If \code{NULL}, defaults are used.
#'
#' @return An object of class 'fdClass' with user-specified or inherited regularization parameters.
#' @export

as.fdClass <- function(x,
                       Smoothing_parameter = NULL,
                       argval = NULL) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "fdClass") ||
        inherits(x, "hdClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', or 'hdClass'")
  }

  # Convert to matrix
  data <- as.matrix(x)

  # Extract sparsity parameter
  sparsity <- attr(x, "Sparsity_parameter")

  # Construct fdClass with overrides
  fdClass(data = data,
          argval = argval,
          Smoothing_parameter = Smoothing_parameter,
          Sparsity_parameter = sparsity)
}
