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

#' Check if an object is of class 'fdClass'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'fdClass', FALSE otherwise.
#' @export
is.fdClass <- function(x) {
  inherits(x, "fdClass")
}


#' Coerce an object of class 'rdClass', 'fdClass', or 'hdClass' to class 'fdClass'
#'
#' @param x An object of class 'rdClass', 'fdClass', or 'hdClass'
#' @return An object of class 'fdClass' with \code{Sparsity_parameter} preserved if present,
#'         and \code{Smoothing_parameter} and \code{argval} reset to default.
#' @export
as.fdClass <- function(x) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "fdClass") ||
        inherits(x, "hdClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', or 'hdClass'")
  }

  # Convert input to matrix (if already matrix-like)
  data <- as.matrix(x)

  # Extract sparsity parameter if available
  sparsity <- attr(x, "Sparsity_parameter")

  # Construct new fd object
  fd(data = data,
     argval = NULL,
     Smoothing_parameter = NULL,
     Sparsity_parameter = sparsity)
}
