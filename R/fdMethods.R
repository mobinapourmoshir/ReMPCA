##################### Print method for functional data #####################
#' @export
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
  data <- x$matrix
  rows_to_show <- min(5, dim(x)[1])
  cols_to_show <- min(5, dim(x)[2])
  print(as.matrix(data)[1:rows_to_show, 1:cols_to_show])

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


#' Coerce an object of class 'rdClass' or 'imgClass' to class 'fdClass'.
#'
#' @param x An object of class 'rdClass' or 'imgClass'.
#' @param Smoothing_parameter Optional smoothing parameter to assign. If \code{NULL}, defaults are used.
#' @param Sparsity_parameter Optional sparsity parameter to assign. If \code{NULL}, defaults are used.
#' @param argval Optional vector of grid points for columns. If \code{NULL}, defaults are used.
#'
#' @return An object of class 'fdClass' with user-specified or inherited regularization parameters.
#'
#' @examples
#' img_object <- imgClass(
#'   image = list(matrix(rnorm(100), nrow = 10), matrix(rnorm(100), nrow = 10)),
#'   argval = NULL,
#'   Smoothing_parameter = NULL,
#'   Sparsity_parameter = 0
#' )
#' newfd <- as.fdClass(img_object)
#' newfd
#'
#' @export
as.fdClass <- function(x,
                       Sparsity_parameter = NULL,
                       Smoothing_parameter = NULL,
                       argval = NULL) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "imgClass"))) {
    stop("Input must be of class 'rdClass', or 'imgClass'!")
  }

  # Convert to matrix
  data <- as.matrix(x)

  # Extract sparsity parameter
  if (is.null(Sparsity_parameter)) {
    Sparsity_parameter <- attr(x, "Sparsity_parameter")
  }else{
    Sparsity_parameter <- Sparsity_parameter}

  # Extract smoothness parameter
  if (is.null(Smoothing_parameter)) {
    Smoothing_parameter <- attr(x, "Smoothing_parameter")
  }else{
    Smoothing_parameter <- Smoothing_parameter}
  if (is.null(argval)) {
    argval <- attr(x, "GridPoints_v")
  }else{
    argval <- argval}

  # Construct fdClass with overrides
  fdClass(data = data,
          argval = argval,
          Smoothing_parameter = Smoothing_parameter,
          Sparsity_parameter = Sparsity_parameter)
}


#' Plot Method for fdClass Objects
#'
#' Generates a line plot of the functional data stored in an object of class \code{fdClass}.
#'
#' @param x An object of class \code{fdClass}.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @details
#' This function uses \code{\link{matplot}} to visualize each observation (row) as a separate curve.
#' It provides a quick overview of the functional data structure stored in the \code{fdClass} object.
#'
#' @return No return value. This function is called for its side effect (plot).
#'
#' @examples
#' fd_obj <- fdClass(matrix(sin(1:100 / 10), nrow = 10, ncol = 10))
#' plot(fd_obj)
#'
#' @export
plot.fdClass <- function(x, ...) {
  matplot(x, type = "l", main = "fd Class Plot", ...)
}

#' Multiply a `fdClass` Object by a Scalar
#'
#' @description Performs element-wise multiplication between a scalar and a `fdClass` object.
#'              All functional data attributes are preserved in the result.
#'
#' @param e1 A scalar numeric value or a `fdClass` object.
#' @param e2 A `fdClass` object or a scalar numeric value.
#'
#' @return A new `fdClass` object with elements scaled by the scalar value.
#'
#' @examples
#' fd <- fdClass(matrix(1:20, 10, 2), Smoothing_parameter = 0.1)
#' scaled_fd <- 3 * fd
#'
#' @export
`*.fdClass` <- function(e1, e2) {
  if (is.numeric(e1) && inherits(e2, "fdClass")) {
    out <- e1 * unclass(e2)
    attributes(out) <- attributes(e2)
    class(out) <- "fdClass"
    return(out)
  } else if (is.numeric(e2) && inherits(e1, "fdClass")) {
    out <- e2 * unclass(e1)
    attributes(out) <- attributes(e1)
    class(out) <- "fdClass"
    return(out)
  } else {
    stop("One operand must be numeric and the other an 'fdClass' object.")
  }
}

#' Indexing operator for fdClass
#'
#' Enables subsetting of an \code{fdClass} object by rows and columns.
#'
#' @param x An object of class \code{fdClass}.
#' @param i Row indices (observations). If \code{NULL}, all rows are included.
#' @param j Column indices (grid points). If \code{NULL}, all columns are included.
#'
#' @return A new \code{fdClass} object with subsetted data and inherited attributes.
#'
#' @export
`[.fdClass` <- function(x, i = NULL, j = NULL) {
  if (is.null(i) && is.null(j)) {
    return(x)
  }
  n <- nrow(x)
  m <- ncol(x)

  # Default to full selection
  if (is.null(i)) i <- seq_len(n)
  if (is.null(j)) j <- seq_len(m)

  # Bounds check
  if (any(i < 1 | i > n)) stop("Row index out of bounds.")
  if (any(j < 1 | j > m)) stop("Column index out of bounds.")

  # Subset data matrix
  data <- x$matrix
  data_sub <- as.matrix(data[i, j])

  # Construct and return new hdClass object
  fdClass(data = data_sub,
          Sparsity_parameter = attr(x, "Sparsity_parameter"),
          Smoothing_parameter = attr(x, "Smoothing_parameter"),
          argval = attr(x, "GridPoints_v"))
}
