##################### Print method for regular data #####################
#' @export
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
  data <- x$matrix
  rows_to_show <- min(5, dim(x)[1])
  cols_to_show <- min(5, dim(x)[2])
  print(as.matrix(data)[1:rows_to_show, 1:cols_to_show])
}

#' Custom `$` operator for rdClass
#' Allows access to the underlying data matrix via `rd$matrix`
#'
#' @param x An object of class 'rdClass'
#' @param name The name of the element to extract
#' @export
`$.rdClass` <- function(x, name) {
  if (name == "matrix") {
    return(as.data.frame(unclass(x)))
  } else {
    stop(sprintf("Unknown field '%s'. Only 'matrix' is supported for rdClass."), call. = FALSE)
  }
}


#' Coerce an object of class 'fdClass' or 'imgClass' to class 'rdClass'
#'
#' This function strips the smoothing and grid-related structure from a functional object
#' and converts it to a regular data object of class \code{rdClass}.
#'
#' @param x An object of class \code{'fdClass'} or \code{'imgClass'}.
#' @param Sparsity_parameter Optional sparsity parameter to override the original.
#'   If \code{NULL}, the existing sparsity parameter (if any) is inherited.
#'
#' @return An object of class \code{rdClass} with smoothing and grid attributes removed,
#'         and sparsity parameter preserved or overridden.
#'
#' @examples
#' img_object <- imgClass(image = list(matrix(rnorm(100), nr = 50),
#'                                     matrix(rnorm(100), nr = 50)),
#'                        argval = NULL,
#'                        Smoothing_parameter = NULL,
#'                        Sparsity_parameter = 0)
#'
#' newrd <- as.rdClass(img_object, Sparsity_parameter = 1:10)
#' attr(newrd, "Sparsity_parameter")
#'
#' @export
as.rdClass <- function(x,
                       Sparsity_parameter = NULL) {
  # Validate input class
  if (!(inherits(x, "fdClass") ||
        inherits(x, "imgClass"))) {
    stop("Input must be of class 'fdClass', or 'imgClass'!")
  }

  # Convert to matrix
  data <- as.matrix(x)

  # Determine which sparsity parameter to use
  if (is.null(Sparsity_parameter)) {
    Sparsity_parameter <- attr(x, "Sparsity_parameter")
  }

  rdClass(data = data,
          Sparsity_parameter = Sparsity_parameter)
}

#' Multiply a `rdClass` Object by a Scalar
#'
#' @description Performs element-wise multiplication between a scalar and a `rdClass` object.
#'              Attributes like sparsity settings are preserved.
#'
#' @param e1 A scalar numeric value or a `rdClass` object.
#' @param e2 A `rdClass` object or a scalar numeric value.
#'
#' @return A `rdClass` object scaled by the numeric scalar.
#'
#' @examples
#' rd <- rdClass(matrix(1:12, nrow = 4))
#' rd_scaled <- rd * 0.5
#'
#' @export
`*.rdClass` <- function(e1, e2) {
  if (is.numeric(e1) && inherits(e2, "rdClass")) {
    out <- e1 * unclass(e2)
    attributes(out) <- attributes(e2)
    class(out) <- "rdClass"
    return(out)
  } else if (is.numeric(e2) && inherits(e1, "rdClass")) {
    out <- e2 * unclass(e1)
    attributes(out) <- attributes(e1)
    class(out) <- "rdClass"
    return(out)
  } else {
    stop("One operand must be numeric and the other an 'rdClass' object.")
  }
}

#' Plot Method for rdClass Objects
#'
#' Produces a scatter plot of regular data stored in an \code{rdClass} object.
#' Each column is plotted as a sequence of solid points.
#'
#' @param obj An object of class \code{rdClass}.
#'
#' @details
#' This method visualizes the regular (non-functional) data in the \code{rdClass} object.
#' It shows the values in each column as solid dots, which is useful for examining patterns across observations or variables.
#'
#' @return No return value. This function is called for its side effect (plot).
#'
#' @examples
#' rd_obj <- rdClass(matrix(rnorm(100), nrow = 10, ncol = 10))
#' plot(rd_obj)
#'
#' @export
plot.rdClass <- function(obj, ...) {
  matplot(obj, type = "p", pch = 16, main = "rd Class Plot", ...)
}

#' Indexing operator for rdClass
#'
#' Enables subsetting of an \code{rdClass} object by rows and columns.
#'
#' @param x An object of class \code{rdClass}.
#' @param i Row indices (observations). If \code{NULL}, all rows are included.
#' @param j Column indices. If \code{NULL}, all columns are included.
#'
#' @return A new \code{rdClass} object with subsetted data and inherited attributes.
#'
#' @export
`[.rdClass` <- function(x, i = NULL, j = NULL) {
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
  rdClass(data = data_sub,
          Sparsity_parameter = attr(x, "Sparsity_parameter"))
}
