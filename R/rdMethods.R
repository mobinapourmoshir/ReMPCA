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


#' Coerce an object of class 'fdClass', 'rdClass', or 'hdClass' to class 'rdClass'
#'
#' @param x An object of class 'fdClass', 'rdClass', or 'hdClass'.
#' @param Sparsity_parameter Optional sparsity parameter to override the original.
#'
#' @return An object of class 'rdClass' with smoothing and grid attributes removed,
#'         and sparsity parameter preserved or overridden.
#' @export

as.rdClass <- function(x,
                       Sparsity_parameter = NULL) {
  # Validate input class
  if (!(inherits(x, "fdClass") ||
        inherits(x, "rdClass") ||
        inherits(x, "hdClass"))) {
    stop("Input must be of class 'fdClass', 'rdClass', or 'hdClass'")
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
