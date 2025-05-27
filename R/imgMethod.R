#' Print Method for Image Data (imgClass)
#'
#' Displays summary information about an imgClass object, including its dimensions,
#' smoothing/sparsity parameters (if any), and a preview of the image data.
#'
#' @param x An object of class \code{imgClass}.
#' @param ... Additional arguments (ignored).
#'
#' @return The object \code{x} is returned invisibly.
print.imgClass <- function(x, ...) {
  cat("Image Data (imgClass) Object\n")
  cat("-----------------------------------\n")

  cat("Dimensions: ", dim(x)[1], " x ", dim(x)[2], "\n", sep = "")

  # Optional attributes
  smoothing_param <- attr(x, "Smoothing_parameter")
  sparsity_param <- attr(x, "Sparsity_parameter")
  grid_points_v <- attr(x, "GridPoints_v")
  nrow_img <- attr(x, "nrow")

  # Show smoothing parameter
  cat("Smoothing Parameter: ")
  if (!is.null(smoothing_param)) {
    if (length(smoothing_param) > 5) {
      cat(head(smoothing_param, 3), "...", tail(smoothing_param, 2), "\n")
    } else {
      cat(smoothing_param, "\n")
    }
  } else {
    cat("None\n")
  }

  # Show sparsity parameter
  cat("Sparsity Parameter: ")
  if (!is.null(sparsity_param)) {
    if (length(sparsity_param) > 5) {
      cat(head(sparsity_param, 3), "...", tail(sparsity_param, 2), "\n")
    } else {
      cat(sparsity_param, "\n")
    }
  } else {
    cat("None\n")
  }

  # Show GridPoints_v
  cat("GridPoints_v: ")
  if (!is.null(grid_points_v)) {
    if (length(grid_points_v) > 5) {
      cat(head(grid_points_v, 3), "...", tail(grid_points_v, 2), "\n")
    } else {
      cat(grid_points_v, "\n")
    }
  } else {
    cat("None\n")
  }

  # Show image size if known
  if (!is.null(nrow_img)) {
    ncol_img <- ncol(x) / nrow_img
    cat("Each image size: ", nrow_img, " x ", ncol_img, "\n", sep = "")
  }

  cat("-----------------------------------\n")
  cat("First few rows and columns of the data:\n")
  rows_to_show <- min(5, nrow(x))
  cols_to_show <- min(5, ncol(x))
  print(x[1:rows_to_show, 1:cols_to_show])

  invisible(x)
}


#' Custom `$` operator for imgClass
#' Allows access to the underlying data matrix via `img$matrix`
#'
#' @param x An object of class 'imgClass'
#' @param name The name of the element to extract
#' @export
`$.imgClass` <- function(x, name) {
  if (name == "matrix") {
    return(as.data.frame(unclass(x)))
  } else {
    stop(sprintf("Unknown field '%s'. Only 'matrix' is supported for imgClass."), call. = FALSE)
  }
}


#' Coerce an Object to imgClass
#'
#' Converts an object of class \code{rdClass} or \code{fdClass} to an \code{imgClass},
#' while preserving or overriding associated attributes such as smoothing, sparsity,
#' and grid points.
#'
#' @param x An object of class \code{rdClass} or \code{fdClass}.
#' @param Sparsity_parameter Optional. A numeric vector of non-negative integers representing
#'        sparsity levels to apply. If \code{NULL}, the parameter is inherited from \code{x}.
#' @param Smoothing_parameter Optional. A numeric value or vector representing smoothing parameters.
#'        If \code{NULL}, the parameter is inherited from \code{x}.
#' @param argval Optional. A numeric vector of grid points (argvals) for functional representation.
#'        If \code{NULL}, the grid is inherited from \code{x}.
#'
#' @return An object of class \code{imgClass}, which also inherits from \code{fdClass} or \code{rdClass},
#' depending on the smoothing parameter.
#'
#' @details
#' This coercion is helpful when an object originally treated as regular or functional data
#' (via \code{rdClass} or \code{fdClass}) should instead be interpreted and processed as image data.
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' fd_obj <- fdClass(mat, Smoothing_parameter = 0.1)
#' img_obj <- as.imgClass(fd_obj)
#' print(class(img_obj))              # "imgClass" "fdClass"
#' print(attr(img_obj, "Smoothing_parameter"))
#'
#' @export
as.imgClass <- function(x,
                        Sparsity_parameter = NULL,
                        Smoothing_parameter = NULL,
                        argval = NULL) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "fdClass"))) {
    stop("Input must be of class 'rdClass', or 'fdClass'!")
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

  # Extract data matrix


  # Construct fdClass with overrides
  imgClass(image  = data,
           argval = argval,
           Smoothing_parameter = Smoothing_parameter,
           Sparsity_parameter = Sparsity_parameter)
}

