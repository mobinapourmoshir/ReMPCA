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

