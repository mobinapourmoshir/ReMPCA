#' @title Image Data Class
#'
#' @description
#' The `imgClass` class is designed to handle images:
#' - Supports either a single image (a matrix) or a list of multiple images (matrices).
#' - When a list is provided, all images are vectorized (row-wise) and combined into one matrix.
#' - The resulting object behaves like a `fdClass` if `Smoothing_parameter` is not zero, otherwise as a `rdClass`.
#' - Grid points (optional) are assigned to the columns (i.e., pixel locations).
#'
#' @param image A matrix representing a single image or a list of matrices representing multiple images.
#'              Each image must be a matrix of the same dimension.
#' @param Smoothing_parameter
#'   - A numeric scalar or vector controlling the level of smoothing (as in `fdClass`).
#'   - Set to 0 for no smoothing (produces an `rdClass` object).
#'   - If NULL, smoothing is tuned from a default sequence.
#' @param Sparsity_parameter
#'   - A fixed number or vector controlling the level of sparsity on pixel columns.
#'   - If NULL, it is tuned automatically based on image size.
#' @param argval Optional numeric vector for grid points over the pixels (columns).
#'
#' @return An object of class `fdClass` or `rdClass`, depending on the `Smoothing_parameter`.
#'
#' @examples
#' img1 <- matrix(rnorm(64), 8, 8)
#' img2 <- matrix(rnorm(64), 8, 8)
#'
#' # Single image
#' fd_img <- imgClass(img1, Smoothing_parameter = 0.5, Sparsity_parameter = 2)
#'
#' # Multiple images
#' rd_img <- imgClass(list(img1, img2), Smoothing_parameter = 0, Sparsity_parameter = NULL)
#'
#' @export
imgClass <- function(image,
                     Smoothing_parameter = 0,
                     Sparsity_parameter = 0,
                     argval = NULL) {

  if (is.list(image)) {
    # Validate: All elements should be matrices of the same dimension
    dims <- lapply(image, dim)
    if (!all(sapply(dims, function(d) all(d == dims[[1]])))) {
      stop("All images in the list must have the same dimensions.")
    }

    # Vectorize each image in row-major order
    data <- t(sapply(image, function(img) c(t(img))))
    nrow_img <- nrow(image[[1]])
  } else if (is.matrix(image)) {
    data <- image
    nrow_img <- NULL
  } else {
    stop("Input must be a matrix or a list of matrices.")
  }

  # Determine whether smoothing is needed
  if (all(Smoothing_parameter != 0)) {
    # Smoothing is requested → treat as fdClass and ignore argval
    x <- fdClass(data = data,
                 argval = NULL,
                 Smoothing_parameter = Smoothing_parameter,
                 Sparsity_parameter = Sparsity_parameter)
  } else {
    # No smoothing → treat as regular data class
    x <- rdClass(data = data,
                 Sparsity_parameter = Sparsity_parameter)
  }

  # Tag as image-specific object
  class(x) <- c("imgClass", class(x))

  if (!is.null(nrow_img)) {
    attr(x, "nrow") <- nrow_img
  }

  return(x)
}
