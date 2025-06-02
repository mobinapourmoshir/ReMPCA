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
  data <- x$matrix
  rows_to_show <- min(5, nrow(x))
  cols_to_show <- min(5, ncol(x))
  print(as.matrix(data)[1:rows_to_show, 1:cols_to_show])

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
  data <- list(as.matrix(x))

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
  imgClass(image  = data,
           argval = argval,
           Smoothing_parameter = Smoothing_parameter,
           Sparsity_parameter = Sparsity_parameter)
}

#' Plot Method for imgClass Objects
#'
#' Visualizes a list of image matrices stored in an \code{imgClass} object.
#' Each image is shown one at a time. The user is prompted to press Enter or
#' click on the graphics window to advance to the next image.
#'
#' @param x An object of class \code{imgClass}, created using \code{imgClass()}.
#' @param col Color palette to use for image rendering (default: grayscale).
#'
#' @return No return value. Called for its side effect (plotting).
#'
#' @examples
#' img_list <- list(matrix(rnorm(100), 10, 10),
#'                  matrix(runif(100), 10, 10))
#' img_obj <- imgClass(img_list)
#' plot(img_obj)
#'
#' @export
plot.imgClass <- function(x) {
  if (!inherits(x, "imgClass")) {
    stop("Input must be of class 'imgClass'.")
  }

  if (!is.list(x) && !is.matrix(x)) {
    stop("imgClass must contain a list of matrices or a single image matrix.")
  }

  # If it's a list of vectorized images, convert back
  if (!is.null(attr(x, "nrow"))) {
    n_imgs <- nrow(x)
    nrow_img <- attr(x, "nrow")
    ncol_img <- ncol(x) / nrow_img
    for (i in 1:n_imgs) {
      img <- matrix(x[i, ], nrow = nrow_img)
      image(img, main = paste("Image", i))
      if(n_imgs>=2 && i<= n_imgs -1){
        readline(prompt = "Press [Enter] or click to continue...")
      }
    }
  } else {
    stop("Unknown image structure. Missing 'nrow' attribute!")
  }

  invisible()
}

#' Multiply a `imgClass` Object by a Scalar
#'
#' @description Multiplies each matrix (image) in an `imgClass` object by a scalar.
#'              All attributes and class structure are retained.
#'
#' @param e1 A scalar numeric value or an `imgClass` object.
#' @param e2 An `imgClass` object or a scalar numeric value.
#'
#' @return A new `imgClass` object with each image scaled by the scalar value.
#'
#' @examples
#' img <- imgClass(image = list(matrix(1:9, 3, 3)))
#' img_scaled <- 2 * img
#'
#' @export
`*.imgClass` <- function(e1, e2) {
  if (is.numeric(e1) && inherits(e2, "imgClass")) {
    out_data <- lapply(unclass(e2), function(mat) e1 * mat)
    attributes(out_data) <- attributes(e2)
    class(out_data) <- class(e2)
    return(out_data)
  } else if (is.numeric(e2) && inherits(e1, "imgClass")) {
    out_data <- lapply(unclass(e1), function(mat) e2 * mat)
    attributes(out_data) <- attributes(e1)
    class(out_data) <- class(e1)
    return(out_data)
  } else {
    stop("One operand must be numeric and the other an 'imgClass' object.")
  }
}

#' Indexing operator for imgClass
#'
#' Enables subsetting of an \code{imgClass} object by rows.
#'
#' @param x An object of class \code{imgClass}.
#' @param i Row indices (Images). If \code{NULL}, all images are included.
#'
#' @return A new \code{imgClass} object with subsetted data and inherited attributes.
#'
#' @export
`[.imgClass` <- function(x, i = NULL) {
  if (is.null(i)) {
    return(x)
  }
  n <- nrow(x)

  # Default to full selection
  if (is.null(i)) i <- seq_len(n)

  # Bounds check
  if (any(i < 1 | i > n)) stop("Row index out of bounds.")

  # Subset data matrix
  data <- x$matrix
  data_sub <- list(as.matrix(data[i, ]))

  # Construct and return new hdClass object
  imgClass(image = data_sub,
          Sparsity_parameter = attr(x, "Sparsity_parameter"),
          Smoothing_parameter = attr(x, "Smoothing_parameter"),
          argval = attr(x, "argval"))
}
