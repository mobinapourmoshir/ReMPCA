#' @title  Define a Set of Multivariate Hybrid Data objects
#'
#' #' @description
#' The `mhd` class represents a set of multivariate hybrid data.
#' Functional data Objects are constructed using matrices, with columns representing grid points and rows indicating observations.
#'
#' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
#'
#' @param mvfd_obj List of matrices or arrays (Multivariate Grid Functional Data)
#' @param centerfns logical. If TRUE  the input data undergoes centralization or demeaning prior to being processed by the function.
#' @param num_pcs  The number of PCs. The default is one (The first principal component only). But the user is able to see higher PCs
#' @description Check for validity of the data in the `mvgfd` object
#'


####################### Define an S3 Class for Hybrid Data #######################

# Constructor for `hd` objects (Hybrid Data)
hd <- function(fd_matrices = list(), nfd_matrices = list()) {
  # Validate inputs
  if (!all(sapply(fd_matrices, is.matrix))) {
    stop("All elements in `fd_matrices` must be matrices.")
  }
  if (!all(sapply(nfd_matrices, is.matrix))) {
    stop("All elements in `nfd_matrices` must be matrices.")
  }

  # Combine functional data matrices side by side
  #combined_fd <- do.call(cbind, fd_matrices)
  # Combine non-functional data matrices side by side
  #combined_nfd <- do.call(cbind, nfd_matrices)

  # Add attributes to the data
  attr(fd_matrices, "label") <- "fd"
  attr(nfd_matrices, "label") <- "nfd"

  # Create the object
  obj <- list(
    fd = fd_matrices,
    nfd = nfd_matrices
  )

  # Assign a class attribute
  class(obj) <- "hd"

  return(obj)
}

# Print method for `hd` class
print.hd <- function(object) {
  cat("Hybrid Data Object:\n")
  cat("Functional Data (fd):\n")
  print(object$fd)
  cat("Attributes:\n")
  print(attributes(object$fd))
  cat("\nNon-Functional Data (nfd):\n")
  if (!is.null(object$nfd)) {
    print(object$nfd)
    cat("Attributes:\n")
    print(attributes(object$nfd))
  } else {
    cat("No non-functional data available.\n")
  }
}

# Example
# Create some example matrices
fd1 <- matrix(1:9, nrow = 3)
fd2 <- matrix(10:18, nrow = 3)
nfd1 <- matrix(19:27, nrow = 3)
nfd2 <- matrix(28:36, nrow = 3)

# Create an `hd` object
hd_obj <- hd(fd_matrices = list(fd1, fd2), nfd_matrices = list(nfd1, nfd2))

# Print the object
print(hd_obj)

