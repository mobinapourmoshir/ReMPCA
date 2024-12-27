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
  combined_fd <- do.call(cbind, fd_matrices)

  # Combine non-functional data matrices side by side
  combined_nfd <- do.call(cbind, nfd_matrices)

  # Add attributes to the data
  attr(combined_fd, "label") <- "fd"
  attr(combined_nfd, "label") <- "nfd"

  # Create the object
  obj <- list(
    fd = combined_fd,
    nfd = combined_nfd
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
  cat("Attributes:", attributes(object$fd), "\n")
  cat("\nNon-Functional Data (nfd):\n")
  if (!is.null(object$nfd)) {
    print(object$nfd)
    cat("Attributes:", attributes(object$nfd), "\n")
  } else {
    cat("No non-functional data available.\n")
  }
}

# Example usage
# Create some example matrices
fd1 <- matrix(1:9, nrow = 3)
fd2 <- matrix(10:18, nrow = 3)
nfd1 <- matrix(19:27, nrow = 3)
nfd2 <- matrix(28:36, nrow = 3)

# Create an `hd` object
hd_obj <- hd(fd_matrices = list(fd1, fd2), nfd_matrices = list(nfd1, nfd2))

# Print the object
print(hd_obj)
