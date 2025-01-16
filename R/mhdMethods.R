################# Method for `hd` class #################

# Print
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
