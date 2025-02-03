################# Method for `hd` class #################
print.hd <- function(x, ...) {
  cat("Hybrid Data (hdClass) Object\n")
  cat("-----------------------------------\n")
  cat("Dimensions: ", dim(x)[1], "x", dim(x)[2], "\n")

  # Print the Smoothing parameter
  smoothing_param <- attr(x, "Smoothing_parameter")
  cat("Smoothing Parameter: ")
  if (length(smoothing_param) > 5) {
    cat(head(smoothing_param, 3), "...", tail(smoothing_param, 2), "\n")
  } else {
    cat(smoothing_param, "\n")
  }

  # Print the Sparsity parameter
  sparsity_param <- attr(x, "Sparsity_parameter")
  cat("Sparsity Parameter: ")
  if (length(sparsity_param) > 5) {
    cat(head(sparsity_param, 3), "...", tail(sparsity_param, 2), "\n")
  } else {
    cat(sparsity_param, "\n")
  }

  cat("-----------------------------------\n")
  cat("First few rows of the data:\n")
  print(head(x))

  invisible(x)  # Return the object invisibly
}
