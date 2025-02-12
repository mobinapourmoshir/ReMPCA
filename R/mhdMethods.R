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

  # Print the Sparsity parameter for column
  Sparsity_parameter_col <- attr(x, "Sparsity_parameter_col")
  cat("Columns Sparsity Parameter: ")
  if (length(Sparsity_parameter_col) > 5) {
    cat(head(Sparsity_parameter_col, 3), "...", tail(Sparsity_parameter_col, 2), "\n")
  } else {
    cat(Sparsity_parameter_col, "\n")
  }
  # Print the Sparsity parameter for row
  Sparsity_parameter_row <- attr(x, "Sparsity_parameter_row")
  cat("Rows Sparsity Parameter: ")
  if (length(Sparsity_parameter_row) > 5) {
    cat(head(Sparsity_parameter_row, 3), "...", tail(Sparsity_parameter_row, 2), "\n")
  } else {
    cat(Sparsity_parameter_row, "\n")
  }

  cat("-----------------------------------\n")
  cat("First few rows of the data:\n")
  print(head(x))

  invisible(x)  # Return the object invisibly
}
