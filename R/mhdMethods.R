################# Method for `hd` class #################
print.hd <- function(x, ...) {
  cat("Hybrid Data (hdClass) Object\n")
  cat("-----------------------------------\n")
  cat("Dimensions: ", dim(x)[1], " x ", dim(x)[2], "\n", sep = "")

  # Print the Smoothing parameter
  smoothing_param <- attr(x, "Smoothing_parameter")
  cat("Smoothing Parameter: ")
  if (!is.null(smoothing_param)) {
    if (length(smoothing_param) > 5) {
      cat(head(smoothing_param, 3), "...", tail(smoothing_param, 2), "\n")
    } else {
      cat(smoothing_param, "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print the Sparsity parameter for columns
  Sparsity_parameter_col <- attr(x, "Sparsity_parameter_col")
  cat("Columns Sparsity Parameter: ")
  if (!is.null(Sparsity_parameter_col)) {
    if (length(Sparsity_parameter_col) > 5) {
      cat(head(Sparsity_parameter_col, 3), "...", tail(Sparsity_parameter_col, 2), "\n")
    } else {
      cat(Sparsity_parameter_col, "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print the Sparsity parameter for rows
  Sparsity_parameter_row <- attr(x, "Sparsity_parameter_row")
  cat("Rows Sparsity Parameter: ")
  if (!is.null(Sparsity_parameter_row)) {
    if (length(Sparsity_parameter_row) > 5) {
      cat(head(Sparsity_parameter_row, 3), "...", tail(Sparsity_parameter_row, 2), "\n")
    } else {
      cat(Sparsity_parameter_row, "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print GridPoints_u (vector)
  GridPoints_u <- attr(x, "GridPoints_u")
  cat("GridPoints_u: ")
  if (!is.null(GridPoints_u)) {
    if (length(GridPoints_u) > 5) {
      cat(head(GridPoints_u, 3), "...", tail(GridPoints_u, 2), "\n")
    } else {
      cat(GridPoints_u, "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print GridPoints_v (vector)
  GridPoints_v <- attr(x, "GridPoints_v")
  cat("GridPoints_v: ")
  if (!is.null(GridPoints_v)) {
    if (length(GridPoints_v) > 5) {
      cat(head(GridPoints_v, 3), "...", tail(GridPoints_v, 2), "\n")
    } else {
      cat(GridPoints_v, "\n")
    }
  } else {
    cat("NULL\n")
  }

  cat("-----------------------------------\n")
  cat("First few rows of the data:\n")
  print(head(x))

  invisible(x)  # Return the object invisibly
}
