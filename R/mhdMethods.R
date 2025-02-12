print.hd <- function(x, ...) {
  cat("Hybrid Data (hdClass) Object\n")
  cat("-----------------------------------\n")
  cat("Dimensions: ", dim(x)[1], " x ", dim(x)[2], "\n", sep = "")

  # Get attributes
  smoothing_param <- attr(x, "Smoothing_parameter")
  two_way_smoothness <- attr(x, "two_way_smoothness")
  Sparsity_parameter_col <- attr(x, "Sparsity_parameter_col")
  Sparsity_parameter_row <- attr(x, "Sparsity_parameter_row")

  # Check if the data is functional or regular
  is_functional <- !identical(smoothing_param, 0) || !identical(two_way_smoothness, 0)

  if (is_functional) {
    cat("Data Type: Functional Data\n")
    # Print GridPoints_u and GridPoints_v
    GridPoints_u <- attr(x, "GridPoints_u")
    GridPoints_v <- attr(x, "GridPoints_v")

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
  } else {
    cat("Data Type: Regular Data\n")
  }

  # Print Smoothing Parameter
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

  # Print Two-Way Smoothing Parameter
  cat("Two-Way Smoothing Parameter: ")
  if (!is.null(two_way_smoothness)) {
    if (length(two_way_smoothness) > 5) {
      cat(head(two_way_smoothness, 3), "...", tail(two_way_smoothness, 2), "\n")
    } else {
      cat(two_way_smoothness, "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print Sparsity Parameters
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

  cat("-----------------------------------\n")
  cat("First few rows and columns of the data:\n")

  # Extract and print only the first 5 rows and 10 columns
  rows_to_show <- min(5, nrow(x))
  cols_to_show <- min(10, ncol(x))
  print(x[1:rows_to_show, 1:cols_to_show])

  invisible(x)  # Return the object invisibly
}
