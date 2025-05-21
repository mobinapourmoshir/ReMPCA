##################### Print method for hybrid data #####################
print.hdClass <- function(x, ...) {
  cat("Hybrid Data (hdClass) Object\n")
  cat("===================================\n")
  cat("Dimensions           : ", dim(x)[1], " rows x ",
      dim(x)[2], " columns\n", sep = "")
  cat("Number of Variables  : ", attr(x, "n_var"), "\n")
  cat("Columns per Variable : ", paste(attr(x, "ncol"), collapse = " | "), "\n")
  cat("-----------------------------------\n")

  # Print GridPoints_u
  cat("Grid Points (Rows)   : ")
  gp_u <- attr(x, "GridPoints_u")
  if (!is.null(gp_u)) {
    if (length(gp_u) > 5) {
      cat(paste0(head(gp_u, 3), collapse = ", "), ", ..., ",
          paste0(tail(gp_u, 2), collapse = ", "), "\n")
    } else {
      cat(paste(gp_u, collapse = ", "), "\n")
    }
  } else {
    cat("NULL\n")
  }

  # Print GridPoints_v
  cat("Grid Points (Columns):\n")
  gp_v <- attr(x, "GridPoints_v")
  if (!is.null(gp_v)) {
    for (i in seq_along(gp_v)) {
      cat(sprintf("  - Variable %d: ", i))
      gpv_i <- gp_v[[i]]
      if (!is.null(gpv_i)) {
        if (length(gpv_i) > 5) {
          cat(paste0(head(gpv_i, 3), collapse = ", "), ", ..., ",
              paste0(tail(gpv_i, 2), collapse = ", "), "\n")
        } else {
          cat(paste(gpv_i, collapse = ", "), "\n")
        }
      } else {
        cat("NULL\n")
      }
    }
  } else {
    cat("  NULL\n")
  }

  # Print Smoothing Parameters
  cat("Smoothing Parameter (Row): ")
  s_row <- attr(x, "Smoothing_parameter")
  if (!is.null(s_row)) {
    if (length(s_row) > 5) {
      cat(paste0(head(s_row, 3), collapse = ", "), ", ..., ",
          paste0(tail(s_row, 2), collapse = ", "), "\n")
    } else {
      cat(paste(s_row, collapse = ", "), "\n")
    }
  } else {
    cat("NULL\n")
  }

  cat("Smoothing Parameters (Col):\n")
  s_col <- attr(x, "Smoothing_parameter_col")
  if (!is.null(s_col)) {
    for (i in seq_along(s_col)) {
      cat(sprintf("  - Variable %d: ", i))
      cat(paste(s_col[[i]], collapse = ", "), "\n")
    }
  } else {
    cat("  NULL\n")
  }

  # Print Sparsity Parameters
  cat("Sparsity Parameter (Row): ")
  sp_row <- attr(x, "Sparsity_parameter")
  if (!is.null(sp_row)) {
    if (length(sp_row) > 5) {
      cat(paste0(head(sp_row, 3), collapse = ", "), ", ..., ",
          paste0(tail(sp_row, 2), collapse = ", "), "\n")
    } else {
      cat(paste(sp_row, collapse = ", "), "\n")
    }
  } else {
    cat("NULL\n")
  }

  cat("Sparsity Parameters (Col):\n")
  sp_col <- attr(x, "Sparsity_parameter_col")
  if (!is.null(sp_col)) {
    for (i in seq_along(sp_col)) {
      cat(sprintf("  - Variable %d: ", i))
      cat(paste(sp_col[[i]], collapse = ", "), "\n")
    }
  } else {
    cat("  NULL\n")
  }

  cat("===================================\n")
  cat("First few rows and columns of the data:\n")
  rows_to_show <- min(5, dim(x)[1])
  cols_to_show <- min(5, dim(x)[2])
  print(x[1:rows_to_show, 1:cols_to_show])
  invisible(x)
}

#' Check if an object is of class 'hdClass'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'hdClass', FALSE otherwise.
#' @export
is.hdClass <- function(x) {
  inherits(x, "hdClass")
}


#' Custom `$` operator for hdClass
#' Returns a clean data.frame when using hd_obj$matrix
#'
#' @param x An object of class 'hdClass'
#' @param name The field to extract
#' @export
`$.hdClass` <- function(x, name) {
  if (name == "matrix") {
    return(as.data.frame(unclass(x)))  # convert matrix to plain data.frame
  } else {
    stop(sprintf("Unknown field '%s'. Only 'matrix' is supported for hdClass."),
         call. = FALSE)
  }
}



#' Coerce an object of class 'rdClass', 'fdClass', or 'hdClass' to class 'hdClass'
#'
#' @param x An object of class 'rdClass', 'fdClass', or 'hdClass'.
#' @param Smoothing_parameter Optional smoothing parameter for the row direction.
#'                            If NULL, taken from input attributes or defaulted in hdClass().
#' @param Sparsity_parameter Optional sparsity parameter for the row direction.
#'                           If NULL, taken from input attributes or defaulted in hdClass().
#' @param argval Optional grid points for rows. If NULL, taken from input attributes or defaulted.
#'
#' @return An object of class 'hdClass'.
#' @export
as.hdClass <- function(x,
                       Smoothing_parameter = NULL,
                       Sparsity_parameter = NULL,
                       argval = NULL) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "fdClass") ||
        inherits(x, "hdClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', or 'hdClass'")
  }

  # Extract or override attributes
  if (is.null(Smoothing_parameter)) {
    Smoothing_parameter <- attr(x, "Smoothing_parameter")
  }else{
    Smoothing_parameter <- Smoothing_parameter}
  if (is.null(Sparsity_parameter)) {
    Sparsity_parameter <- attr(x, "Sparsity_parameter")
  }else{
    Sparsity_parameter <- Sparsity_parameter}
  if (is.null(argval)) {
    argval <- attr(x, "GridPoints_u")
  }else{
    argval <- argval}

  # Construct hdlist from $matrix using hdClass attributes
  if (inherits(x, "hdClass")) {
    data <- x$matrix
    ncol <- attr(x, "ncol")
    n_var <- attr(x, "n_var")
    Smoothing_parameter_col <- attr(x, "Smoothing_parameter_col")
    Sparsity_parameter_col <- attr(x, "Sparsity_parameter_col")
    GridPoints_v <- attr(x, "GridPoints_v")
    attr(x, "GridPoints_u") <- argval
    attr(x, "Sparsity_parameter") <- Sparsity_parameter
    attr(x, "Smoothing_parameter") <- Smoothing_parameter
    x

  } else {
    datalist <- list(x)
    # Call hdClass using the list
    hdClass(hdlist = datalist,
            argval = argval,
            Smoothing_parameter = Smoothing_parameter,
            Sparsity_parameter = Sparsity_parameter)
  }
}


#' Compute Scaling Weights for `hdClass` Object
#'
#' This function computes scaling weights for each component (column) in a `hdClass` object.
#' The goal is to normalize the contributions of different functional variables by
#' accounting for their scale, using the inverse of total variance (integrated over domain).
#'
#' @param hd_obj A hybrid data object of class `hdClass`.
#'
#' @return A numeric vector of weights (length equals the number of variables in the hybrid data).
#'
#' @examples
#' weights <- get_hd_scaling_weights(hd_obj)
scale_hd <- function(hd_obj) {
  if (!inherits(hd_obj, "hdClass")) stop("Input must be of class 'hdClass'")

  n_var <- attr(hd_obj, "n_var")
  ncol_vec <- attr(hd_obj, "ncol")  # number of columns per variable
  weights <- numeric(n_var)

  scaled_matrix <- hd_obj$matrix
  start_idx <- 1

  for (i in 1:n_var) {
    end_idx <- start_idx + ncol_vec[i] - 1
    mat <- scaled_matrix[, start_idx:end_idx, drop = FALSE]
    weights[i] <- 1 / sqrt(mean(diag(var(mat))))
    scaled_matrix[, start_idx:end_idx] <- weights[i] * mat
    start_idx <- end_idx + 1
  }

  hd_scaled <- hd_obj
  hd_scaled$matrix <- scaled_matrix

  return(list(
    scaled_hd = hd_scaled,
    weights = weights
  ))
}
