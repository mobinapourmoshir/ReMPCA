##################### Print method for hybrid data #####################
#' @export
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



#' Coerce an object of class 'rdClass', 'fdClass', or 'imgClass' to class 'hdClass'
#'
#' @param x An object of class 'rdClass', 'fdClass', or 'imgClass'.
#' @param Smoothing_parameter Optional smoothing parameter for the row direction.
#'                            If NULL, taken from input attributes or defaulted in hdClass().
#' @param Sparsity_parameter Optional sparsity parameter for the row direction.
#'                           If NULL, taken from input attributes or defaulted in hdClass().
#' @param argval Optional grid points for rows. If NULL, taken from input attributes or defaulted.
#'
#' #' @details
#' This function returns a modified version of the input object with an updated \code{"Sparsity_parameter"} attribute.
#' Due to R's copy-on-modify semantics for S3 objects, the user must reassign the object
#'
#' @return An object of class 'hdClass'.
#' #' @examples
#' fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
#' hd_obj <- as.hdClass(fd_obj, Sparsity_parameter = 0:5)
#' attr(hd_obj, "Sparsity_parameter")
#' @export
as.hdClass <- function(x,
                       Smoothing_parameter = NULL,
                       Sparsity_parameter = NULL,
                       argval = NULL) {
  # Validate class
  if (!(inherits(x, "rdClass") ||
        inherits(x, "fdClass") ||
        inherits(x, "imgClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', or 'imgClass'!")
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


  datalist <- list(x)
  # Call hdClass using the list
  hdClass(hdlist = datalist,
          argval = argval,
          Smoothing_parameter = Smoothing_parameter,
          Sparsity_parameter = Sparsity_parameter)
}

#' Set Sparsity Tuning Parameter
#'
#' Assigns a sparsity tuning parameter to an object of class \code{rdClass}, \code{fdClass}, \code{imgClass}, or \code{hdClass}.
#'
#' @param obj An object of class \code{rdClass}, \code{fdClass}, \code{imgClass}, or \code{hdClass}.
#' @param Sparsity_parameter A numeric vector of non-negative integers indicating sparsity tuning levels.
#'   If \code{NULL}, a default sequence will be generated based on the object's dimensions.
#'
#' @details
#' This function returns a modified version of the input object with an updated \code{"Sparsity_parameter"} attribute.
#' Due to R's copy-on-modify semantics for S3 objects, the user must reassign the object:
#' \preformatted{
#'   obj <- setSparsityParameter(obj, c(0, 2, 4))
#' }
#' The object will not be updated in-place unless reassigned.
#'
#' @return The modified object with updated sparsity parameters.
#'
#' @examples
#' fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
#' fd_obj <- setSparsityParameter(fd_obj, 0:5)
#' attr(fd_obj, "Sparsity_parameter")
#'
#' @export

setSparsityParameter <- function(obj, Sparsity_parameter) {
  if (!(inherits(obj, "rdClass") ||
        inherits(obj, "fdClass") ||
        inherits(obj, "hdClass") ||
        inherits(obj, "imgClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', 'imgClass', or 'hdClass'.")
  }

  # Infer dimension for validation
  dim_target <- if (inherits(obj, "hdClass")) nrow(obj) else ncol(obj)

  if (is.null(Sparsity_parameter)) {
    if (dim_target <= 15) {
      Sparsity_parameter <- 0:(dim_target - 1)
    } else {
      extra_vals <- unique(c(0:3, 2^(0:floor(log2(dim_target - 1))),
                             dim_target - 1))
      Sparsity_parameter <- sort(unique(extra_vals[extra_vals <=
                                                     (dim_target - 1)]))
    }
  }

  if (!is.numeric(Sparsity_parameter) ||
      any(Sparsity_parameter < 0) ||
      any(Sparsity_parameter != floor(Sparsity_parameter))) {
    stop("Sparsity_parameter must be a vector of non-negative integers.")
  }

  attr(obj, "Sparsity_parameter") <- as.vector(Sparsity_parameter)
  return(obj)
}

#' Set Smoothness Tuning Parameter
#'
#' Assigns a smoothing tuning parameter to an object of class \code{rdClass}, \code{fdClass}, \code{imgClass}, or \code{hdClass}.
#'
#' @param obj An object of class \code{rdClass}, \code{fdClass}, \code{imgClass}, or \code{hdClass}.
#' @param Smoothing_parameter A numeric value or vector representing the smoothing parameter(s)
#'   to assign to the object.
#'
#' @details
#' This function updates the \code{"Smoothing_parameter"} attribute of the given object.
#' Since R uses copy-on-modify semantics for S3 objects, users must reassign the object after calling this function:
#' \preformatted{
#'   obj <- setSmoothnessParameter(obj, c(0.01, 0.1, 1))
#' }
#' Without reassignment, the original object remains unchanged.
#'
#' @return The modified object with updated \code{Smoothing_parameter} attribute.
#'
#' @examples
#' fd_obj <- fdClass(matrix(rnorm(100), 10, 10))
#' fd_obj <- setSmoothnessParameter(fd_obj, c(0.01, 0.1, 1))
#' attr(fd_obj, "Smoothing_parameter")
#'
#' @export
setSmoothnessParameter <- function(obj, Smoothing_parameter) {
  if (!(inherits(obj, "rdClass") ||
        inherits(obj, "fdClass") ||
        inherits(obj, "hdClass") ||
        inherits(obj, "imgClass"))) {
    stop("Input must be of class 'rdClass', 'fdClass', 'imgClass', or 'hdClass'.")
  }

  attr(obj, "Smoothing_parameter") <- as.vector(Smoothing_parameter)
  return(obj)
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
    end_idx <- as.numeric(start_idx + ncol_vec[i] - 1)
    mat <- scaled_matrix[, start_idx:end_idx, drop = FALSE]
    weights[i] <- 1 / sqrt(mean(diag(var(mat))))
    scaled_matrix[, start_idx:end_idx] <- weights[i] * mat
    start_idx <- end_idx + 1
  }
  return(list(
    scaled_hd = scaled_matrix,
    weights = weights
  ))
}
