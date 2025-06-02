#' @title Hybrid Data Constructor (`hdClass`)
#'
#' @description Constructs an object of class `hdClass`, representing hybrid data composed of multiple functional and/or raw variables.
#'
#' @details
#' The `hdClass` is an S3 class for representing *hybrid data*, consisting of a list of functional (`fdClass`), raw (`rdClass`), or image (`imgClass`) data objects.
#'
#' - The input list (`hdlist`) can contain a mix of `fdClass`, `rdClass`, and `imgClass` objects.
#' - All objects must have the **same number of observations** (i.e., same number of rows).
#' - Image objects (`imgClass`) are interpreted as either functional or raw based on their internal class inheritance.
#'
#' The resulting object is a matrix with column-wise concatenation of the input matrices and includes the following attributes:
#'
#' - `GridPoints_u`: Row grid points.
#' - `Smoothing_parameter`: Row smoothing parameter(s).
#' - `Sparsity_parameter`: Row sparsity parameter(s).
#' - `n_var`: Number of variables (i.e., number of objects in `hdlist`).
#' - `ncol`: Number of columns in each variable.
#' - `Smoothing_parameter_col`: List of smoothing parameters for each variable.
#' - `Sparsity_parameter_col`: List of sparsity parameters for each variable.
#' - `GridPoints_v`: List of grid points (columns) for each variable.
#' - `variable_types`: Character vector of type labels for each variable: `"hd"` (functional) or `"rd"` (raw).
#'
#' @param hdlist A list of `fdClass`, `rdClass`, or `imgClass` objects.
#' @param argval Optional numeric vector of grid points along the rows. If `NULL`, it defaults to a uniform grid over [0, 1].
#' @param Smoothing_parameter Smoothing parameter(s) for the rows:
#'   - If 0, no smoothing is applied.
#'   - If a numeric value or vector, it is used directly.
#'   - If `NULL`, defaults to `2^seq(-30, 5, length.out = 10)` for tuning.
#' @param Sparsity_parameter Sparsity parameter(s) for the rows:
#'   - If 0, no sparsity is applied.
#'   - If a numeric vector, it is used for tuning.
#'   - If `NULL`, a suitable default sequence is generated.
#'
#' @return An object of class `hdClass` (a matrix) with several hybrid-aware attributes.
#'
#' @examples
#' fd_obj <- fdClass(matrix(rnorm(100), nrow = 10))
#' rd_obj <- rdClass(matrix(rnorm(100), nrow = 10))
#'
#' hd_obj <- hdClass(
#'   hdlist = list(fd_obj, rd_obj),
#'   argval = seq(0, 1, length.out = 10),
#'   Smoothing_parameter = 0.5,
#'   Sparsity_parameter = 2
#' )
#'
#' print(hd_obj)
#' attr(hd_obj, "variable_types")  # Shows "hd", "rd", etc.
#'
#' @export

hdClass <- function(hdlist,
                    argval = NULL,
                    Smoothing_parameter = 0,
                    Sparsity_parameter = 0) {

  if (!is.list(hdlist)) {
    hdlist <- list(hdlist)
  }

  # Validate input: must be fdClass, rdClass, or imgClass
  if (!all(sapply(hdlist, function(obj) any(class(obj) %in% c("rdClass", "fdClass", "imgClass"))))) {
    stop("All elements in the list must be of class 'fdClass', 'rdClass', or 'imgClass'.")
  }

  # Ensure all matrices have the same number of observations
  nrows <- sapply(hdlist, function(obj) nrow(obj))
  if (!all(nrows == nrows[1])) {
    stop("Error: Not all matrices have the same number of rows!")
  }

  hd <- do.call(cbind, hdlist)

  ####### Grid Points for rows #######
  if (!is.null(argval)) {
    if (length(argval) != nrow(hd)) {
      warning("There should be an equal number of grid points and rows for hybrid data. 'argval' is set to NULL!")
      argval <- NULL
    }
  }

  GridPoints_u <- if (!is.null(argval)) argval else seq(from = 1 / nrow(hd), to = 1, length.out = nrow(hd))
  attr(hd, "GridPoints_u") <- GridPoints_u

  ####### Smoothing_parameter #######
  if (!is.null(Smoothing_parameter) &&
      !is.numeric(Smoothing_parameter) &&
      !identical(Smoothing_parameter, 0)) {
    stop("Smoothing_parameter must be a numeric value, numeric vector, 0, or NULL.")
  }

  if (is.null(Smoothing_parameter)) {
    Smoothing_parameter <- 2^seq(-30, 5, length.out = 10)
  }
  attr(hd, "Smoothing_parameter") <- Smoothing_parameter

  ####### Sparsity_parameter #######
  if (!is.null(Sparsity_parameter)) {
    if (!is.numeric(Sparsity_parameter) ||
        any(Sparsity_parameter < 0) ||
        any(Sparsity_parameter != floor(Sparsity_parameter))) {
      stop("Sparsity_parameter must be a vector of non-negative integers.")
    }
    if (any(Sparsity_parameter > nrow(hd) - 1)) {
      stop("All elements of Sparsity_parameter must be between 0 and nrow(hd) - 1.")
    }
  } else {
    if (nrow(hd) <= 15) {
      Sparsity_parameter <- 0:(nrow(hd) - 1)
    } else {
      extra_vals <- unique(c(0:3, 2^(0:floor(log2(nrow(hd) - 1))), nrow(hd) - 1))
      Sparsity_parameter <- sort(unique(extra_vals[extra_vals <= (nrow(hd) - 1)]))
    }
  }

  attr(hd, "Sparsity_parameter") <- Sparsity_parameter
  attr(hd, "n_var") <- length(hdlist)
  attr(hd, "ncol") <- as.data.frame(sapply(hdlist, dim))[2, ]

  ####### Column-wise parameter attributes #######
  attr(hd, "Smoothing_parameter_col") <- lapply(hdlist, function(obj) attr(obj, "Smoothing_parameter"))
  attr(hd, "GridPoints_v") <- lapply(hdlist, function(obj) attr(obj, "GridPoints_v"))
  attr(hd, "Sparsity_parameter_col") <- lapply(hdlist, function(obj) attr(obj, "Sparsity_parameter"))

  ####### Column type: hd or rd #######
  column_class <- sapply(hdlist, function(obj) {
    class_type <- class(obj)
    if ("imgClass" %in% class_type) {
      if ("fdClass" %in% class_type) {
        return("hd")
      } else if ("rdClass" %in% class_type) {
        return("rd")
      } else {
        stop("imgClass must also inherit either fdClass or rdClass.")
      }
    } else if ("fdClass" %in% class_type) {
      return("hd")
    } else {
      return("rd")
    }
  })
  attr(hd, "variable_types") <- column_class  # vector of "hd"/"rd" types

  class(hd) <- "hdClass"
  return(hd)
}
