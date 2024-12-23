#' @title  Define a Set of Multivariate Functional Data objects
#'
#' #' @description
#' The `mfd` class represents a set of multivariate functional data.
#' Functional data Objects are constructed using matrices, with columns representing grid points and rows indicating observations.
#'
#'
#'


mfd <- function(argval = NULL, mvgfd_obj, centerfns = TRUE, num_pcs = 1) {

  #' @param argval A list of numeric vectors of argument values at which the `mfd` object is to be evaluated
  if (!is.list(argval) && !is.null(argval)) {
    stop("`argval` must be a list of numeric vectors or NULL")
  }

  #' @param mvfd_obj List of matrices or arrays (Multivariate Grid Functional Data )
  if (!(is.list(mvgfd_obj) || is.matrix(mvgfd_obj))) {
    stop("Input data must be a list of matrices and arrays or just one matrix!")
  }

  #' @param centerfns logical. If TRUE  the input data undergoes centralization or demeaning prior to being processed by the function.
  #' @param num_pcs  The number of PCs. The default is one (The first principal component only). But the user is able to see higher PCs
  structure(
    list(
      argval = argval,
      mvfd_obj = mvgfd_obj,
      centerfns = centerfns,
      num_pcs = num_pcs
    ),
    class = "mvgfd"
  )
}

# Check for a valid data
check_data.mvgfd <- function(object) {
  #' @description Check for validity of the data in the `mvgfd` object

  mvfd_obj <- object$mvfd_obj

  # Check if the list contains only matrices or arrays
  if (is.list(mvfd_obj)) {
    valid_entries <- all(sapply(mvfd_obj, function(x) is.matrix(x) || is.array(x)))
    if (!valid_entries) {
      stop("The list contains entries other than matrices or arrays.")
    }
  }

  # Check if the list contains arrays or handles univariate case
  if (is.matrix(mvfd_obj)) {
    all_matrices <- "Univariate case!"
  } else {
    all_matrices <- all(sapply(mvfd_obj, is.matrix))
  }

  # Check for NA or non-numeric entries in matrices or arrays
  contains_na_non_numeric <- any(sapply(mvfd_obj, function(x) any(is.na(x) | !is.numeric(x))))

  # Display report on the data
  cat("Report on the data:\n")
  cat("All matrices: ", all_matrices, "\n")
  cat("Contains NA or non-numeric entries: ", contains_na_non_numeric, "\n")

  return(list(all_matrices = all_matrices, contains_na_non_numeric = contains_na_non_numeric))
}

# Create getter functions for fields
get_argval <- function(object) {
  object$argval
}

get_mvfd_obj <- function(object) {
  object$mvfd_obj
}

get_centerfns <- function(object) {
  object$centerfns
}

get_num_pcs <- function(object) {
  object$num_pcs
}

############### TEST ###############

library(fda)
df <- gait
datagait <- list(t(df[,,1]), t(df[,,2]))

# Create an `mvgfd` object
data_checker_new <- mvgfd(mvgfd_obj = datagait, centerfns = FALSE)

# Check the data
result2 <- check_data.mvgfd(data_checker_new)

# Test with invalid data
data_checker_new2 <- mvgfd(mvgfd_obj = c(2, 5), centerfns = FALSE)
check_data.mvgfd(data_checker_new2)
