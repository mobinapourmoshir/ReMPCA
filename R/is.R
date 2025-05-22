#' Check if an object is of class 'fdClass'
#'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'fdClass', FALSE otherwise.
#'
#' @export
is.fdClass <- function(x) {
  inherits(x, "fdClass")
}

#' Check if an object is of class 'rdClass'
#'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'rdClass', FALSE otherwise.
#'
#' @export
is.rdClass <- function(x) {
  inherits(x, "rdClass")
}

#' Check if an object is of class 'hdClass'
#'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'hdClass', FALSE otherwise.
#'
#' @export
is.hdClass <- function(x) {
  inherits(x, "hdClass")
}

#' Check if an object is of class 'imgClass'
#'
#' @param x An object to test.
#' @return Logical; TRUE if the object inherits from class 'imgClass', FALSE otherwise.
#'
#' @export
is.imgClass <- function(x) {
  inherits(x, "imgClass")
}

