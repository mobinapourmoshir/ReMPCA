#' Read ReMPCA Example Data
#'
#' Loads one of the datasets included with the ReMPCA package.
#'
#' @param data Character string specifying the dataset to load.
#'   Available options are `"bike_day"` and `"bike_hour"`.
#'
#' @return The requested dataset as a data frame.
#'
#' @examples
#' bike_day <- read.data("bike_day")
#' bike_hour <- read.data("bike_hour")
#'
#' @export
read.data <- function(data = c("bike_day", "bike_hour")) {

  data <- match.arg(data)

  e <- new.env()

  utils::data(
    list = data,
    package = "ReMPCA",
    envir = e
  )

  e[[data]]
}
