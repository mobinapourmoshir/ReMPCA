############################### Calculating the norm of a vector ###############################
norm_vec <- function(x) sqrt(sum(x^2))


############################### Process bar indexing ###############################
ordinal <- function(i) {
  if (i == 1) {
    return(paste0(i, "st"))
  } else if (i == 2) {
    return(paste0(i, "nd"))
  } else if (i == 3) {
    return(paste0(i, "rd"))
  }
  else {
    return(paste0(i, "th"))
  }
}
