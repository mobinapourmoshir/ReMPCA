############################### Sparse penalty for coefficients ###############################
# Lemma 2 (Sparse PCA via regularized low rank matrix approximation by Huang)
# y is the coefficients (u's)
csparse_pen_fun <- function(y,tuning_parameter, type,alpha = 3.7) {
  y_sorted <- sort(abs(y))
  lambda = y_sorted[tuning_parameter]
  if (tuning_parameter == 0) {
    return(y)
  }
  if (type == "soft") {
    return(sign(y) * pmax(abs(y) - lambda, 0))
  }
  else if (type == "hard") {
    return(ifelse(abs(y) > lambda, y, 0))
  }
  else if (type == "SCAD") {
    res <- ifelse(abs(y) <= 2 * lambda,
                  sign(y) * pmax(abs(y) - lambda, 0),
                  ifelse(abs(y) <= alpha * lambda,
                         ((alpha - 1) * y - sign(y) * alpha * lambda) / (alpha - 2),
                         y))
    return(res)
  }
}
