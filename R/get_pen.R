############################### get.pen function ###############################
get.pen <- function(td, alpha, type = "Second_order") {

  if(type == "Second_order"){
    m = length(td);
    h = td[2:m] - td[1:(m-1)];
    Q = matrix(0, m, m-1);
    R = matrix(0, m-1, m-1);

    for(k in 2:(m-1))
    {
      Q[k-1,k] = 1/h[k-1];
      Q[k,k] = -1/h[k-1] - 1/h[k];
      Q[k+1,k] = 1/h[k]
    }

    for(j in 2:(m-2))
    {
      R[j,j] = 1/3 * (h[j-1] + h[j]);
      R[j,j+1] = 1/6 * h[j];
      R[j+1,j] = 1/6 * h[j]
    }
    R[m-1,m-1] = 1/3 * (h[m-2] + h[m-1]);
    s <- solve(R[2:(m-1), 2:(m-1)]) %*% t(Q[1:m, 2:(m-1)]);
    OMEGA = Q[1:m, 2:(m-1)] %*% s;
    EIG.O <- eigen(OMEGA); GAMMA=EIG.O$vectors; LAMBDA=diag(EIG.O$values);
    S.alpha <- GAMMA%*%diag((1/(1+alpha*diag(LAMBDA))))%*%t(GAMMA);

    if(alpha == 0){
      S.alpha = diag(nrow(S.alpha))}

    return(S.alpha=S.alpha)

  } else if(type == "First_order") {
    m = length(td)
    h = td[2:m] - td[1:(m-1)]

    Q = matrix(0, m, m-1)
    for (k in 2:m) {
      Q[k-1, k-1] = -1/h[k-1]
      Q[k, k-1] = 1/h[k-1]
    }

    R = matrix(0, m-1, m-1)

    for (j in 1:(m-1)) {
      if (j == 1) {
        R[j, j] = h[j] / 3
      } else if (j == (m-1)) {
        R[j, j] = h[j-1] / 3
      } else {
        R[j, j] = (h[j-1] + h[j]) / 3
        R[j, j+1] = h[j] / 6
        R[j+1, j] = h[j] / 6
      }
    }

    s <- solve(R) %*% t(Q)
    OMEGA = Q %*% s
    EIG.O <- eigen(OMEGA)
    GAMMA = EIG.O$vectors
    LAMBDA = diag(EIG.O$values)
    S.alpha = GAMMA %*% diag((1/(1 + alpha * diag(LAMBDA)))) %*% t(GAMMA)

    if(alpha == 0){
      S.alpha = diag(nrow(S.alpha))}

    return(S.alpha=S.alpha)

  } else if (type == "Indicator"){
    m = length(td)
    h = td[2:m] - td[1:(m-1)]

    # Q matrix for first-order differences with threshold penalty
    Q = matrix(0, m, m-1)
    for (k in 2:m) {
      diff_val = 1 / h[k-1]
      Q[k-1, k-1] = -diff_val
      Q[k, k-1] = diff_val
    }

    # R matrix for the integration weights
    R = matrix(0, m-1, m-1)

    for (j in 1:(m-1)) {
      if (j == 1) {
        R[j, j] = h[j] / 3
      } else if (j == (m-1)) {
        R[j, j] = h[j-1] / 3
      } else {
        R[j, j] = (h[j-1] + h[j]) / 3
        R[j, j+1] = h[j] / 6
        R[j+1, j] = h[j] / 6
      }
    }

    s <- solve(R) %*% t(Q)
    OMEGA = Q %*% s
    EIG.O <- eigen(OMEGA)
    GAMMA = EIG.O$vectors
    LAMBDA = diag(EIG.O$values)

    # Indicator
    penalty = diag((1 / (1 + alpha * diag(LAMBDA))))
    for (i in 2:m) {
      if (abs(penalty[i-1, i-1]) > thrs) {
        penalty[i-1, i-1] <- penalty[i-1, i-1]^2
      } else {
        penalty[i-1, i-1] <- 0
      }
    }

    S.alpha = GAMMA %*% penalty %*% t(GAMMA)

    if(alpha == 0){
      S.alpha = diag(nrow(S.alpha))}

    return(S.alpha=S.alpha)
  }
}
