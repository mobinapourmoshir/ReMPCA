############################### get.pen function ###############################
get.pen <- function(td, alpha=0) {
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
  return(S.alpha=S.alpha)
}
