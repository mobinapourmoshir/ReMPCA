############################### Calculating the optimal alpha for u using GCV ###############################
opt_alpha_u <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas, # A vector
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- length(alphas)
  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    return(list(GCV_u = Inf, opt.alpha_u = 0, opt_s.alpha_u = diag(n), GCVdf_u = data.frame(alphas, rep(Inf, n_iter))))
  } else {
    for (i in 1:n_iter) {
      S <- S_alphas_u[[i]]
      GCV_alpha <- 0

      power_result <- power_algo(data = X,
                                 n_var,
                                 ncol,
                                 sparse_tuning_result_u, # A fixed number
                                 sparse_tuning_result_v, # A vector (length = p)
                                 S_alpha_v = S_alphas_v, # A list (lenght = p)
                                 S_alpha_u = S, # A matrix
                                 sparse_tuning_type)
      u <- power_result[[2]]
      v <- t(X) %*% u
      GCV_alpha <- ( ((1/n) * (diag(nrow(S) - S))%*% (X %*% v) / norm_vec(v))^2 ) / ( 1 -  (1/n)* sum(diag(S)))^2
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV)]
    opt_s.alpha <- S_alphas_u[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV_u = GCV, opt.alpha_u = opt.alpha, opt_s.alpha_u = opt_s.alpha, GCVdf = data.frame(alphas, GCV))

    return(result)
  }
}


############################### Calculating the optimal alpha for v using conditional GCV ###############################
opt_alpha_v <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas, # A matrix
                        alpha_u,
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        Omegas_u,
                        sparse_tuning_type) {

  n_iter <- length(alphas)
  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    return(list(GCV_v = Inf, opt.alpha_v = 0, opt_s.alpha_v = diag(n), GCVdf_v = data.frame(alphas, rep(Inf, n_iter))))
  } else {
    for (i in 1:n_iter) {
      S <- S_alphas_u[[i]]
      Omega <- Omegas_u[[i]]
      alpha <- alphas[i]
      GCV_alpha <- 0
      m <- ncol(X)

      power_result <- power_algo(data = X,
                                 n_var,
                                 ncol,
                                 sparse_tuning_result_u, # A fixed number
                                 sparse_tuning_result_v, # A vector (length = p)
                                 S_alpha_v = S_alphas_v, # A list (lenght = p)
                                 S_alpha_u = S, # A matrix
                                 sparse_tuning_type)
      u <- power_result[[2]]
      v <- t(X) %*% u
      R <- (t(u) %*% Omega %*% u) / norm_vec(u)

      GCV_alpha <- ((1/m)*(norm_vec(( t(X)%*%u / norm_vec(u)) - v))^2) / ( 1 - (1/m) * (sum(diag(S)) / (1 + alpha_u*R) ))^2
      #setTxtProgressBar(pb, i + (k - 1) / n_var)
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV_v = GCV, opt.alpha_v = opt.alpha, opt_s.alpha_v = opt_s.alpha, GCVdf_v = data.frame(alphas, GCV))

    return(result)
  }
}

