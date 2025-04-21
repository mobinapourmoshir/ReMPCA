############################### Calculating the optimal alpha for u using GCV ###############################
opt_alpha_u <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas, # A vector
                        Omegas_u,
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
    Omega_u <- Omega_u[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV_u = GCV, opt.alpha_u = opt.alpha, opt_s.alpha_u = opt_s.alpha, GCVdf = data.frame(alphas, GCV), Omega_u = Omega_u)

    return(result)
  }
}


############################### Calculating the optimal alpha for v using conditional GCV ###############################
opt_alpha_v <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v,  # List of length = n_iter, each element is list of length = n_var
                        S_alphas_u,  # Matrix (n × n)
                        alphas,      # Matrix of alpha_v values (n_iter × n_var)
                        alpha_u,
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        Omega_u,
                        sparse_tuning_type) {

  n_iter <- nrow(alphas)
  n <- nrow(X)
  m <- ncol(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    return(list(GCV_v = Inf,
                opt.alpha_v = rep(0, n_var),
                opt_s.alpha_v = replicate(n_var, diag(m / n_var), simplify = FALSE),
                GCVdf_v = data.frame(alphas, rep(Inf, n_iter))))
  } else {
    for (i in 1:n_iter) {
      S_alpha_v <- S_alphas_v[[i]]  # list of S_k
      GCV_alpha <- 0
      trace_total <- 0

      power_result <- power_algo(data = X,
                                 n_var = n_var,
                                 ncol = ncol,
                                 sparse_tuning_result_u = sparse_tuning_result_u,
                                 sparse_tuning_result_v = sparse_tuning_result_v,
                                 S_alpha_v = S_alpha_v,
                                 S_alpha_u = S_alphas_u,
                                 sparse_tuning_type = sparse_tuning_type)

      u <- power_result[[2]]  # fixed left singular vector
      df <- X

      for (k in 1:n_var) {
        m_k <-as.integer(ncol[k])
        X_k <- df[, 1:m_k]
        df <- df[, -(1:m_k), drop = FALSE]
        S_k <- S_alpha_v[[k]]

        # Residual norm after smoothing
        residual <- (diag(m_k) - S_k) %*% (t(X_k) %*% u)
        GCV_alpha <- GCV_alpha + (1 / m_k) * (norm_vec(residual)^2)
        trace_total <- trace_total + sum(diag(S_k))
      }

      # Final GCV with global denominator
      GCV[i] <- GCV_alpha / (1 - (trace_total / m))^2
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]

    return(list(GCV_v = GCV,
                opt.alpha_v = opt.alpha,
                opt_s.alpha_v = opt_s.alpha,
                GCVdf_v = data.frame(alphas, GCV)))
  }
}

