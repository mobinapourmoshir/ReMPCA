#' @importFrom Matrix bdiag
#' @importFrom utils  txtProgressBar setTxtProgressBar
utils::globalVariables(c("GridPoints_v", "alpha_v", "smoothness_type", "pb"))
############### Calculating the optimal alpha for u using GCV ###############
opt_alpha_u <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas_u,   # A vector of alpha_u
                        alpha_v,    # A vector
                        Omega_v,    # A list associated with alpha_v
                        Omegas_u,   # List of length, each element is matrix
                        thresh,
                        maxit,
                        conditional,
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- length(alphas_u)
  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas_u == 0)) {
    return(list(GCV_u = Inf,
                opt.alpha_u = 0,
                opt_s.alpha_u = diag(n),
                GCVdf = data.frame(alphas_u, rep(Inf, n_iter)),
                Omega_u = diag(n),
                opt_alpha_Omega_u = 0 * diag(n)))
  } else {
    alpha_Omega_v <- bdiag(Map(function(a, M) a * M, alpha_v, Omega_v))
    for (i in 1:n_iter) {
      GCV_alpha <- 0
      S <- S_alphas_u[[i]]
      Omega_u <- Omegas_u[[i]]
      alpha_Omega_u <- as.matrix(alphas_u[i] * Omega_u)

      power_result <- power_algo(data = X,
                                 n_var,
                                 ncol,
                                 thresh = thresh,
                                 maxit = maxit,
                                 sparse_tuning_result_u, # A fixed number
                                 sparse_tuning_result_v, # A vector (length = p)
                                 S_alpha_v = S_alphas_v, # A list (lenght = p)
                                 S_alpha_u = S, # A matrix
                                 conditional = conditional,
                                 alpha_Omega_v = alpha_Omega_v,
                                 alpha_Omega_u = alpha_Omega_u,
                                 sparse_tuning_type = sparse_tuning_type)

      u <- power_result[[2]]
      v <- t(X) %*% u
      v <- as.vector(v)
      #Omega_v <- as.matrix(Omega_v)
      #storage.mode(Omega_v) <- "numeric"

      norm_v2 <- as.numeric(norm_vec(v)^2)
      alphaR_v <- as.numeric(t(v) %*% alpha_Omega_v %*% v / norm_v2)

      GCV_alpha <- ( (1/n) *(norm_vec(((as.matrix(X) %*%v) / norm_vec(v)) - u ))^2 ) /
        ( 1 -  (1/n)* (sum(diag(S)))/(1 + alphaR_v))^2
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas_u[which.min(GCV)]
    opt_s.alpha <- S_alphas_u[[which.min(GCV)]]
    Omega_u <- Omegas_u[[which.min(GCV)]]
    opt_alpha_Omega_u <- opt.alpha * Omega_u

    #close(pb)
    result <- list(GCV_u = GCV,
                   opt.alpha_u = opt.alpha,
                   opt_s.alpha_u = opt_s.alpha,
                   GCVdf = data.frame(alphas_u, GCV),
                   Omega_u = Omega_u,
                   opt_alpha_Omega_u = opt_alpha_Omega_u)

    return(result)
  }
}

########## Calculating the optimal alpha for v using conditional GCV ##########
opt_alpha_v <- function(X,
                        n_var,
                        ncol,
                        thresh,
                        maxit,
                        conditional,
                        S_alphas_v,   # List of length = n_iter, each element is list of length = n_var
                        S_alphas_u,   # Matrix (n × n)
                        alphas_v,     # Matrix of alpha_v values (n_iter × n_var)
                        alpha_u,      # scalar
                        Omega_u,      # A matrix associated with alpha_u
                        Omegas_v,     # List of length = n_iter, each element is list of length = n_var
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- nrow(alphas_v)
  n <- nrow(X)
  m <- ncol(X)
  GCV <- numeric(n_iter)

  if (all(alphas_v == 0)) {
    # S_alpha for v
    S_alpha_v0 <- Omega_v0 <- list()
    for (i in 1:n_var) {
      tds <- GridPoints_v[[i]]
      if(is.null(tds)){
        S_alpha_v0[[i]] <- Omega_v0[[i]] <- diag(ncol[,i])
      }else{
        getpenresult <- get.pen(td = tds,
                                alpha = as.numeric(alpha_v[i]),
                                type = smoothness_type)
        S_alpha_v0[[i]] <- getpenresult$S.alpha
        Omega_v0[[i]] <- getpenresult$Omega
        step <- step + 1
        setTxtProgressBar(pb, step)
      }
    }
    opt_alpha_Omega_v <- as.matrix(bdiag(lapply(1:n_var, function(j)
      as.numeric(rep(0, n_var)) * S_alpha_v0[[j]])))

    return(list(GCV_v = Inf,
                opt.alpha_v = rep(0, n_var),
                opt_s.alpha_v = S_alpha_v0,
                GCVdf_v = data.frame(alphas_v, rep(Inf, n_iter)),
                Omega_v = S_alpha_v0,
                opt_alpha_Omega_v = opt_alpha_Omega_v))
  } else {
    alpha_Omega_u <- alpha_u * Omega_u
    for (i in 1:n_iter) {
      S_alpha_v <- S_alphas_v[[i]]
      Omega_v <- Omegas_v[[i]]
      alpha_Omega_v <- as.matrix(bdiag(lapply(1:n_var, function(j)
        alphas_v[i, j] * Omega_v[[j]])))

      # Power Algorithm
      power_result <- power_algo(data = X,
                                 n_var = n_var,
                                 ncol = ncol,
                                 thresh = thresh,
                                 maxit = maxit,
                                 sparse_tuning_result_u =sparse_tuning_result_u,
                                 sparse_tuning_result_v =sparse_tuning_result_v,
                                 S_alpha_v = S_alpha_v,
                                 S_alpha_u = S_alphas_u,
                                 sparse_tuning_type = sparse_tuning_type,
                                 conditional = conditional,
                                 alpha_Omega_v = alpha_Omega_v,
                                 alpha_Omega_u = alpha_Omega_u)

      v <- power_result[[1]]
      u <- power_result[[2]]

      norm_u2 <- as.numeric(norm_vec(u)^2)
      Xu_proj <- as.vector(t(X) %*% u / norm_u2)
      R_u <- as.numeric(t(u) %*% Omega_u %*% u / norm_u2)

      GCV_alpha <- 0
      idx <- 0
      for (k in 1:n_var) {
        m_k <- as.integer(ncol[k])
        v_k <- v[(idx + 1):(idx + m_k)]
        xuk <- Xu_proj[(idx + 1):(idx + m_k)]
        S_k <- S_alpha_v[[k]]

        trace_k <- sum(diag(S_k))
        denom_k <- (1 - (1 / m_k) * (trace_k / (1 + alpha_u * R_u)))^2
        if(denom_k == 0){denom_k<-1}
        numer_k <- norm_vec(xuk - v_k)^2

        GCV_alpha <- GCV_alpha + numer_k / denom_k
        idx <- idx + m_k
      }

      GCV[i] <- GCV_alpha / m
    }

    opt.alpha <- alphas_v[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]
    Omega_v <- Omegas_v[[which.min(GCV)]]
    opt_alpha_Omega_v <- as.matrix(bdiag(lapply(1:n_var, function(j)
      as.numeric(opt.alpha) * Omega_v[[j]])))

    return(list(GCV_v = GCV,
                opt.alpha_v = opt.alpha,
                opt_s.alpha_v = opt_s.alpha,
                GCVdf_v = data.frame(alphas_v, GCV),
                Omega_v = Omega_v,
                opt_alpha_Omega_v = opt_alpha_Omega_v))
  }
}
