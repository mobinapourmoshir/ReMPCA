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
      GCV_alpha <- (1/m) * (norm_vec((diag(m) - S) %*% (t(X) %*% u))^2 / (norm_vec(u)^2 * (1 - (1/m) * sum(diag(S)))^2))

      #setTxtProgressBar(pb, i + (k - 1) / n_var)
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV))

    ########### Smoothness on u ###########
    if(two_way_smoothness != 0){
      GCV_u <- numeric(length(two_way_smoothness))
      for(j in 1:length(two_way_smoothness)){
        alpha_u <-  two_way_smoothness[j]
        S_u <-  S_alphas_u[[j]]
        GCV_alpha_u <- 0
        m <- ncol(t(X))

        power_result <- power_algo(data = t(X),
                                   S_alpha = S_u,
                                   sparse_tuning_result_u = sparse_tuning_result_u,
                                   sparse_tuning_result_v = sparse_tuning_result_v,
                                   sparse_tuning_type = sparse_tuning_type)
        u_u <- power_result[[2]]
        GCV_alpha_u <- (1/m) * (norm_vec((diag(m) - S_u) %*% (X %*% u_u))^2 / (norm_vec(u_u)^2 * (1 - (1/m) * sum(diag(S_u)))^2))
        GCV_u[j] <- GCV_alpha_u
      }
      opt.alpha_u <- two_way_smoothness[which.min(GCV_u), ]
      opt_s.alpha_u <- S_alphas_u[[which.min(GCV_u)]]

      result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV),
                     GCV_u = GCV_u, opt.alpha_u = opt.alpha_u, opt_s.alpha_u = opt_s.alpha_u, GCVdf_u = data.frame(two_way_smoothness, GCV_u))
    }

    return(result)
  }
}




############################### Calculating the optimal alpha for v using GCV ###############################
opt_alpha_u <- function(X,
                        n_var,
                        ncol,
                        S_alphas_v, # A list (length = n_var)
                        S_alphas_u, # A list (length = length(alphas))
                        alphas, # A vector
                        sparse_tuning_result_u,
                        sparse_tuning_result_v,
                        sparse_tuning_type) {

  n_iter <- nrow(alphas)  # Update to get the number of iterations
  #pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  #                     max = n_iter, # Maximum value of the progress bar
  #                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
  #                     width = 50,   # Progress bar width. Defaults to getOption("width")
  #                     char = "=")   # Character used to create the bar

  n <- nrow(X)
  GCV <- numeric(n_iter)

  if (all(alphas == 0)) {
    #close(pb)
    return(list(GCV = Inf, opt.alpha = 0, opt_s.alpha = diag(n), GCVdf = data.frame(alphas, rep(Inf, n_iter))))
  } else {
    for (i in 1:n_iter) {
      S <- S_alphas_u[[i]]
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
      GCV_alpha <- (1/m) * (norm_vec((diag(m) - S) %*% (t(X) %*% u))^2 / (norm_vec(u)^2 * (1 - (1/m) * sum(diag(S)))^2))

      #setTxtProgressBar(pb, i + (k - 1) / n_var)
      GCV[i] <- GCV_alpha
    }

    opt.alpha <- alphas[which.min(GCV), ]
    opt_s.alpha <- S_alphas_v[[which.min(GCV)]]
    #close(pb)
    result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV))

    ########### Smoothness on u ###########
    if(two_way_smoothness != 0){
      GCV_u <- numeric(length(two_way_smoothness))
      for(j in 1:length(two_way_smoothness)){
        alpha_u <-  two_way_smoothness[j]
        S_u <-  S_alphas_u[[j]]
        GCV_alpha_u <- 0
        m <- ncol(t(X))

        power_result <- power_algo(data = t(X),
                                   S_alpha = S_u,
                                   sparse_tuning_result_u = sparse_tuning_result_u,
                                   sparse_tuning_result_v = sparse_tuning_result_v,
                                   sparse_tuning_type = sparse_tuning_type)
        u_u <- power_result[[2]]
        GCV_alpha_u <- (1/m) * (norm_vec((diag(m) - S_u) %*% (X %*% u_u))^2 / (norm_vec(u_u)^2 * (1 - (1/m) * sum(diag(S_u)))^2))
        GCV_u[j] <- GCV_alpha_u
      }
      opt.alpha_u <- two_way_smoothness[which.min(GCV_u), ]
      opt_s.alpha_u <- S_alphas_u[[which.min(GCV_u)]]

      result <- list(GCV = GCV, opt.alpha = opt.alpha, opt_s.alpha = opt_s.alpha, GCVdf = data.frame(alphas, GCV),
                     GCV_u = GCV_u, opt.alpha_u = opt.alpha_u, opt_s.alpha_u = opt_s.alpha_u, GCVdf_u = data.frame(two_way_smoothness, GCV_u))
    }

    return(result)
  }
}

