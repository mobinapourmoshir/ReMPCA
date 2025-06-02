#' Plot GCV Scores for u Components
#'
#' This function visualizes the GCV scores for the left singular vectors (u)
#' across different smoothing parameter values (alpha_u), typically obtained
#' during regularization tuning.
#'
#' @param ReMPCA_obj An object returned from ReMPCA().
#' @param ... Additional plotting options.
#'
#' @return A plot object showing GCV scores across alpha_u values.
#' @export

plot_gcv_u <- function(ReMPCA_obj,...) {
  GCVlist <- ReMPCA_obj$GCVResultsU
  smooth_result_u <- ReMPCA_obj$OptimalAlphaU
  n_pc <- length(GCVlist)

  # Setup plotting window
  par(mfrow = c(1, n_pc), mar = c(4, 4, 2, 1))

  for (i in seq_len(n_pc)) {
    df <- GCVlist[[i]]
    alpha_vals <- df$alphas_u
    gcv_vals <- df$GCV

    # Plot GCV curve
    matplot(alpha_vals, gcv_vals, type = "l", log = "x",
            lty = 1, col = "black",
            xlab = expression(alpha[u]),
            ylab = "GCV",
            main = paste("Component", i))

    # Highlight optimal alpha
    opt_alpha <- smooth_result_u[[i]]
    opt_gcv <- gcv_vals[which.min(abs(alpha_vals - opt_alpha))]

    points(opt_alpha, opt_gcv, pch = 19, col = "red", cex = 1.5)
  }
}


#' Plot GCV Results for Functional Loadings (`v`) in ReMPCA
#'
#' This function visualizes the Generalized Cross-Validation (GCV) curves or surfaces
#' for the functional loading vectors (`v`) across components in a ReMPCA object.
#' It determines the number of functional variables being tuned and chooses the
#' appropriate plot type:
#'
#' - **Line plot** if only one variable is tuned
#' - **Filled contour plot** if exactly two variables are tuned
#' - **Error** if more than two variables are tuned simultaneously (not supported)
#'
#' The optimal smoothing parameter (selected \eqn{\alpha_v}) is highlighted in each plot.
#'
#' @param ReMPCA_obj A list-like ReMPCA object, expected to contain:
#'   \describe{
#'     \item{\code{GCVResultsV}}{A list of data frames, one per principal component (PC),
#'     each containing columns \code{Var1}, \code{Var2}, ..., and a \code{GCV} column.}
#'     \item{\code{OptimalAlphaV}}{A list of length equal to the number of PCs; each element
#'     is a named numeric vector with optimal smoothing values for the corresponding variables.}
#'   }
#'
#' @return Produces plots side by side in a loop, one for each component, and highlights the optimal \eqn{\alpha_v}.
#' Does not return a value.
#'
#' @details
#' If all variables for a PC have zero-valued tuning parameters, no plot is drawn.
#' @export
plot_gcv_v <- function(ReMPCA_obj) {
  GCVlist <- ReMPCA_obj$GCVResultsV
  smooth_result_v <- ReMPCA_obj$OptimalAlphaV
  n_pc <- length(GCVlist)

  layout(matrix(1:n_pc, nrow = 1))  # Arrange plots side by side

  for (i in seq_len(n_pc)) {
    df <- GCVlist[[i]]
    var_cols <- grep("^Var", names(df), value = TRUE)

    # Identify non-zero variance dimensions
    non_zero_vars <- var_cols[sapply(var_cols,
                                     function(vname) any(df[[vname]] != 0))]
    num_nonzero <- length(non_zero_vars)

    if (num_nonzero == 0) {
      message(sprintf("PC %d: All variables zero, skipping plot.", i))
      next
    }

    main_title <- paste("Component", i)

    if (num_nonzero == 1) {
      x_vals <- df[[non_zero_vars[1]]]
      y_vals <- df$GCV

      plot(x_vals, y_vals, type = "l", log = "x", lty = 1, col = "black",
           xlab = non_zero_vars[1], ylab = "GCV", main = main_title)

      # Highlight optimal alpha
      opt_alpha <- smooth_result_v[[i]][[which(var_cols == non_zero_vars[1])]]
      opt_gcv <- y_vals[which.min(abs(x_vals - opt_alpha))]

      points(opt_alpha, opt_gcv, pch = 19, col = "red", cex = 1.5)

    } else if (num_nonzero == 2) {
      x_var <- non_zero_vars[1]
      y_var <- non_zero_vars[2]

      # Create grid and matrix of GCV values
      x_vals <- log10(unique(df[[x_var]]))
      y_vals <- log10(unique(df[[y_var]]))
      z_mat <- matrix(df$GCV, nrow = length(x_vals),
                      ncol = length(y_vals), byrow = TRUE)

      opt_x <- log10(smooth_result_v[[i]][[which(var_cols == x_var)]])
      opt_y <- log10(smooth_result_v[[i]][[which(var_cols == y_var)]])

      filled.contour(x = x_vals, y = y_vals, z = z_mat,
                     main = main_title,
                     plot.axes = {
                       box()
                       # Add contour lines
                       contour(x_vals, y_vals, z_mat, add = TRUE, drawlabels = TRUE, col = "black")
                       # Add axis labels
                       mtext(x_var, side = 1, line = 1, cex = 1.1)
                       mtext(y_var, side = 2, line = 1, cex = 1.1)
                       # Highlight optimal (alpha1, alpha2)
                       points(opt_x, opt_y, pch = 19, col = "black", cex = 1.5)
                     })
    } else {
      stop(paste("GCV plotting for v is only supported for up to 2 functional variables.",
                 "Component", i, "has", num_nonzero, "non-zero variables."))
    }
  }
}

#' Plot CV Error with 1-SE Rule for u
#'
#' For each PC component, plots CV error vs. sparsity tuning parameter (`gamma`)
#' and highlights the optimal value selected via the 1-SE rule.
#'
#' @param ReMPCA_obj Output list from the ReMPCA routine.
#'
#' @return Side-by-side base R plots for each component showing the CV error and 1-SE threshold.
#' @export
plot_cv_u <- function(ReMPCA_obj) {
  CVlist <- ReMPCA_obj$CVResultsU
  opt_gamma <- ReMPCA_obj$OptimalGammaU
  n_pc <- length(CVlist)

  par(mfrow = c(1, n_pc), mar = c(4, 4, 2, 1))

  for (i in seq_len(n_pc)) {
    df <- CVlist[[i]]
    gamma_vals <- df$sparse_tuning_result_u
    cv_means <- df$CV_errors
    cv_ses <- df$SE_errors
    threshold <- min(cv_means) + cv_ses[which.min(cv_means)]

    plot(gamma_vals, cv_means, type = "b", pch = 19,
         xlab = bquote(gamma[u]), ylab = "CV Error",
         main = paste("Component", i))

    abline(h = threshold, col = "red", lty = 2)

    points(opt_gamma[[i]],
           cv_means[which(gamma_vals == opt_gamma[[i]])],
           pch = 19, col = "red", cex = 1.5)
  }
}


#' Plot CV Scores for v Components
#'
#' This function visualizes the cross-validation (CV) errors for each functional variable
#' across a range of sparsity parameters (`gamma_v`), for each principal component (PC).
#' It is designed for use with hybrid PCA models where multiple variables are tuned.
#'
#' @param ReMPCA_obj A list-like ReMPCA object, expected to contain:
#'   \describe{
#'     \item{\code{CVResultsV}}{A nested list of data frames, one per principal component (PC),
#'     each containing a list of data frames (one per variable) with columns:
#'     \code{gamma_Xi}, \code{cv_means}, and \code{cv_ses}.}
#'     \item{\code{OptimalGammaV}}{A list of numeric vectors with optimal gamma values
#'     for each variable in each PC.}
#'   }
#'
#' @return A grid of CV plots for variables and PCs. For each variable-PC pair with more than one tuning value,
#'   the function plots CV error vs. gamma, highlights the selected gamma (in blue), and the 1-SE threshold (red line).
#' @export
plot_cv_v <- function(ReMPCA_obj) {
  CV_v <- ReMPCA_obj$CVResultsV
  OptimalGammaV <- ReMPCA_obj$OptimalGammaV
  n_pc <- length(CV_v)
  n_var <- length(CV_v[[1]])

  # Determine which variable-PC pairs will be plotted
  to_plot <- matrix(FALSE, nrow = n_var, ncol = n_pc)
  for (i in seq_len(n_pc)) {
    for (j in seq_len(n_var)) {
      if (nrow(CV_v[[i]][[j]]) > 1) {
        to_plot[j, i] <- TRUE
      }
    }
  }

  # Count how many variables have at least one valid plot
  vars_with_plot <- which(rowSums(to_plot) > 0)
  n_plot_vars <- length(vars_with_plot)

  # Set layout: rows = # of variables with plots, cols = # of PCs
  par(mfrow = c(n_plot_vars, n_pc), mar = c(4, 4, 2, 1))

  for (j in vars_with_plot) {
    for (i in seq_len(n_pc)) {
      if (!to_plot[j, i]) {
        plot.new()
        next
      }

      df <- CV_v[[i]][[j]]
      gammas <- df$gamma_Xi
      cv_means <- df$cv_means
      cv_ses <- df$cv_ses
      opt_gamma <- OptimalGammaV[[i]][j]

      # 1-SE rule threshold
      j_min <- which.min(cv_means)
      threshold <- cv_means[j_min] + cv_ses[j_min]

      # Plot
      plot(gammas, cv_means, type = "b", pch = 19, col = "black",
           xlab = bquote(gamma[.(j)]),
           ylab = "CV Scores",
           main = paste(ordinal(j),"Functional Variable", " - PC", i))

      abline(h = threshold, lty = 2, col = "red")
      points(opt_gamma, cv_means[which.min(abs(gammas - opt_gamma))],
             pch = 19, col = "red", cex = 1.5)
    }
  }
}

#' Plot Principal Component Scores from ReMPCA Object
#'
#' This function visualizes the principal component scores from a ReMPCA object.
#' Each principal component (i.e., each column in `PCScores`) is plotted in a
#' separate dot plot, displaying the score values across observations.
#'
#' @param ReMPCA_obj A list-like ReMPCA object, expected to contain a component
#'        named `PCScores` which is a matrix or data frame where each column
#'        corresponds to the scores of one principal component.
#'
#' @return No return value. The function produces a series of dot plots, one for
#'         each principal component.
#' @export

plot_pc_scores <- function(ReMPCA_obj) {
  scores <- ReMPCA_obj$PCScores

  if (!is.data.frame(scores) && !is.matrix(scores)) {
    stop("'ReMPCA_obj$PCScores' must be a data frame or matrix.")
  }

  n_pc <- ncol(scores)
  par(mfrow = c(1, n_pc), mar = c(4, 4, 2, 1))

  for (i in seq_len(n_pc)) {
    plot(scores[, i], pch = 19, cex = 0.6, col = "black",
         xlab = "Observation", ylab = "Score",
         main = paste("Component", i))
    abline(h = 0, lty = 2, col = "gray")
  }
}

#' Plot Principal Component Functions
#'
#' This function visualizes the estimated principal component functions (v) from a
#' ReMPCA object. Each component is plotted across its variables using either a
#' line plot (for functional variables) or a dot plot (for regular variables).
#'
#' @param ReMPCA_obj A ReMPCA object containing:
#'   \describe{
#'     \item{\code{PCFunctions}}{A list of length equal to number of PCs. Each element is a list of PC functions for each variable.}
#'     \item{\code{variable_types}}{A character vector indicating type of each variable: either `"hd"` (functional) or `"rd"` (regular).}
#'   }
#'
#' @details The function arranges plots in a matrix layout with rows corresponding
#' to variables and columns to principal components. A light gray horizontal line
#' at 0 is added for reference unless the minimum value in the plot is ≥ 5.
#'
#' @return No return value. This function is called for its side effect of plotting.
#'
#' @export
plot_pc_functions <- function(ReMPCA_obj) {
  PCFunctions <- ReMPCA_obj$PCFunctions
  variable_types <- ReMPCA_obj$variable_types

  n_pc <- length(PCFunctions)
  n_var <- length(variable_types)

  par(mfrow = c(n_var, n_pc), mar = c(3, 3, 2, 1))

  for (i in seq_len(n_pc)) {
    pc_funcs <- PCFunctions[[i]]
    for (j in seq_len(n_var)) {
      func <- pc_funcs[[j]]
      y_min <- min(func, na.rm = TRUE)
      y_max <- max(func, na.rm = TRUE)

      plot(func, type = "n", main = paste("PC", i, "- Var", j),
           xlab = "", ylab = "", ylim = c(y_min, y_max))

      # Add gray zero line if min < 5
      if (y_min < 5) abline(h = 0, col = "gray80", lty = 2)

      if (variable_types[j] == "rd") {
        points(func, pch = 16, col = "black")
      } else {
        lines(func, col = "black")
      }
    }
  }
}
