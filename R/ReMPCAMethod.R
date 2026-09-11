#' Plot GCV Scores for U
#'
#' This function visualizes the GCV scores for the left singular vectors (U)
#' across different regularization values.
#'
#' @param ReMPCA_obj Output list from ReMPCA.
#' @param ... Additional plotting parameters.
#' @return A base R plot.
#' @export
#' @importFrom graphics text arrows
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

    # If all GCV values are Inf or NA
    if (all(is.infinite(gcv_vals)) || all(is.na(gcv_vals))) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = paste("Component", i))
      text(1, 1, "GCV scores are Inf!", cex = 1.2)
      next
    }

    # Plot GCV curve
    matplot(alpha_vals, gcv_vals, type = "l", log = "x",
            lty = 1, col = "black",
            xlab = expression(alpha[u]),
            ylab = "GCV",
            main = paste("Component", i),
            ...)

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
#'
#' @param ReMPCA_obj A list-like ReMPCA object with GCVResultsV and OptimalAlphaV.
#' @param ... Additional plotting options.
#' @export
plot_gcv_v <- function(ReMPCA_obj, ...) {
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

      if (all(is.infinite(y_vals)) || all(is.na(y_vals))) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = main_title)
        text(1, 1, "GCV scores are Inf!", cex = 1.2)
        next
      }

      plot(x_vals, y_vals, type = "l", log = "x",
           lty = 1, col = "black",
           xlab = non_zero_vars[1],
           ylab = "GCV", main = main_title,
           ...)

      # Highlight optimal alpha
      opt_alpha <- smooth_result_v[[i]][[which(var_cols == non_zero_vars[1])]]
      opt_gcv <- y_vals[which.min(abs(x_vals - opt_alpha))]

      points(opt_alpha, opt_gcv, pch = 19, col = "red", cex = 1.5)

    } else if (num_nonzero == 2) {
      x_var <- non_zero_vars[1]
      y_var <- non_zero_vars[2]

      x_vals_raw <- df[[x_var]]
      y_vals_raw <- df[[y_var]]
      z_vals <- df$GCV

      if (all(is.infinite(z_vals)) || all(is.na(z_vals))) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = main_title)
        text(1, 1, "GCV scores are Inf!", cex = 1.2)
        next
      }

      x_vals <- log10(unique(x_vals_raw))
      y_vals <- log10(unique(y_vals_raw))
      z_mat <- matrix(z_vals, nrow = length(x_vals),
                      ncol = length(y_vals), byrow = TRUE)

      opt_x <- log10(smooth_result_v[[i]][[which(var_cols == x_var)]])
      opt_y <- log10(smooth_result_v[[i]][[which(var_cols == y_var)]])

      filled.contour(x = x_vals, y = y_vals, z = z_mat,
                     main = main_title,
                     plot.axes = {
                       box()
                       contour(x_vals, y_vals, z_mat, add = TRUE,
                               drawlabels = TRUE, col = "black")
                       mtext(x_var, side = 1, line = 1, cex = 1.1)
                       mtext(y_var, side = 2, line = 1, cex = 1.1)
                       points(opt_x, opt_y, pch = 19, col = "black",
                              cex = 1.5)
                     },
                     ...)
    } else {
      stop(paste("GCV plotting for v is only supported for up to 2 functional variables.",
                 "Component", i, "has", num_nonzero, "non-zero variables."))
    }
  }
}

#' Plot CV Error with 1-SE Rule for u
#'
#' For each principal component (PC), plots cross-validation (CV) error versus the
#' sparsity tuning parameter (`gamma`) for the u-direction. The red dashed line
#' indicates the 1-SE rule threshold, and the selected optimal gamma is highlighted.
#' Optionally, standard error bars (±1 SE) can be displayed.
#'
#' @param ReMPCA_obj Output list from the ReMPCA routine.
#' @param show_se Logical. If TRUE, adds standard error bars to each point. Default is TRUE.
#' @param ... Additional plotting options passed to `plot()`.
#' @return Generates base R plots side-by-side for each PC.
#' @export
plot_cv_u <- function(ReMPCA_obj, show_se = TRUE, ...) {
  CVlist <- ReMPCA_obj$CVResultsU
  opt_gamma <- ReMPCA_obj$OptimalGammaU
  n_pc <- length(CVlist)

  par(mfrow = c(1, n_pc), mar = c(4, 4, 2, 1))

  for (i in seq_len(n_pc)) {
    df <- CVlist[[i]]
    gamma_vals <- df$sparse_tuning_result_u
    cv_means <- df$CV_errors
    cv_ses <- df$SE_errors

    # Handle Inf or NA
    if (all(is.infinite(cv_means)) || all(is.na(cv_means))) {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = paste("Component", i))
      text(1, 1, "CV scores are Inf!", cex = 1.2)
      next
    }

    threshold <- min(cv_means) + cv_ses[which.min(cv_means)]

    # Compute local ylim for this component
    local_ymin <- min(cv_means - cv_ses, na.rm = TRUE)
    local_ymax <- max(cv_means + cv_ses, na.rm = TRUE)

    plot(gamma_vals, cv_means, type = "b", pch = 19,
         xlab = bquote(gamma[u]), ylab = "CV Error",
         main = paste("Component", i),
         ylim = c(local_ymin, local_ymax), ...)

    # Optional: show SE bars
    if (show_se) {
      arrows(gamma_vals, cv_means - cv_ses,
             gamma_vals, cv_means + cv_ses,
             angle = 90, code = 3, length = 0.05, col = "gray40")
    }

    abline(h = threshold, col = "red", lty = 2)

    points(opt_gamma[[i]],
           cv_means[which.min(abs(gamma_vals - opt_gamma[[i]]))],
           pch = 19, col = "red", cex = 1.5)
  }
}

#' Plot Cross-Validation Scores for ReMPCA v-direction
#'
#' This function creates a grid of plots showing cross-validation (CV) scores
#' for each principal component (PC) and functional variable based on different
#' values of the gamma penalty parameter in the v-direction. The red dashed line
#' indicates the 1-SE rule threshold, and the selected optimal gamma is highlighted.
#'
#' @param ReMPCA_obj A list object returned by ReMPCA containing `CVResultsV` and `OptimalGammaV`.
#' @param show_se Logical. If TRUE, adds standard error bars (±1 SE) to each point. Default is TRUE.
#' @param ... Additional arguments passed to the `plot()` function.
#' @export
plot_cv_v <- function(ReMPCA_obj, show_se = TRUE, ...) {
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

  vars_with_plot <- which(rowSums(to_plot) > 0)
  n_plot_vars <- length(vars_with_plot)

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

      if (all(is.infinite(cv_means)) || all(is.na(cv_means))) {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
             main = paste("Variable", j, "- PC", i))
        text(1, 1, "CV scores are Inf!", cex = 1.2)
        next
      }

      # 1-SE rule threshold
      j_min <- which.min(cv_means)
      threshold <- cv_means[j_min] + cv_ses[j_min]

      # Local ylim for this plot
      local_ymin <- min(cv_means - cv_ses, na.rm = TRUE)
      local_ymax <- max(cv_means + cv_ses, na.rm = TRUE)

      plot(gammas, cv_means, type = "b", pch = 19, col = "black",
           xlab = bquote(gamma[.(j)]), ylab = "CV Scores",
           main = paste("Variable", j, "- PC", i),
           ylim = c(local_ymin, local_ymax), ...)

      if (show_se) {

        keep <- is.finite(cv_ses) & cv_ses > 0

        if (any(keep)) {
          arrows(
            gammas[keep],
            cv_means[keep] - cv_ses[keep],
            gammas[keep],
            cv_means[keep] + cv_ses[keep],
            angle = 90,
            code = 3,
            length = 0.05,
            col = "gray40"
          )
        }
      }

      abline(h = threshold, lty = 2, col = "red")

      points(opt_gamma,
             cv_means[which.min(abs(gammas - opt_gamma))],
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
#' @param ... Additional plotting options.
#' @return No return value. The function produces a series of dot plots, one for
#'         each principal component.
#' @export

plot_pc_scores <- function(ReMPCA_obj, ...) {
  scores <- ReMPCA_obj$PCScores

  if (!is.data.frame(scores) && !is.matrix(scores)) {
    stop("'ReMPCA_obj$PCScores' must be a data frame or matrix.")
  }

  n_pc <- ncol(scores)
  par(mfrow = c(1, n_pc), mar = c(4, 4, 2, 1))

  for (i in seq_len(n_pc)) {
    plot(scores[, i], pch = 19, cex = 0.6, col = "black",
         xlab = "Observation", ylab = "Score",
         main = paste("Component", i), ...)
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
#' @param ... Additional plotting options.
#' @details The function arranges plots in a matrix layout with rows corresponding
#' to variables and columns to principal components. A light gray horizontal line
#' at 0 is added for reference unless the minimum value in the plot is ≥ 5.
#'
#' @return No return value. This function is called for its side effect of plotting.
#'
#' @export
plot_pc_functions <- function(ReMPCA_obj, ...) {
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
           xlab = "", ylab = "", ylim = c(y_min, y_max), ...)

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
