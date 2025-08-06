library(ggplot2)
library(gridExtra)

# Normalize a vector
norm_vec <- function(x) sqrt(sum(x^2))

# Helper to set y-axis limits based on both true and estimated values
get_ylim <- function(true1, est1, true2, est2) {
  y_combined <- c(true1, est1, true2, est2)
  range(y_combined, na.rm = TRUE)
}

# Plot function with shared y-limits
plot_function <- function(true, est, title, show_x = FALSE, show_y = FALSE, ylim = NULL) {
  df <- data.frame(x = seq_along(true), True = true, Estimated = est)
  p <- ggplot(df, aes(x)) +
    geom_line(aes(y = True), color = "black", linetype = "dashed") +
    geom_line(aes(y = Estimated), color = "red", alpha = 0.7) +
    labs(title = title) +
    coord_cartesian(ylim = ylim) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.title = element_blank(),
      axis.text.x = if (show_x) element_text(size = 7) else element_blank(),
      axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank()
    )
  return(p)
}

# Panel constructor: 3 rows × 2 columns
make_panel <- function(u1, u2, v11, v12, v21, v22,
                       u1_est, u2_est, v11_est, v12_est, v21_est, v22_est,
                       panel_title) {
  # Set shared y-axis limits for each row
  ylim_u  <- get_ylim(u1, u1_est, u2, u2_est)
  ylim_v1 <- get_ylim(v11, v11_est, v12, v12_est)
  ylim_v2 <- get_ylim(v21, v21_est, v22, v22_est)

  # Row 1: u1, u2
  p1 <- plot_function(u1, u1_est, "u1", show_y = TRUE, ylim = ylim_u)
  p2 <- plot_function(u2, u2_est, "u2", ylim = ylim_u)

  # Row 2: v11, v12
  p3 <- plot_function(v11, v11_est, "v11", show_y = TRUE, ylim = ylim_v1)
  p4 <- plot_function(v12, v12_est, "v12", ylim = ylim_v1)

  # Row 3: v21, v22
  p5 <- plot_function(v21, v21_est, "v21", show_x = TRUE, show_y = TRUE, ylim = ylim_v2)
  p6 <- plot_function(v22, v22_est, "v22", show_x = TRUE, ylim = ylim_v2)

  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6,
                          nrow = 3, ncol = 2,
                          top = grid::textGrob(panel_title, gp = grid::gpar(fontsize = 14, fontface = "bold")))
}

library(gridExtra)

with(result1, {
  # Panel 1
  p1 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_ss_uv, u2_est_ss_uv, v11_est_ss_uv, v12_est_ss_uv, v21_est_ss_uv, v22_est_ss_uv,
               "Two-way Sparsity and Smoothness")
  )
  ggsave("Two-way_Sparsity_and_Smoothness.pdf", p1, width = 6, height = 6)

  # Panel 2
  p2 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_SVD, u2_est_SVD, v11_est_SVD, v12_est_SVD, v21_est_SVD, v22_est_SVD,
               "SVD")
  )
  ggsave("SVD.pdf", p2, width = 6, height = 6)

  # Panel 3
  p3 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_ss_u, u2_est_ss_u, v11_est_ss_u, v12_est_ss_u, v21_est_ss_u, v22_est_ss_u,
               "Smoothness & Sparsity on u")
  )
  ggsave("Smoothness_and_Sparsity_on_u.pdf", p3, width = 6, height = 6)

  # Panel 4
  p4 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_sm_uv, u2_est_sm_uv, v11_est_sm_uv, v12_est_sm_uv, v21_est_sm_uv, v22_est_sm_uv,
               "Two-Way Smoothness")
  )
  ggsave("Two-Way_Smoothness.pdf", p4, width = 6, height = 6)

  # Panel 5
  p5 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_sp_uv, u2_est_sp_uv, v11_est_sp_uv, v12_est_sp_uv, v21_est_sp_uv, v22_est_sp_uv,
               "Two-Way Sparsity")
  )
  ggsave("Two-Way_Sparsity.pdf", p5, width = 6, height = 6)

  # Panel 6
  p6 <- arrangeGrob(
    make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
               u1_est_ss_v, u2_est_ss_v, v11_est_ss_v, v12_est_ss_v, v21_est_ss_v, v22_est_ss_v,
               "Smoothness & Sparsity on v")
  )
  ggsave("Smoothness_and_Sparsity_on_v.pdf", p6, width = 6, height = 6)
})
