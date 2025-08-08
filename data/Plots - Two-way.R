# library(ggplot2)
# library(gridExtra)
#
# # Normalize a vector
# norm_vec <- function(x) sqrt(sum(x^2))
#
# # Helper to set y-axis limits based on both true and estimated values
# get_ylim <- function(true1, est1, true2, est2) {
#   y_combined <- c(true1, est1, true2, est2)
#   range(y_combined, na.rm = TRUE)
# }
#
# # Plot function with shared y-limits
# plot_function <- function(true, est, title, show_x = FALSE, show_y = FALSE, ylim = NULL) {
#   df <- data.frame(x = seq_along(true), True = true, Estimated = est)
#   p <- ggplot(df, aes(x)) +
#     geom_line(aes(y = True), color = "black", linetype = "dashed") +
#     geom_line(aes(y = Estimated), color = "red", alpha = 0.7) +
#     labs(title = title) +
#     coord_cartesian(ylim = ylim) +
#     theme_minimal(base_size = 10) +
#     theme(
#       plot.title = element_text(size = 9, hjust = 0.5),
#       axis.title = element_blank(),
#       axis.text.x = if (show_x) element_text(size = 7) else element_blank(),
#       axis.text.y = if (show_y) element_text(size = 7) else element_blank(),
#       axis.ticks.x = if (show_x) element_line() else element_blank(),
#       axis.ticks.y = if (show_y) element_line() else element_blank()
#     )
#   return(p)
# }
#
# # Panel constructor: 3 rows × 2 columns
# make_panel <- function(u1, u2, v11, v12, v21, v22,
#                        u1_est, u2_est, v11_est, v12_est, v21_est, v22_est,
#                        panel_title) {
#   # Set shared y-axis limits for each row
#   ylim_u  <- get_ylim(u1, u1_est, u2, u2_est)
#   ylim_v1 <- get_ylim(v11, v11_est, v12, v12_est)
#   ylim_v2 <- get_ylim(v21, v21_est, v22, v22_est)
#
#   # Row 1: u1, u2
#   p1 <- plot_function(u1, u1_est, "u1", show_y = TRUE, ylim = ylim_u)
#   p2 <- plot_function(u2, u2_est, "u2", ylim = ylim_u)
#
#   # Row 2: v11, v12
#   p3 <- plot_function(v11, v11_est, "v11", show_y = TRUE, ylim = ylim_v1)
#   p4 <- plot_function(v12, v12_est, "v12", ylim = ylim_v1)
#
#   # Row 3: v21, v22
#   p5 <- plot_function(v21, v21_est, "v21", show_x = TRUE, show_y = TRUE, ylim = ylim_v2)
#   p6 <- plot_function(v22, v22_est, "v22", show_x = TRUE, ylim = ylim_v2)
#
#   gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6,
#                           nrow = 3, ncol = 2,
#                           top = grid::textGrob(panel_title, gp = grid::gpar(fontsize = 14, fontface = "bold")))
# }
#
# library(gridExtra)
#
# with(result1, {
#   # Panel 1
#   p1 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_ss_uv, u2_est_ss_uv, v11_est_ss_uv, v12_est_ss_uv, v21_est_ss_uv, v22_est_ss_uv,
#                "Two-way Sparsity and Smoothness")
#   )
#   ggsave("Two-way_Sparsity_and_Smoothness.pdf", p1, width = 6, height = 6)
#
#   # Panel 2
#   p2 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_SVD, u2_est_SVD, v11_est_SVD, v12_est_SVD, v21_est_SVD, v22_est_SVD,
#                "SVD")
#   )
#   ggsave("SVD.pdf", p2, width = 6, height = 6)
#
#   # Panel 3
#   p3 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_ss_u, u2_est_ss_u, v11_est_ss_u, v12_est_ss_u, v21_est_ss_u, v22_est_ss_u,
#                "Smoothness & Sparsity on u")
#   )
#   ggsave("Smoothness_and_Sparsity_on_u.pdf", p3, width = 6, height = 6)
#
#   # Panel 4
#   p4 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_sm_uv, u2_est_sm_uv, v11_est_sm_uv, v12_est_sm_uv, v21_est_sm_uv, v22_est_sm_uv,
#                "Two-Way Smoothness")
#   )
#   ggsave("Two-Way_Smoothness.pdf", p4, width = 6, height = 6)
#
#   # Panel 5
#   p5 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_sp_uv, u2_est_sp_uv, v11_est_sp_uv, v12_est_sp_uv, v21_est_sp_uv, v22_est_sp_uv,
#                "Two-Way Sparsity")
#   )
#   ggsave("Two-Way_Sparsity.pdf", p5, width = 6, height = 6)
#
#   # Panel 6
#   p6 <- arrangeGrob(
#     make_panel(u1 / norm_vec(u1), u2 / norm_vec(u2), v11, v12, v21, v22,
#                u1_est_ss_v, u2_est_ss_v, v11_est_ss_v, v12_est_ss_v, v21_est_ss_v, v22_est_ss_v,
#                "Smoothness & Sparsity on v")
#   )
#   ggsave("Smoothness_and_Sparsity_on_v.pdf", p6, width = 6, height = 6)
# })
#
#
# ############### Plot and table of ISE, R_ISE
# # Combine all ResultsTabels into one big data frame
# all_results <- do.call(rbind, lapply(results_list, function(x) x$ResultsTabel))
#
# # Compute mean ISE and mean R_ISE for each combination of parameter and method
# colnames(all_results)[colnames(all_results) == "R ISE"] <- "R_ISE"
#
#
# library(tidyr)
# library(dplyr)
#
# # First compute summary if not already done
# summary_df <- all_results %>%
#   group_by(param, method) %>%
#   summarise(
#     mean_ISE = mean(ISE, na.rm = TRUE),
#     mean_R_ISE = mean(R_ISE, na.rm = TRUE),
#     .groups = "drop"
#   )
#
# # Wide table for ISE
# ise_table <- summary_df %>%
#   select(method, param, mean_ISE) %>%
#   pivot_wider(names_from = param, values_from = mean_ISE)
#
# View(ise_table)
#
# # Wide table for R_ISE
# r_ise_table <- summary_df %>%
#   select(method, param, mean_R_ISE) %>%
#   pivot_wider(names_from = param, values_from = mean_R_ISE)
#
# View(r_ise_table)
#
# # LateX code
# library(knitr)
# ise_table_formatted <- ise_table %>%
#   mutate(across(-method, ~ sprintf("%.5f", .)))
#
# kable(ise_table_formatted, format = "latex", booktabs = TRUE,
#       caption = "Mean ISE for each method and parameter", label = "ISE")
#
# kable(ise_table, format = "latex", digits = 4, booktabs = TRUE,
#       caption = "Mean ISE for each method and parameter")
#
# kable(r_ise_table[-3,], format = "latex", digits = 4, booktabs = TRUE,
#       caption = "Mean Relative ISE for each method and parameter")
#
#
# library(dplyr)
# library(tidyr)
# library(knitr)
#
# # Reshape and compute v1 and v2
# ise_table <- summary_df %>%
#   select(method, param, mean_ISE) %>%
#   pivot_wider(names_from = param, values_from = mean_ISE) %>%
#   mutate(
#     v1 = (v11 + v21) / 2,
#     v2 = (v12 + v22) / 2
#   ) %>%
#   select(method, u1, u2, v1, v2) %>%
#   arrange(desc(u1))  # Sort in decreasing order of u1
#
# # Format for LaTeX with 5 decimal places
# ise_table_formatted <- ise_table %>%
#   mutate(across(-method, ~ sprintf("%.5f", .)))
#
# # Create LaTeX table
# kable(ise_table_formatted, format = "latex", booktabs = TRUE,
#       caption = "Mean ISE for each method and parameter",
#       label = "table: ISE")
#
#
#
# r_ise_table <- summary_df %>%
#   select(method, param, mean_R_ISE) %>%
#   pivot_wider(names_from = param, values_from = mean_R_ISE) %>%
#   mutate(
#     v1 = (v11 + v21) / 2,
#     v2 = (v12 + v22) / 2
#   ) %>%
#   select(method, u1, u2, v1, v2) %>%
#   arrange(desc(u1))
#
# # Format and show LaTeX table, excluding row 3
# kable(r_ise_table[-6, ], format = "latex", digits = 4, booktabs = TRUE,
#       caption = "Mean Relative ISE for each method and parameter",
#       label = "table: R")
#
#
# ############## Box Plot of ISE
# library(dplyr)
# library(tidyr)
# library(ggplot2)
#
# # 1) Bind all simulations + compute v1, v2 per sim/method
# all_results <- dplyr::bind_rows(
#   lapply(seq_along(results_list), function(i)
#     dplyr::mutate(results_list[[i]]$ResultsTabel, sim = i))
# )
#
# wide <- all_results %>%
#   filter(param %in% c("v11","v12","v21","v22")) %>%
#   select(sim, method, param, ISE) %>%
#   pivot_wider(names_from = param, values_from = ISE)
#
# long_v <- wide %>%
#   mutate(v1 = (v11 + v21)/2,
#          v2 = (v12 + v22)/2) %>%
#   select(sim, method, v1, v2) %>%
#   pivot_longer(cols = c(v1, v2), names_to = "component", values_to = "ISE")
#
# # Desired method order (edit strings to match yours exactly)
# method_levels <- c("SVD",
#                    "Smooth & Sparse u",
#                    "Smooth & Sparse v",
#                    "Two-way Smoothness",
#                    "Two-way Sparsity",
#                    "Smooth & Sparse u & v")
#
# long_v$method    <- factor(long_v$method, levels = method_levels)
# long_v$component <- factor(long_v$component, levels = c("v1","v2"))
#
# # 2) Boxplot: two boxes per method (v1, v2)
# ggplot(long_v, aes(x = method, y = ISE, fill = component)) +
#   geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier_size = 0.8) +
#   labs(x = "Method", y = "ISE", fill = "Component") +
#   theme_minimal(base_size = 12) +
#   theme(axis.text.x = element_text(angle = 25, hjust = 1))
#
#
#
