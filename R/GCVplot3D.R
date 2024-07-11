# library(dplyr)
# #library(akima)
#
# # GCV 3D Plot with increased jitter
# GCV_plot <- function(GCVdf, jitter_amount = 1e-6) {
#   names(GCVdf)[1] <- "alpha1"
#   names(GCVdf)[2] <- "alpha2"
#
#   # Add jitter to alpha1 and alpha2
#   alpha1_jittered <- GCVdf$alpha1 + runif(nrow(GCVdf), -jitter_amount, jitter_amount)
#   alpha2_jittered <- GCVdf$alpha2 + runif(nrow(GCVdf), -jitter_amount, jitter_amount)
#
#   interp_data <- with(GCVdf, interp(x = alpha1_jittered, y = alpha2_jittered, z = GCV))
#
#   fig <- plot_ly(z = ~interp_data$z, x = ~interp_data$x, y = ~interp_data$y, type = "surface")
#   fig <- fig %>% layout(scene = list(
#     xaxis = list(title = 'alpha 1'),
#     yaxis = list(title = 'alpha 2'),
#     zaxis = list(title = 'GCV')
#   ))
#
#   return(fig)
# }
