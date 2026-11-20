
required <- c("ggplot2", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required R packages: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# Locate this script whether called using Rscript or source().
script_directory <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f) && file.exists(f)) {
      return(dirname(normalizePath(f)))
    }
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    f <- sub("^--file=", "", file_arg[[1]])
    if (file.exists(f)) return(dirname(normalizePath(f)))
  }
  pkg_dir <- system.file("reproducibility", "figure1", package = "ReMPCA")
  if (nzchar(pkg_dir)) return(pkg_dir)
  getwd()
}

# All output goes beside the R file, regardless of the working directory.
script_dir <- script_directory()

# Read archived plotting inputs. For offline reproducibility, keep both CSV.gz
# files beside fig1_plot.R. Missing inputs are downloaded into R's temporary
# directory and deleted after reading, leaving only the PDF as an output.
files <- c("fig1_heatmaps_from_pdf.csv.gz", "fig1_vectors_from_pdf.csv.gz")
upstream <- paste0(
  "https://raw.githubusercontent.com/mmadoliat/ReMPCA/master/",
  "inst/reproducibility/figure1/"
)
read_plot_data <- function(name) {
  path <- file.path(script_dir, name)
  if (file.exists(path)) {
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }

  url <- paste0(upstream, name)
  tmp <- tempfile(pattern = "fig1_", fileext = ".csv.gz")
  on.exit(unlink(tmp), add = TRUE)
  message("Reading missing input from original repository: ", name)
  tryCatch({
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
    utils::read.csv(tmp, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop("Cannot read ", name, " from ", url, ": ", conditionMessage(e),
         call. = FALSE)
  })
}

H <- read_plot_data(files[[1]])
L <- read_plot_data(files[[2]])
stopifnot(all(c("image", "row", "col", "gray") %in% names(H)),
          all(c("vector", "method", "role", "index", "value_std") %in% names(L)),
          all(is.finite(H$gray)), all(H$gray >= 0 & H$gray <= 1),
          all(is.finite(L$value_std)))

img_lev <- c("True", "Noisy", "SVD", "Smooth", "Group lasso",
             "Smooth + group lasso")
img_lab <- c("True", "Noisy", "SVD", "Smooth", "Group lasso",
             "Smooth +\ngroup lasso")
meth_lev <- c("SVD", "Smooth", "Group lasso", "Smooth + group lasso")
meth_col <- c("SVD" = "#E69F00", "Smooth" = "#009E73",
              "Group lasso" = "#CC79A7", "Smooth + group lasso" = "#0072B2")
truth_col <- "grey70"
SHOW_KEY <- TRUE
# Set REMPCA_FIG1_FONT=Arial to use Arial on systems where it is installed.
FONT <- Sys.getenv("REMPCA_FIG1_FONT", unset = "sans")
FS_TITLE <- 8.5
FS_SMALL <- 7.5
LW_TRUTH <- 1.6 / .pt
LW_EST <- 0.8 / .pt

base_theme <- theme_void(base_family = FONT, base_size = FS_TITLE) +
  theme(strip.text = element_text(size = FS_TITLE, colour = "black",
                                  margin = margin(0, 0, 2, 0)),
        plot.margin = margin(1, 2, 1, 2))

# Top: six reconstructed chessboard heatmaps.
H$panel <- factor(img_lab[match(H$image, img_lev)], levels = img_lab)
H$fill <- grDevices::rgb(H$gray, H$gray, H$gray)
stopifnot(!anyNA(H$panel))
# Each statement starts with the plot object, avoiding a dangling '+' at the
# R console when lines are executed separately.
p_img <- ggplot(H, aes(x = col, y = row, fill = fill))
p_img <- p_img + geom_raster(interpolate = FALSE)
p_img <- p_img + scale_fill_identity()
p_img <- p_img + scale_x_continuous(expand = c(0, 0))
p_img <- p_img + scale_y_reverse(expand = c(0, 0))
p_img <- p_img + facet_wrap(~panel, nrow = 1)
p_img <- p_img + base_theme
p_img <- p_img + theme(
  aspect.ratio = 1,
  strip.text = element_text(
    size = FS_TITLE, vjust = 0, lineheight = 0.9,
    margin = margin(0, 0, 2.5, 0)
  ),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
  panel.spacing.x = grid::unit(0.075, "in"),
  plot.margin = margin(-4, 2, 7, 2)
)

# Bottom: u (row vector) and v (column vector), four methods each.
# The plotting values are standardized *within each original PDF panel* so
# the two truth levels are 0 and 1; this is not an estimate from ReMPCA.
L$method <- factor(L$method, levels = meth_lev)
L$role <- factor(L$role, levels = c("truth", "estimate"))
stopifnot(!anyNA(L$method), !anyNA(L$role))
row_lab <- c(u = 'bolditalic(u)~"(rows)"',
             v = 'bolditalic(v)~"(columns)"')

vec_plot <- function(vec, show_titles, breaks) {
  d <- L[L$vector == vec, , drop = FALSE]
  if (!nrow(d)) stop("No vector data found for ", vec)
  d$rowlab <- row_lab[[vec]]
  n <- max(d$index)

  p <- ggplot(d, aes(x = index, y = value_std))
  p <- p + geom_line(
    data = d[d$role == "truth", , drop = FALSE],
    colour = truth_col, linewidth = LW_TRUTH,
    lineend = "round", linejoin = "round"
  )
  p <- p + geom_line(
    data = d[d$role == "estimate", , drop = FALSE],
    mapping = aes(colour = method), linewidth = LW_EST,
    lineend = "round", linejoin = "round"
  )
  p <- p + scale_colour_manual(values = meth_col, guide = "none")
  p <- p + scale_x_continuous(
    breaks = breaks, limits = c(1, n), expand = expansion(mult = 0.03)
  )
  p <- p + scale_y_continuous(expand = expansion(mult = 0.06))
  p <- p + facet_grid(
    rowlab ~ method, switch = "y",
    labeller = labeller(rowlab = label_parsed, method = label_value)
  )
  p <- p + base_theme
  p <- p + theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.x = grid::unit(0.13, "in"),
    strip.placement = "outside",
    strip.text.x = if (show_titles) {
      element_text(size = FS_TITLE, margin = margin(0, 0, 2.5, 0))
    } else {
      element_blank()
    },
    strip.text.y.left = element_text(
      size = FS_TITLE, angle = 90, margin = margin(0, 3, 0, 0)
    ),
    axis.text.x = element_text(
      size = FS_SMALL, colour = "grey20", margin = margin(1.2, 0, 0, 0)
    ),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.25),
    axis.ticks.length.x = grid::unit(1.6, "pt")
  )
  p
}

p_u <- vec_plot("u", TRUE, c(1, 20, 40, 60, 80))
p_v <- vec_plot("v", FALSE, c(1, 20, 40, 60))

# Key below the vector panels.
p_key <- ggplot()
p_key <- p_key + annotate(
  "segment", x = 0.00, xend = 0.045, y = 0, yend = 0,
  colour = truth_col, linewidth = LW_TRUTH, lineend = "round"
)
p_key <- p_key + annotate(
  "text", x = 0.055, y = 0, label = "True vector", hjust = 0,
  size = FS_SMALL / .pt, family = FONT
)
p_key <- p_key + annotate(
  "segment", x = 0.24 + 0.03 * (0:3),
  xend = 0.24 + 0.03 * (0:3) + 0.022,
  y = 0, yend = 0, colour = unname(meth_col),
  linewidth = LW_EST * 1.6, lineend = "round"
)
p_key <- p_key + annotate(
  "text", x = 0.365, y = 0,
  label = "Estimated vector (colour = method)", hjust = 0,
  size = FS_SMALL / .pt, family = FONT
)
p_key <- p_key + coord_cartesian(
  xlim = c(-0.02, 1), ylim = c(-1, 1),
  expand = FALSE, clip = "off"
)
p_key <- p_key + theme_void()
p_key <- p_key + theme(plot.margin = margin(3, 2, 0, 2))

plots <- list(p_img, p_u, p_v)
hts <- c(65.6, 1, 1)
hun <- c("pt", "null", "null")
if (SHOW_KEY) {
  plots <- c(plots, list(p_key))
  hts <- c(hts, 0.10)
  hun <- c(hun, "null")
}
fig <- patchwork::wrap_plots(plots, ncol = 1,
                             heights = grid::unit(hts, hun))

# Write only a publication-quality PDF in the same folder as fig1_plot.R.
W <- 6.3
Ht <- 3.75
out_pdf <- file.path(script_dir, "fig1_chessboard.pdf")
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_pdf, width = W, height = Ht, family = FONT)
} else {
  grDevices::pdf(out_pdf, width = W, height = Ht, family = "sans")
}
print(fig)
invisible(grDevices::dev.off())
message("Saved: ", normalizePath(out_pdf))
