# Figure 3: True / Noisy / Fitted hybrid rank-one signals
# Reference-matched colors, 3D viewpoint, and panel layout
# Dependencies: base R only. Run: Rscript figure3_reference_style.R
# Output: figure3_reference_style.pdf (in the script directory)

script_directory <- function() {
  fr <- sys.frames()
  for (i in rev(seq_along(fr))) {
    f <- fr[[i]]$ofile
    if (!is.null(f) && length(f) == 1L && file.exists(f)) {
      return(dirname(normalizePath(f, winslash = "/")))
    }
  }
  args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(args)) {
    f <- sub("^--file=", "", args[1L])
    if (file.exists(f)) return(dirname(normalizePath(f, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

# ---------------------------------------------------------------------------
# 1. Global settings
# ---------------------------------------------------------------------------
FIG_WIDTH  <- 6.50
FIG_HEIGHT <- 6.90
N          <- 150L
BLOCK_SIZE <- c(40L, 20L, 30L)
NOISE_SD   <- 0.205

# Colors sampled/approximated from the provided image.
# Their lighter/darker faces are produced by the lighting function below.
SURFACE_COLORS <- c(True = "#050578", Noisy = "#CA971D", Fitted = "#C65A46")
DOT_COLORS     <- c(True = "#090944", Noisy = "#9B6C0B", Fitted = "#993F35")
FRAME_COLOR    <- "#303030"

# Position/width of columns, as fractions of the full figure width.
# Column 1 = u, columns 2 and 4 = 3D surfaces, column 3 = grayscale block.
COLS <- list(
  c(.040, .125),
  c(.180, .230),
  c(.450, .255),
  c(.716, .230)
)

# Row entries: [surface bottom, surface top, loading bottom, loading top]
# in units of overall figure height.  The taller surface regions reproduce
# the tall rectangular 3D frames seen in the supplied image.
GROUPS <- list(
  c(.666, .884, .889, .979),   # True
  c(.342, .560, .565, .660),   # Noisy
  c(.027, .247, .252, .347)    # Fitted
)

# This is a reference-matched *perspective-style* projection, not the original
# affine (parallel) projection.  The top and bottom of each 3D box have
# different slopes, giving the same view seen in the example.  Variables 1
# and 3 intentionally share this viewpoint.
CAMERA <- list(
  x0       = .105,
  u_dx     = .615,
  v_dx     = .195,
  bottom_y = .235,
  bottom_u = -.140,
  bottom_v =  .170,
  top_y    = .885,
  top_u    =  .110,
  top_v    = -.150
)

# ---------------------------------------------------------------------------
# 2. Generate the data (unchanged statistical construction)
# ---------------------------------------------------------------------------
s  <- seq(0, 1, length.out = N)
t1 <- seq(0, 1, length.out = BLOCK_SIZE[1L])
t3 <- seq(0, 1, length.out = BLOCK_SIZE[3L])

active <- pmin(1, pmax(0, (s - .58) / .42))
u_true <- ifelse(s > .58, sin(pi * active)^1.3, 0)
v1_true <- ifelse(t1 <= .49, 1.15 * sin(pi * pmin(1, t1 / .49))^2, 0) -
  .72 * exp(-((t1 - .52) / .075)^2)
v2_true <- numeric(BLOCK_SIZE[2L])
v2_true[c(2L, 6L, 10L, 14L, 18L)] <- c(1.2, 1.5, 1.1, 1.4, 1.6)
v3_true <- -1.08 * sin(pi * t3)
v_true <- c(v1_true, v2_true, v3_true)

indices <- list(
  seq_len(BLOCK_SIZE[1L]),
  BLOCK_SIZE[1L] + seq_len(BLOCK_SIZE[2L]),
  sum(BLOCK_SIZE[1:2]) + seq_len(BLOCK_SIZE[3L])
)
X_true <- tcrossprod(u_true, v_true)

# Deterministic pseudo-normal values for a reproducible illustrative figure.
normal_hash <- function(nr, nc) {
  ii <- matrix(rep(seq_len(nr), nc), nrow = nr, ncol = nc)
  jj <- matrix(rep(seq_len(nc), each = nr), nrow = nr, ncol = nc)
  a <- sin(ii * 127.1 + jj * 311.7) * 43758.5453
  b <- sin(ii * 269.5 + jj * 183.3) * 18731.3672
  a <- a - floor(a)
  b <- b - floor(b)
  sqrt(-2 * log(pmax(a, 1e-9))) * cos(2 * pi * b)
}

X_noisy <- X_true + NOISE_SD * normal_hash(N, sum(BLOCK_SIZE))
u_noisy <- u_true + .18 * normal_hash(N, 1L)[, 1L]
v_noisy <- v_true + .12 * normal_hash(1L, sum(BLOCK_SIZE))[1L, ]

# ---------------------------------------------------------------------------
# 3. Regularized rank-one fit (same method and tuning as original)
# ---------------------------------------------------------------------------
second_difference_smoother <- function(m, lambda) {
  I <- diag(m)
  D <- I[3:m, , drop = FALSE] - 2 * I[2:(m - 1L), , drop = FALSE] +
    I[1:(m - 2L), , drop = FALSE]
  solve(I + lambda * crossprod(D))
}

keep_largest <- function(z, k) {
  z <- as.numeric(z)
  k <- as.integer(min(length(z), max(0L, k)))
  if (k == 0L) return(z * 0)
  if (k >= length(z)) return(z)
  z[order(abs(z), decreasing = FALSE)[seq_len(length(z) - k)]] <- 0
  z
}

unit <- function(z) {
  z <- as.numeric(z)
  nrm <- sqrt(sum(z^2))
  if (!is.finite(nrm) || nrm < 1e-12) stop("Rank-one fit has zero norm")
  z / nrm
}

Su <- second_difference_smoother(N, 18)
Sv <- list(
  second_difference_smoother(BLOCK_SIZE[1L], 6),
  diag(BLOCK_SIZE[2L]),
  second_difference_smoother(BLOCK_SIZE[3L], 5)
)
v <- unit(svd(X_noisy, nu = 0, nv = 1)$v[, 1L])

for (iter in seq_len(350L)) {
  u <- unit(keep_largest(drop(Su %*% (X_noisy %*% v)), 61L))
  next_v <- numeric(length(v))
  for (j in seq_along(indices)) {
    ix <- indices[[j]]
    next_v[ix] <- keep_largest(
      drop(Sv[[j]] %*% crossprod(X_noisy[, ix, drop = FALSE], u)),
      c(23L, 5L, 30L)[j]
    )
  }
  next_v <- unit(next_v)
  if (sum(next_v * v) < 0) next_v <- -next_v
  err <- sum((next_v - v)^2)
  v <- next_v
  if (err < 1e-11) break
}

sigma <- as.numeric(crossprod(u, X_noisy %*% v))
if (sigma < 0) { v <- -v; sigma <- -sigma }
if (sum(u * u_true) < 0) { u <- -u; v <- -v }
u_fit <- u * sqrt(sigma)
v_fit <- v * sqrt(sigma)
X_fit <- tcrossprod(u_fit, v_fit)

rows <- list(
  list(name = "True", u = u_true, v = v_true, X = X_true),
  list(name = "Noisy", u = u_noisy, v = v_noisy, X = X_noisy),
  list(name = "Fitted", u = u_fit, v = v_fit, X = X_fit)
)

# ---------------------------------------------------------------------------
# 4. Drawing functions, using base R only
# ---------------------------------------------------------------------------
new_panel_factory <- function() {
  first <- TRUE
  function(x, y, w, h) {
    graphics::par(
      fig = c(x, x + w, y, y + h), mar = c(0, 0, 0, 0),
      oma = c(0, 0, 0, 0), xpd = FALSE, new = !first
    )
    first <<- FALSE
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1),
                          xaxs = "i", yaxs = "i")
  }
}

panel_title <- function(label, x = .5, y = .84, cex = 1.05) {
  graphics::text(x, y, label, cex = cex, font = 2)
}

plot_u <- function(d) {
  max_u <- max(vapply(rows, function(r) max(abs(r$u)), numeric(1))) * 1.10
  xx <- .36 + .39 * d$u / max_u
  yy <- seq(.055, .685, length.out = length(xx))
  graphics::lines(xx, yy, col = "#111111", lwd = .7)
  panel_title(paste(d$name, "u"), x = .45, y = .735, cex = 1.05)
}

plot_v <- function(d, j) {
  ix <- indices[[j]]
  values <- d$v[ix]
  a <- max(vapply(rows, function(r) max(abs(r$v[ix])), numeric(1))) * 1.06
  xx <- seq(.01, .99, length.out = length(values))
  yy <- .42 + .32 * values / a

  if (j == 2L) {
    graphics::segments(0, .42, 1, .42, col = "#B5B5B5", lwd = .35)
    graphics::points(xx, yy, pch = 16, cex = .42,
                     col = DOT_COLORS[[d$name]])
  } else {
    graphics::lines(xx, yy, col = "#141414", lwd = .68)
  }
  panel_title(paste0(d$name, " v", j), y = .85, cex = 1.06)
  # Put the small labels ABOVE the surface, as in the reference image.
  prefix <- if (d$name == "Fitted") "Reconstructed" else d$name
  graphics::text(.5, .12, paste(prefix, "Variable", j),
                 cex = .74, font = 2, col = "#1B1B1B")
}

plot_discrete <- function(d) {
  Z <- abs(d$X[, indices[[2L]], drop = FALSE])
  max_ref <- max(abs(X_true[, indices[[2L]], drop = FALSE]),
                 abs(X_fit[, indices[[2L]], drop = FALSE]))
  val <- (Z / (max_ref * 1.04))^.78
  val[val > 1] <- 1  # Preserve matrix dimensions for image().

  # In base R, image() puts matrix row 1 at the bottom.  Thus the active
  # observations (last rows of X) appear at the TOP, as in the supplied image.
  graphics::image(
    x = seq(.025, .975, length.out = ncol(val)),
    y = seq(.025, .975, length.out = nrow(val)),
    z = t(val), col = grDevices::gray(seq(1, .08, length.out = 256)),
    xlab = "", ylab = "", axes = FALSE, add = TRUE, useRaster = TRUE
  )
  graphics::rect(.025, .025, .975, .975, border = FRAME_COLOR, lwd = .55)
}

# Project points in the plotting box.  x = row (u) coordinate, y = variable
# (v_j) coordinate, z = surface height in the original units.
# The interpolation between projected bottom and top edges reproduces the
# tilted rectangular boxes of the supplied figure much more closely than the
# old single-affine formula (where both edges had identical slopes).
project_3d <- function(x, y, z, zlim) {
  zz <- (z - zlim[1L]) / diff(zlim)
  bot <- CAMERA$bottom_y + CAMERA$bottom_u * x + CAMERA$bottom_v * y
  top <- CAMERA$top_y + CAMERA$top_u * x + CAMERA$top_v * y
  c(CAMERA$x0 + CAMERA$u_dx * x + CAMERA$v_dx * y,
    (1 - zz) * bot + zz * top)
}

# Compute face color from surface normal, with reference-like dark and light
# faces.  Flat regions stay muted/light; steep portions are more saturated.
# This creates a blue/gold/coral surface rather than a grey wired mesh.
face_color <- function(base, slope_u, slope_v, noisy = FALSE) {
  rgb0 <- as.numeric(grDevices::col2rgb(base)) / 255
  normal <- c(-slope_u, -slope_v, 1)
  normal <- normal / sqrt(sum(normal^2))
  light <- c(-.32, .38, .87)
  light <- light / sqrt(sum(light^2))
  lit <- max(0, sum(normal * light))
  intensity <- .48 + .54 * lit

  flatness <- 1 / (1 + .10 * (slope_u^2 + slope_v^2))
  white_mix <- if (noisy) .025 * flatness else .31 * flatness
  col <- rgb0 * intensity * (1 - white_mix) + white_mix
  col <- pmax(0, pmin(1, col))
  grDevices::rgb(col[1], col[2], col[3])
}

plot_surface <- function(d, j) {
  stopifnot(j %in% c(1L, 3L))

  # Only the ROW dimension is reversed, to make active u appear left/front.
  Z <- d$X[N:1L, indices[[j]], drop = FALSE]
  if (d$name == "Noisy") {
    zlim <- if (j == 1L) c(-1.70, 1.70) else c(-1.45, .85)
  } else {
    zlim <- if (j == 1L) c(-1.22, 1.22) else c(-1.19, .14)
  }

  nr <- if (d$name == "Noisy") 95L else 45L
  nc <- if (d$name == "Noisy") {
    if (j == 1L) 32L else 28L
  } else 25L
  ir <- unique(as.integer(round(seq(1L, nrow(Z), length.out = nr))))
  ic <- unique(as.integer(round(seq(1L, ncol(Z), length.out = nc))))
  z <- Z[ir, ic, drop = FALSE]
  z[z < zlim[1L] + .006] <- zlim[1L] + .006
  z[z > zlim[2L] - .006] <- zlim[2L] - .006
  xs <- seq(0, 1, length.out = nrow(Z))[ir]
  ys <- seq(0, 1, length.out = ncol(Z))[ic]

  edge <- function(a, b, col = FRAME_COLOR, lwd = .48) {
    pa <- project_3d(a[1L], a[2L], a[3L], zlim)
    pb <- project_3d(b[1L], b[2L], b[3L], zlim)
    graphics::segments(pa[1L], pa[2L], pb[1L], pb[2L],
                       col = col, lwd = lwd)
  }

  # Back edges (drawn before the surface).
  for (zz in zlim) {
    for (yy in c(0, 1)) edge(c(0, yy, zz), c(1, yy, zz),
                             col = "#686868", lwd = .35)
    for (xx in c(0, 1)) edge(c(xx, 0, zz), c(xx, 1, zz),
                             col = "#686868", lwd = .35)
  }
  for (xx in c(0, 1)) for (yy in c(0, 1)) {
    edge(c(xx, yy, zlim[1L]), c(xx, yy, zlim[2L]),
         col = "#686868", lwd = .35)
  }

  # Painter's algorithm: farther patches drawn before nearer ones.
  cells <- expand.grid(i = seq_len(nrow(z) - 1L),
                       h = seq_len(ncol(z) - 1L))
  order_faces <- order(xs[cells$i] + ys[cells$h], decreasing = TRUE)
  for (k in order_faces) {
    i <- cells$i[k]
    h <- cells$h[k]
    corners <- rbind(
      project_3d(xs[i],     ys[h],     z[i, h],         zlim),
      project_3d(xs[i + 1], ys[h],     z[i + 1L, h],    zlim),
      project_3d(xs[i + 1], ys[h + 1], z[i + 1L, h + 1L], zlim),
      project_3d(xs[i],     ys[h + 1], z[i, h + 1L],    zlim)
    )
    # Derivatives on normalized axes, for controlled lighting.
    du <- ((z[i + 1L, h] + z[i + 1L, h + 1L]) -
             (z[i, h] + z[i, h + 1L])) /
      (2 * (xs[i + 1L] - xs[i]) * diff(zlim))
    dv <- ((z[i, h + 1L] + z[i + 1L, h + 1L]) -
             (z[i, h] + z[i + 1L, h])) /
      (2 * (ys[h + 1L] - ys[h]) * diff(zlim))
    fill <- face_color(SURFACE_COLORS[[d$name]], du, dv,
                       noisy = identical(d$name, "Noisy"))
    graphics::polygon(corners[, 1L], corners[, 2L],
                      col = fill, border = NA)
  }

  # Thin, color-matched grid lines add detail without the dark wireframe look.
  mesh_col <- grDevices::adjustcolor(SURFACE_COLORS[[d$name]],
                                     alpha.f = if (d$name == "Noisy") .33 else .25)
  for (ii in unique(c(seq(1L, nrow(z), by = 4L), nrow(z)))) {
    if (diff(range(z[ii, ])) < .025 * diff(zlim)) next
    pp <- vapply(seq_along(ys), function(h) {
      project_3d(xs[ii], ys[h], z[ii, h], zlim)
    }, numeric(2L))
    graphics::lines(pp[1L, ], pp[2L, ], col = mesh_col, lwd = .20)
  }
  for (hh in unique(c(seq(1L, ncol(z), by = 3L), ncol(z)))) {
    if (diff(range(z[, hh])) < .025 * diff(zlim)) next
    pp <- vapply(seq_along(xs), function(i) {
      project_3d(xs[i], ys[hh], z[i, hh], zlim)
    }, numeric(2L))
    graphics::lines(pp[1L, ], pp[2L, ], col = mesh_col, lwd = .20)
  }

  # Crisp outlines of the rectangular frame, similar to the supplied figure.
  for (zz in zlim) {
    for (yy in c(0, 1)) edge(c(0, yy, zz), c(1, yy, zz), lwd = .53)
    for (xx in c(0, 1)) edge(c(xx, 0, zz), c(xx, 1, zz), lwd = .53)
  }
  for (xx in c(0, 1)) for (yy in c(0, 1)) {
    edge(c(xx, yy, zlim[1L]), c(xx, yy, zlim[2L]), lwd = .53)
  }

  # Axis labels and the small z-direction arrow seen in the reference.
  graphics::text(.39, .034, "u", cex = .68, srt = -13)
  graphics::text(.936, .172, paste0("v", j), cex = .68, srt = 56)
  graphics::arrows(.068, .31, .068, .68, length = .055,
                   code = 2L, col = FRAME_COLOR, lwd = .55)
}

# ---------------------------------------------------------------------------
# 5. Assemble the 3 x 4 figure and save PDF
# ---------------------------------------------------------------------------
make_figure3 <- function(filename) {
  grDevices::pdf(filename, width = FIG_WIDTH, height = FIG_HEIGHT,
                 paper = "special", pointsize = 10, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  open_panel <- new_panel_factory()

  for (r in seq_along(rows)) {
    d <- rows[[r]]
    gy <- GROUPS[[r]]

    # First column: left vertical u-curve and its title.
    x <- COLS[[1L]]
    open_panel(x[1L], gy[1L], x[2L], gy[4L] - gy[1L])
    plot_u(d)

    # Three variable columns: loading above surface/heatmap.
    for (j in 1:3) {
      x <- COLS[[j + 1L]]
      open_panel(x[1L], gy[3L], x[2L], gy[4L] - gy[3L])
      plot_v(d, j)
      open_panel(x[1L], gy[1L], x[2L], gy[2L] - gy[1L])
      if (j == 2L) plot_discrete(d) else plot_surface(d, j)
    }
  }
  invisible(filename)
}

pdf_file <- file.path(script_directory(), "figure3.pdf")
make_figure3(pdf_file)
message("Saved PDF: ", normalizePath(pdf_file, winslash = "/"))
