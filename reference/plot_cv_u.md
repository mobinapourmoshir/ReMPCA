# Plot CV Error with 1-SE Rule for u

For each principal component (PC), plots cross-validation (CV) error
versus the sparsity tuning parameter (`gamma`) for the u-direction. The
red dashed line indicates the 1-SE rule threshold, and the selected
optimal gamma is highlighted. Optionally, standard error bars (±1 SE)
can be displayed.

## Usage

``` r
plot_cv_u(ReMPCA_obj, show_se = TRUE, ...)
```

## Arguments

- ReMPCA_obj:

  Output list from the ReMPCA routine.

- show_se:

  Logical. If TRUE, adds standard error bars to each point. Default is
  TRUE.

- ...:

  Additional plotting options passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Generates base R plots side-by-side for each PC.
