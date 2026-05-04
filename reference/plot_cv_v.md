# Plot Cross-Validation Scores for ReMPCA v-direction

This function creates a grid of plots showing cross-validation (CV)
scores for each principal component (PC) and functional variable based
on different values of the gamma penalty parameter in the v-direction.
The red dashed line indicates the 1-SE rule threshold, and the selected
optimal gamma is highlighted.

## Usage

``` r
plot_cv_v(ReMPCA_obj, show_se = TRUE, ...)
```

## Arguments

- ReMPCA_obj:

  A list object returned by ReMPCA containing `CVResultsV` and
  `OptimalGammaV`.

- show_se:

  Logical. If TRUE, adds standard error bars (±1 SE) to each point.
  Default is TRUE.

- ...:

  Additional arguments passed to the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) function.
