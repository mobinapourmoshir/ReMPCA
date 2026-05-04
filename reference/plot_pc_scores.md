# Plot Principal Component Scores from ReMPCA Object

This function visualizes the principal component scores from a ReMPCA
object. Each principal component (i.e., each column in `PCScores`) is
plotted in a separate dot plot, displaying the score values across
observations.

## Usage

``` r
plot_pc_scores(ReMPCA_obj, ...)
```

## Arguments

- ReMPCA_obj:

  A list-like ReMPCA object, expected to contain a component named
  `PCScores` which is a matrix or data frame where each column
  corresponds to the scores of one principal component.

- ...:

  Additional plotting options.

## Value

No return value. The function produces a series of dot plots, one for
each principal component.
