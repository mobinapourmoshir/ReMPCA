# Plot Principal Component Functions

This function visualizes the estimated principal component functions (v)
from a ReMPCA object. Each component is plotted across its variables
using either a line plot (for functional variables) or a dot plot (for
regular variables).

## Usage

``` r
plot_pc_functions(ReMPCA_obj, ...)
```

## Arguments

- ReMPCA_obj:

  A ReMPCA object containing:

  `PCFunctions`

  :   A list of length equal to number of PCs. Each element is a list of
      PC functions for each variable.

  `variable_types`

  :   A character vector indicating type of each variable: either `"hd"`
      (functional) or `"rd"` (regular).

- ...:

  Additional plotting options.

## Value

No return value. This function is called for its side effect of
plotting.

## Details

The function arranges plots in a matrix layout with rows corresponding
to variables and columns to principal components. A light gray
horizontal line at 0 is added for reference unless the minimum value in
the plot is ≥ 5.
