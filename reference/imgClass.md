# Image Data Class

The `imgClass` class is designed to handle images:

- Supports either a single image (a matrix) or a list of multiple images
  (matrices).

- When a list is provided, all images are vectorized (row-wise) and
  combined into one matrix.

- The resulting object behaves like a `fdClass` if `Smoothing_parameter`
  is not zero, otherwise as a `rdClass`.

- Grid points (optional) are assigned to the columns (i.e., pixel
  locations).

## Usage

``` r
imgClass(image, Smoothing_parameter = 0, Sparsity_parameter = 0, argval = NULL)
```

## Arguments

- image:

  A matrix representing a single image or a list of matrices
  representing multiple images. Each image must be a matrix of the same
  dimension.

- Smoothing_parameter:

  - A numeric scalar or vector controlling the level of smoothing (as in
    `fdClass`).

  - Set to 0 for no smoothing (produces an `rdClass` object).

  - If NULL, smoothing is tuned from a default sequence.

- Sparsity_parameter:

  - A fixed number or vector controlling the level of sparsity on pixel
    columns.

  - If NULL, it is tuned automatically based on image size.

- argval:

  Optional numeric vector for grid points over the pixels (columns).

## Value

An object of class `fdClass` or `rdClass`, depending on the
`Smoothing_parameter`.

## Examples

``` r
img1 <- matrix(rnorm(64), 8, 8)
img2 <- matrix(rnorm(64), 8, 8)

# Multiple images
rd_img <- imgClass(list(img1, img2), Smoothing_parameter = 0, Sparsity_parameter = NULL)
```
