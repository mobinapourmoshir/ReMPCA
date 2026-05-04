# Coerce an Object to imgClass

Converts an object of class `rdClass` or `fdClass` to an `imgClass`,
while preserving or overriding associated attributes such as smoothing,
sparsity, and grid points.

## Usage

``` r
as.imgClass(
  x,
  Sparsity_parameter = NULL,
  Smoothing_parameter = NULL,
  argval = NULL
)
```

## Arguments

- x:

  An object of class `rdClass` or `fdClass`.

- Sparsity_parameter:

  Optional. A numeric vector of non-negative integers representing
  sparsity levels to apply. If `NULL`, the parameter is inherited from
  `x`.

- Smoothing_parameter:

  Optional. A numeric value or vector representing smoothing parameters.
  If `NULL`, the parameter is inherited from `x`.

- argval:

  Optional. A numeric vector of grid points (argvals) for functional
  representation. If `NULL`, the grid is inherited from `x`.

## Value

An object of class `imgClass`, which also inherits from `fdClass` or
`rdClass`, depending on the smoothing parameter.

## Details

This coercion is helpful when an object originally treated as regular or
functional data (via `rdClass` or `fdClass`) should instead be
interpreted and processed as image data.

## Examples

``` r
mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
fd_obj <- fdClass(mat, Smoothing_parameter = 0.1)
img_obj <- as.imgClass(fd_obj)
print(class(img_obj))              # "imgClass" "fdClass"
#> [1] "imgClass" "fdClass" 
print(attr(img_obj, "Smoothing_parameter"))
#> [1] 0.1
```
